# coding: utf-8

"""
Test model definition. (unfinished)
"""

from __future__ import annotations

from typing import Any, Sequence
from collections import defaultdict

import law
import order as od

from columnflow.ml import MLModel
from columnflow.util import maybe_import, dev_sandbox
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, remove_ak_column


ak = maybe_import("awkward")
np = maybe_import("numpy")
tf = maybe_import("tensorflow")
nprec = maybe_import("numpy.lib.recfunctions")

law.contrib.load("tensorflow")

logger = law.logger.get_logger(__name__)

hhh_couplings = ((0, 0), (1, 0), (4, 9))

class HHHModel(MLModel):

    dataset_names = [
        "tth_hbb_powheg",
        "tth_hnonbb_powheg",
        "hh_ggf_hbb_htt_kl1_kt1_powheg",
        "tt_sl_powheg",
        "tt_dl_powheg",
        "tt_fh_powheg",
    ] + ["hhh_4b2tau_c3{c3}_d4{d4}_amcatnlo".format(
            c3=str(c3).replace("-", "m").replace(".", "p"),
            d4=str(d4).replace("-", "m").replace(".", "p"),
            ) for c3,d4 in hhh_couplings]


    def __init__(
        self,
        *args,
        folds: int | None = None,
        model_name: str | None = None,
        activation_func: str | None = None,
        batch_norm: bool | None = None,
        input_features: set[str] | None = None,
        nodes: list | None = None,
        L2: bool | None = None,
        validation_split: float | None = None,
        epochs: int = 10,
        **kwargs,
    ):
        """
        Parameters that need to be set by derived model:
        folds, layers, learningrate, batchsize, epochs, eqweight, dropout,
        processes, ml_process_weights, dataset_names, input_features, store_name,
        """
        single_config = True  # noqa

        super().__init__(*args, **kwargs)

        # class- to instance-level attributes
        # (before being set, self.folds refers to a class-level attribute)
        self.folds = folds or self.folds
        self.model_name = model_name or self.model_name
        self.input_features = input_features or self.input_features
        self.batch_norm = batch_norm or self.batch_norm
        self.activation_func = activation_func or self.activation_func
        self.L2 = L2 or self.L2
        self.nodes = nodes or self.nodes
        self.validation_split = validation_split or self.validation_split
        self.epochs = epochs

        # store parameters of interest in the ml_model_inst, e.g. via the parameters attribute
        self.parameters = self.parameters | {
            "batchsize": int(self.parameters.get("batchsize", 512)),
            "layers": tuple(int(layer) for layer in self.parameters.get("layers", self.nodes)),
            "epochs": int(self.parameters.get("epochs", self.epochs)),
        }

    def setup(self):
        self.dataset_insts = {
            dataset_name: self.config_inst.get_dataset(dataset_name)
            for dataset_name in self.dataset_names
        }
        self.process_insts = {}

        for id, dataset in enumerate(self.dataset_insts.values()):
            proc = dataset.processes.get_first()
            proc.x.ml_id = id
            proc.x.process_weight = 1  # TODO: set process weights
            self.process_insts[dataset.name] = proc

        # dynamically add variables for the quantities produced by this model
        for c3,d4 in hhh_couplings:
            if f"{self.cls_name}.score_kl{c3+1}_kc_{d4+1}" not in self.config_inst.variables:
                self.config_inst.add_variable(
                    name=f"{self.cls_name}.score_kl{c3+1}_kc_{d4+1}",
                    null_value=-1,
                    binning=(10, 0, 1),
                    x_title=f"HHH Model output ($\kappa_{{\lambda}}$={c3+1}, $\kappa_{{\lambda}}$={d4+1})",
                )

    def sandbox(self, task: law.Task) -> str:
        return dev_sandbox("bash::$HHH4B2TAU_BASE/sandboxes/example.sh")

    def datasets(self, config_inst: od.Config) -> set[od.Dataset]:
        return {config_inst.get_dataset(dataset_name) for dataset_name in self.dataset_names}

    def uses(self, config_inst: od.Config) -> set[Route | str]:
        columns = set(self.input_features) | {"normalization_weight", "process_ids"}
        return columns

    def training_producers(self, config_inst: od.Config, requested_producers: Sequence[str]) -> list[str]:
        producers = {"default", "deterministic_seeds", "res_pdnn", "hh_mass"} | set(requested_producers)
        # fix MLTraining Phase Space
        return sorted(list(producers))

    def produces(self, config_inst: od.Config) -> set[Route | str]:
        # mark columns that you don't want to be filtered out
        # preserve the networks prediction of a specific feature for each fold
        # cls_name would be the name of your model
        ml_predictions = {f"{self.cls_name}.fold{fold}.{feature}"
            for fold in range(self.folds)
            for feature in self.target_columns}

        # save indices used to create the folds
        util_columns = {f"{self.cls_name}.fold_indices"}
        additional_columns = {"process_ids"}
        # combine all columns to a unique set
        preserved_columns = ml_predictions | util_columns | additional_columns
        return preserved_columns

    def output(self, task: law.Task) -> law.FileSystemDirectoryTarget:
        # needs to be given via config
        max_folds = self.folds
        current_fold = task.fold

        # create directory at task.target, if it does not exist
        target = task.target(f"mlmodel_f{current_fold}of{max_folds}.keras", dir=False)
        return target

    def build_field_names(self, dtype: np.typing.DTypeLike) -> list[str]:
        fields = []
        for (field, typ) in dtype.descr:
            if isinstance(typ, list):
                fields.extend([f"{field}.{subfield[0]}" for subfield in typ])
            else:
                fields.append(field)
        return fields

    def merge_weight_and_class_weights(self, data: dict, class_weight: dict[int, float]) -> np.ndarray:
        indecies = np.argmax(data["target"], axis=-1)
        weight = np.zeros_like(data["weight"])

        for i, w in class_weight.items():
            weight[indecies == i] = w * data["weight"][indecies == i]

        return weight

    def open_model(self, target: law.FileSystemDirectoryTarget):
        return target.load(formatter="tf_keras_model")

    def build_model(self) -> Any:
        # helper function to handle model building
        x_inp = tf.keras.Input(shape=(len(self.input_features),))
        x = x_inp
        for node in self.nodes:
            if self.batch_norm:
                x = tf.keras.layers.BatchNormalization()(x)
            if self.L2:
                x = tf.keras.layers.Dense(
                    node,
                    activation=self.activation_func,
                    kernel_regularizer=tf.keras.regularizers.l2(0.01),
                )(x)
            else:
                x = tf.keras.layers.Dense(node, activation=self.activation_func)(x)
        y = tf.keras.layers.Dense(4, activation="softmax")(x)
        model = tf.keras.Model(inputs=x_inp, outputs=y)
        return model

    def prepare_input(self, input: ak.Array) -> dict[str, tf.Tensor]:
        weight_sum: dict[str, float] = {}
        dataset_proc_idx: dict[str, int] = {}
        training: defaultdict[str, list] = defaultdict(list)
        valid: defaultdict[str, list] = defaultdict(list)
        fields: list = None

        for dataset, files in input["events"][self.config_inst.name].items():
            dataset_inst = self.dataset_insts[dataset]

            if len(dataset_inst.processes) != 1:
                raise Exception("only 1 process inst is expected for each dataset")

            # calculate the sum of weights for each dataset
            weight_sum[dataset] = sum(ak.sum(file["mlevents"].load().normalization_weight) for file in files)

            for proc in self.process_insts.values():
                leaf_procs = [p for p, _, _ in self.config_inst.get_process(proc).walk_processes(include_self=True)]
                if dataset_inst.processes.get_first() in leaf_procs:
                    logger.info(f"the dataset *{dataset}* is used for training the *{proc.name}* output node")
                    dataset_proc_idx[dataset] = proc.x.ml_id
                    continue

            if dataset_proc_idx.get(dataset, -1) == -1:
                raise Exception(f"dataset {dataset} is not matched to any of the given processes")

            for i, inp in enumerate(files):
                events = inp["mlevents"].load()
                weights = ak.to_numpy(events.normalization_weight)

                events = remove_ak_column(events, "normalization_weight")
                events = events.to_numpy()

                # create field names and check if they are matching
                if not fields:
                    fields = self.build_field_names(events.dtype)
                else:
                    if fields != self.build_field_names(events.dtype):
                        raise Exception("fields are not matching")

                events = nprec.structured_to_unstructured(events)

                # set EMPTY_FLOAT to -10
                events[events == EMPTY_FLOAT] = -10

                # check for infinite values in weights
                if np.any(~np.isfinite(weights)):
                    raise Exception(f"Infinite values found in weights from dataset {dataset}")

                # create target array
                target = np.zeros((len(events), len(self.process_insts)))
                target[:, dataset_proc_idx[dataset]] = 1

                # split into training and validation set
                if self.validation_split:
                    split = int(len(events) * self.validation_split)
                    choice = np.random.choice(range(events.shape[0]), size=(split,), replace=False)
                    ind = np.zeros(events.shape[0], dtype=bool)
                    ind[choice] = True
                    valid["events"].append(events[ind])
                    events = events[~ind]
                    valid["target"].append(target[ind])
                    target = target[~ind]
                    valid["weight"].append(weights[ind])
                    weights = weights[~ind]

                    logger.info(f"file {i} of dataset *{dataset}* split into {len(events)} training and {split}"
                                " validation events")

                training["events"].append(events)
                training["target"].append(target)
                training["weight"].append(weights)

        mean_weight = np.mean(list(weight_sum.values()))
        class_weight = {dataset_proc_idx[d]: mean_weight / w for d, w in weight_sum.items()}

        # concatenate all events and targets
        training["events"] = np.concatenate(training["events"])
        training["target"] = np.concatenate(training["target"])
        training["weight"] = np.concatenate(training["weight"])
        training["weight"] = self.merge_weight_and_class_weights(training, class_weight)
        # helper function to get tensor tuple
        get_tensor_tuple = lambda x: tuple(x[field] for field in ["events", "target", "weight"])

        # create tf tensors
        train_tensor = tf.data.Dataset.from_tensor_slices(get_tensor_tuple(training))
        train_tensor = train_tensor.shuffle(buffer_size=train_tensor.cardinality())
        train_tensor = train_tensor.batch(self.parameters["batchsize"])
        train_tensor = train_tensor.prefetch(tf.data.experimental.AUTOTUNE)

        # create an output for the fit function of tf
        result = {"x": train_tensor}  # ,"class_weight": class_weight} included in sample weights

        if self.validation_split:
            valid["events"] = np.concatenate(valid["events"])
            valid["target"] = np.concatenate(valid["target"])
            valid["weight"] = np.concatenate(valid["weight"])
            valid["weight"] = self.merge_weight_and_class_weights(valid, class_weight)

            valid_tensor = tf.data.Dataset.from_tensor_slices(get_tensor_tuple(valid))
            valid_tensor = valid_tensor.shuffle(buffer_size=valid_tensor.cardinality())
            valid_tensor = valid_tensor.batch(self.parameters["batchsize"])
            valid_tensor = valid_tensor.prefetch(tf.data.experimental.AUTOTUNE)
            result["validation_data"] = valid_tensor

        return result, fields

    def train(
        self,
        task: law.Task,
        input: dict[str, list[law.FileSystemFileTarget]],
        output: law.FileSystemDirectoryTarget,
    ) -> None:
        physical_devices = tf.config.list_physical_devices("GPU")
        try:
            tf.config.experimental.set_memory_growth(physical_devices[0], True)
        except:
            # Invalid device or cannot modify virtual devices once initialized.
            pass
        # use helper functions to define model, open input parquet files and prepare events
        # init a model structure
        model = self.build_model()

        # get data tensors
        data_config, feautres = self.prepare_input(input)
        model_config = {
            "epochs": self.parameters["epochs"],
            "steps_per_epoch": 100,
            "validation_freq": 5,
        }

        # setup everything needed for training
        optimizer = tf.keras.optimizers.SGD()
        # optimizer = tf.keras.optimizers.Adam(
        #     learning_rate=self.learningrate, beta_1=0.9, beta_2=0.999,
        #     epsilon=1e-6, amsgrad=False,
        # )

        model.compile(
            optimizer,
            loss="categorical_crossentropy",
            steps_per_execution=10,
        )

        # print model summary
        model.summary()

        # train, throw model_history away
        _ = model.fit(
            **data_config,
            **model_config,
        )

        # save your model and everything you want to keep
        output.dump(model, formatter="tf_keras_model")

    def evaluate(
        self,
        task: law.Task,
        events: ak.Array,
        models: list[Any],
        fold_indices: ak.Array,
        events_used_in_training: bool = False,
    ) -> ak.Array:
        # prepare ml_models input features, target is not important
        inputs_tensor, _ = self.prepare_events(events)

        # do evaluation on all data
        # one can get test_set of fold using: events[fold_indices == fold]
        for fold, model in enumerate(models):
            # convert tf.tensor -> np.array -> ak.array
            # to_regular is necessary to make array contigous (otherwise error)
            prediction = ak.to_regular(model(inputs_tensor).numpy())

            # update events with predictions, sliced by feature, and fold_indices for identification purpose
            for index_feature, target_feature in enumerate(self.target_features):
                events = set_ak_column(
                    events,
                    f"{self.cls_name}.fold{fold}.{target_feature}",
                    prediction[:, index_feature],
                )

        events = set_ak_column(
            events,
            f"{self.cls_name}.{fold_indices}",
            fold_indices,
        )

        return events


# usable derivations
hhh_model = HHHModel.derive("hhh_model", cls_dict={
    "folds": 3,
    "model_name": "hhh_model",
    "activation_func": "relu",
    "batch_norm": False,
    "nodes": [16, 16, 16],
    "validation_split": 0.2,
    "L2": True,
    "epochs": 50,
    "input_features": {
        f"res_pdnn_s0_m500_{node}" for node in ["hhh", "dy", "tt"]
    } | {
        f"hhh.{var}" for var in ["pt", "eta", "mass"]
    },
    "target_features": {
        f"hhh_model.score_kl{c3+1}_kc_{d4+1}" for c3,d4 in hhh_couplings
    },
})