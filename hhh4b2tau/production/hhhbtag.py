# coding: utf-8

"""
Producers for the HHHBtag score.
See https://github.com/hh-italian-group/HHbtag.
"""

import law

from columnflow.production import Producer, producer
from columnflow.util import maybe_import, dev_sandbox, InsertableDict
from columnflow.columnar_util import EMPTY_FLOAT, layout_ak_array

from hhh4b2tau.util import IF_RUN_2

np = maybe_import("numpy")
ak = maybe_import("awkward")

logger = law.logger.get_logger(__name__)


@producer(
    uses={
        # custom columns created upstream, probably by a selector
        "channel_id",
        # nano columns
        "event",
        "Jet.{pt,eta,phi,mass,jetId,btagDeepFlavB}", 
        IF_RUN_2("Jet.puId"),
        "MET.{pt,phi}",
    },
    sandbox=dev_sandbox("bash::$HHH4B2TAU_BASE/sandboxes/venv_columnar_tf.sh"),
)
def hhhbtag(
    self: Producer,
    events: ak.Array,
    jet_mask: ak.Array,
    lepton_pair: ak.Array,
    **kwargs,
) -> ak.Array:
    """
    Returns the HHHBtag score per passed jet.
    """
    # get a mask of events where there are at least two tau candidates and at least two jets
    # and only get the scores for jets in these events
    event_mask = (
        (ak.num(lepton_pair, axis=1) >= 2) &
        (ak.sum(jet_mask, axis=1) >= 2)
    )

    # prepare objects
    n_jets_max = 10
    jets = events.Jet[jet_mask][event_mask][..., :n_jets_max]
    leps = lepton_pair[event_mask][..., [0, 1]]
    htt = leps[..., 0] + leps[..., 1]
    met = events[event_mask].MET
    jet_shape = abs(jets.pt) >= 0
    n_jets_capped = ak.num(jets, axis=1)

    # get input features
    input_features = [
        jet_shape * 1,
        jets.pt,
        jets.eta,
        jets.mass / jets.pt,
        jets.energy / jets.pt,
        abs(jets.eta - htt.eta),
        jets.btagDeepFlavB,
        jets.delta_phi(htt),
        jet_shape * (self.config_inst.campaign.x.year),
        jet_shape * (events[event_mask].channel_id - 1),
        jet_shape * htt.pt,
        jet_shape * htt.eta,
        jet_shape * htt.delta_phi(met),
        jet_shape * (met.pt / htt.pt),
        jet_shape * ak.sum(leps.pt, axis=1),
    ]

    # helper to split events, cast to float32, concatenate across new axis,
    # then pad with zeros for up to n_jets_max jets
    def split(where):
        features = ak.concatenate(
            [
                ak.values_astype(f[where][..., None], np.float32)
                for f in input_features
            ],
            axis=2,
        )
        # fill
        features = ak.fill_none(
            ak.pad_none(features, n_jets_max, axis=1),
            np.zeros(len(input_features), dtype=np.float32),
            axis=1,
        )
        # fix the dimension of the last axis to the known number of input features
        features = features[..., list(range(len(input_features)))]
        return ak.to_numpy(features)

    # reserve an output score array
    scores = np.ones((ak.sum(event_mask), n_jets_max), dtype=np.float32) * EMPTY_FLOAT

    # fill even and odd events if there are any
    even_mask = ak.to_numpy((events[event_mask].event % 2) == 0)
    if ak.sum(even_mask):
        input_features_even = split(even_mask)
        scores_even = self.hhhbtag_model_even(input_features_even)[0].numpy()
        scores[even_mask] = scores_even
    if ak.sum(~even_mask):
        input_features_odd = split(~even_mask)
        scores_odd = self.hhhbtag_model_odd(input_features_odd)[0].numpy()
        scores[~even_mask] = scores_odd

    # remove the scores of padded jets
    where = ak.from_regular(ak.local_index(scores) < n_jets_capped[..., None], axis=1)
    scores = ak.from_regular(scores, axis=1)[where]

    # add scores to events that had more than n_jets_max selected jets
    # (use zero here as this is also what the hhhbtag model does for missing jets)
    layout_ext = events.Jet.pt[jet_mask][event_mask][..., n_jets_max:]
    # when there are no selected events, we can reuse layout_ext and consider it to be scores_ext
    if len(layout_ext) == 0:
        scores_ext = layout_ext
    else:
        scores_ext = layout_ak_array(np.zeros(len(ak.flatten(layout_ext)), dtype=np.int32), layout_ext)
    scores = ak.concatenate([scores, scores_ext], axis=1)

    # prevent Memory Corruption Error
    jet_mask = ak.fill_none(jet_mask, False, axis=-1)

    # insert scores into an array with same shape as input jets (without jet_mask and event_mask)
    all_scores = ak.fill_none(ak.full_like(events.Jet.pt, EMPTY_FLOAT, dtype=np.float32), EMPTY_FLOAT, axis=-1)
    np.asarray(ak.flatten(all_scores))[ak.flatten(jet_mask & event_mask, axis=1)] = np.asarray(ak.flatten(scores))

    return all_scores


@hhhbtag.requires
def hhhbtag_requires(self: Producer, reqs: dict) -> None:
    """
    Add the external files bundle to requirements.
    """
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)


@hhhbtag.setup
def hhhbtag_setup(self: Producer, reqs: dict, inputs: dict, reader_targets: InsertableDict) -> None:
    """
    Sets up the two HHHBtag TF models.
    """
    tf = maybe_import("tensorflow")

    # unpack the external files bundle, create a subdiretory and unpack the hhbtag repo in it
    bundle = reqs["external_files"]
    arc = bundle.files.hh_btag_repo
    repo_dir = bundle.files_dir.child("hh_btag_repo", type="d")
    arc.load(repo_dir, formatter="tar")
    repo_dir = repo_dir.child(repo_dir.listdir(pattern="HHbtag-*")[0])

    # get the version of the external file
    self.hhhbtag_version = self.config_inst.x.external_files["hh_btag_repo"][1]

    # define the model path
    model_path = f"models/HHbtag_{self.hhhbtag_version}_par"

    # save both models (even and odd event numbers)
    with self.task.publish_step("loading hhhbtag models ..."):
        self.hhhbtag_model_even = tf.saved_model.load(repo_dir.child(f"{model_path}_0").path)
        self.hhhbtag_model_odd = tf.saved_model.load(repo_dir.child(f"{model_path}_1").path)