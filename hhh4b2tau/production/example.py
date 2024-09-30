# coding: utf-8

"""
Column production methods related to higher-level features.
"""


from columnflow.production import Producer, producer
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.cms.muon import muon_weights
from columnflow.production.cms.electron import electron_weights
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

# my own variables
from hhh4b2tau.production.newvariables import jet_angle_difference
from hhh4b2tau.production.newvariables import hhh_decay_invariant_mass
from hhh4b2tau.production.newvariables import final_state_variables

np = maybe_import("numpy")
ak = maybe_import("awkward")


@producer(
    uses={
        # nano columns
        "Jet.pt",
    },
    produces={
        # new columns
        "ht", "n_jet",
    },
)
def features(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    events = set_ak_column(events, "ht", ak.sum(events.Jet.pt, axis=1))
    events = set_ak_column(events, "n_jet", ak.num(events.Jet.pt, axis=1), value_type=np.int32)

    return events


@producer(
    uses={
        mc_weight, category_ids,
        # nano columns
        "Jet.pt",
    },
    produces={
        mc_weight, category_ids,
        # new columns
        "cutflow.jet1_pt",
    },
)
def cutflow_features(
    self: Producer,
    events: ak.Array,
    object_masks: dict[str, dict[str, ak.Array]],
    **kwargs,
) -> ak.Array:
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)

    # apply object masks and create new collections
    reduced_events = events
    if object_masks:
        reduced_events = create_collections_from_masks(events, object_masks)

    # create category ids per event and add categories back to the
    events = self[category_ids](reduced_events, target_events=events, **kwargs)

    # add cutflow columns
    events = set_ak_column(
        events,
        "cutflow.jet1_pt",
        Route("Jet.pt[:,0]").apply(events, EMPTY_FLOAT),
    )

    return events


@producer(
    uses={
        features, category_ids, normalization_weights, muon_weights, deterministic_seeds,
    },
    produces={
        features, category_ids, normalization_weights, muon_weights, deterministic_seeds,
    },
)
def example(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    # features
    events = self[features](events, **kwargs)

    # category ids
    events = self[category_ids](events, **kwargs)

    # deterministic seeds
    events = self[deterministic_seeds](events, **kwargs)

    # mc-only weights
    if self.dataset_inst.is_mc:
        # normalization weights
        events = self[normalization_weights](events, **kwargs)

        # muon weights
        events = self[muon_weights](events, **kwargs)

    return events

@producer(
    uses={
        features, category_ids, normalization_weights, deterministic_seeds, jet_angle_difference, hhh_decay_invariant_mass, final_state_variables
    },
    produces={
        features, category_ids, normalization_weights, deterministic_seeds, jet_angle_difference, hhh_decay_invariant_mass, final_state_variables
    },
)
def empty(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events = self[normalization_weights](events, **kwargs)
    # features
    events = self[features](events, **kwargs)

    # category ids
    events = self[category_ids](events, **kwargs)

    # deterministic seeds
    events = self[deterministic_seeds](events, **kwargs)

    # adding new variables
    events = self[jet_angle_difference](events, **kwargs)

    if (self.dataset_inst.is_mc and
        any(self.dataset_inst.name.lower().startswith(x)
            for x in ("hhh",))
    ):
        events = self[hhh_decay_invariant_mass](events, **kwargs)
        
    if (self.dataset_inst.is_mc and
        any(self.dataset_inst.name.lower().startswith(x)
            for x in ("tth_hbb_powheg",))
    ):
        events = self[final_state_variables](events, **kwargs)

    return events
