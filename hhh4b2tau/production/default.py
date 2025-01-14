# coding: utf-8

"""
Wrappers for some default sets of producers.
"""

from columnflow.production import Producer, producer
from columnflow.production.normalization import stitched_normalization_weights
from columnflow.production.categories import category_ids
from columnflow.production.cms.electron import electron_weights
from columnflow.production.cms.muon import muon_weights
from columnflow.util import maybe_import

from hhh4b2tau.production.features import features
from hhh4b2tau.production.weights import (
    normalized_pu_weight, normalized_pdf_weight, normalized_murmuf_weight,
)
# from hhh4b2tau.production.btag import normalized_btag_weights
from hhh4b2tau.production.btag import normalized_btag_weights_deepjet, normalized_btag_weights_pnet
from hhh4b2tau.production.tau import tau_weights, trigger_weights
from hhh4b2tau.production.newvariables import jet_angle_difference
from hhh4b2tau.production.newvariables import dectector_variables
from hhh4b2tau.util import IF_DATASET_HAS_LHE_WEIGHTS, IF_RUN_3


ak = maybe_import("awkward")


@producer(
    uses={
        category_ids, features, stitched_normalization_weights, normalized_pu_weight,
        # tau_weights, trigger_weights,
        normalized_btag_weights_deepjet, IF_RUN_3(normalized_btag_weights_pnet),
        electron_weights, muon_weights, jet_angle_difference, dectector_variables, 
        IF_DATASET_HAS_LHE_WEIGHTS(normalized_pdf_weight, normalized_murmuf_weight),
    },
    produces={
        category_ids, features, stitched_normalization_weights, normalized_pu_weight,
        # tau_weights, trigger_weights,
        normalized_btag_weights_deepjet, IF_RUN_3(normalized_btag_weights_pnet),
        electron_weights, muon_weights, jet_angle_difference, dectector_variables, 
        IF_DATASET_HAS_LHE_WEIGHTS(normalized_pdf_weight, normalized_murmuf_weight),
    },
    produce_weights=True,
)
def default(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    # from IPython import embed; embed(header="default producer")
    # category ids
    events = self[category_ids](events, **kwargs)

    # features
    events = self[features](events, **kwargs)

    # mc-only weights
    if self.dataset_inst.is_mc:
        # normalization weights
        events = self[stitched_normalization_weights](events, **kwargs)

        # normalized pdf weight
        if self.has_dep(normalized_pdf_weight):
            events = self[normalized_pdf_weight](events, **kwargs)

        # normalized renorm./fact. weight
        if self.has_dep(normalized_murmuf_weight):
            events = self[normalized_murmuf_weight](events, **kwargs)

        # normalized pu weights
        events = self[normalized_pu_weight](events, **kwargs)

        # btag weights
        events = self[normalized_btag_weights_deepjet](events, **kwargs)
        if self.has_dep(normalized_btag_weights_pnet):
            events = self[normalized_btag_weights_pnet](events, **kwargs)

        # tau weights
        if self.has_dep(tau_weights):
            events = self[tau_weights](events, **kwargs)

        # electron weights ## work in progress
        # if self.has_dep(electron_weights):
        #     events = self[electron_weights](events, **kwargs)

        # muon weights
        if self.has_dep(muon_weights):
            events = self[muon_weights](events, **kwargs)

        # trigger weights
        if self.has_dep(trigger_weights):
            events = self[trigger_weights](events, **kwargs)
        
    events = self[jet_angle_difference](events, **kwargs)
    events = self[dectector_variables](events, **kwargs)

    return events

@default.init
def default_init(self: Producer) -> None:
    if self.produce_weights:
        weight_producers = {tau_weights, electron_weights, muon_weights, trigger_weights}
        self.uses |= weight_producers
        self.produces |= weight_producers


empty = default.derive("empty", cls_dict={"produce_weights": False})
