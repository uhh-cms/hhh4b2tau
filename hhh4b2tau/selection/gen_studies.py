# coding: utf-8

"""
Exemplary selection methods.
"""

from collections import defaultdict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.columnar_util import sorted_indices_from_mask
from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.util import maybe_import

from hhh4b2tau.production.example import cutflow_features
from hhh4b2tau.production.gen_higgs_decay_products import gen_higgs_decay_products
from hhh4b2tau.production.gen_higgs_decay_products import gen_tth_decay_products
from hhh4b2tau.production.gen_higgs_decay_products import gen_Hadron_products
np = maybe_import("numpy")
ak = maybe_import("awkward")


#
# other unexposed selectors
# (not selectable from the command line but used by other, exposed selectors)
#

#
# exposed selectors
# (those that can be invoked from the command line)
#

@selector(
    uses={
        # selectors / producers called within _this_ selector
        mc_weight, cutflow_features, process_ids, increment_stats,
        gen_higgs_decay_products, gen_tth_decay_products, gen_Hadron_products,
        
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        mc_weight, cutflow_features, process_ids, 
        gen_higgs_decay_products, gen_tth_decay_products, gen_Hadron_products,
    },
    exposed=True,
)
def gen_studies(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # get higgs decay products

    if (self.dataset_inst.is_mc and
        any(self.dataset_inst.name.lower().startswith(x)
            for x in ("hhh",))
    ):
        events = self[gen_higgs_decay_products](events, **kwargs)
    
    if (self.dataset_inst.is_mc and
        any(self.dataset_inst.name.lower().startswith(x)
            for x in ("tth_hbb_powheg",))
    ):
        events = self[gen_tth_decay_products](events, **kwargs)

    events = self[gen_Hadron_products](events, **kwargs)

    # select events with at least 4 gen b jets
    n_gen_b_jet = ak.num(events.gen_b_jet)

    results.steps["one_b_jets"] = n_gen_b_jet >= 1
    results.steps["two_b_jets"] = n_gen_b_jet >= 2
    results.steps["three_b_jets"] = n_gen_b_jet >= 3
    # results.steps["four_b_jets"] = n_gen_b_jet >= 4

    # select events with at least 2 GenVisTau
    n_GenVisTau =ak.num(events.GenVisTau)
    
    results.steps["one_tau"] = n_GenVisTau >= 1
    results.steps["two_tau"] = n_GenVisTau >= 2
        
    # # get tau decay products

    # if self.dataset_inst.is_mc and self.dataset_inst.name.startswith("tau"):
    #     events = self[gen_higgs_decay_products](events, **kwargs)
    
    results.event = ak.ones_like(events.event, dtype=bool)

    # create process ids
    events = self[process_ids](events, **kwargs)
    
    # add the mc weight
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)

    # add cutflow features, passing per-object masks
    events = self[cutflow_features](events, results.objects, **kwargs)

    # increment stats
    weight_map = {
        "num_events": Ellipsis,
        "num_events_selected": results.event,
    }
    group_map = {}
    if self.dataset_inst.is_mc:
        weight_map = {
            **weight_map,
            # mc weight for all events
            "sum_mc_weight": (events.mc_weight, Ellipsis),
            "sum_mc_weight_selected": (events.mc_weight, results.event),
        }
        group_map = {
            # per process
            "process": {
                "values": events.process_id,
                "mask_fn": (lambda v: events.process_id == v),
            },
            # per jet multiplicity
            "njet": {
                "values": results.x.n_jets,
                "mask_fn": (lambda v: results.x.n_jets == v),
            },
        }
    events, results = self[increment_stats](
        events,
        results,
        stats,
        weight_map=weight_map,
        group_map=group_map,
        **kwargs,
    )

    # from IPython import embed; embed(header="end selector")

    return events, results

