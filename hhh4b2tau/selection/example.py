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

np = maybe_import("numpy")
ak = maybe_import("awkward")


#
# other unexposed selectors
# (not selectable from the command line but used by other, exposed selectors)
#


@selector(
    uses={f"Muon.{var}" 
          for var in ("pt", "eta", )},
)
def muon_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    # example muon selection: exactly one muon
    muon_mask = (events.Muon.pt >= 20.0) & (abs(events.Muon.eta) < 2.1)
    muon_sel = ak.sum(muon_mask, axis=1) == 1

    # build and return selection results
    # "objects" maps source columns to new columns and selections to be applied on the old columns
    # to create them, e.g. {"Muon": {"MySelectedMuon": indices_applied_to_Muon}}
    return events, SelectionResult(
        steps={
            "muon": muon_sel,
        },
        objects={
            "Muon": {
                "Muon": muon_mask,
            },
        },
    )

@selector(
    uses={f"Electron.{var}" 
          for var in ("pt", "eta", )},
)
def electron_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    # example electron selection: exactly one electron
    electron_mask = (events.Electron.pt >= 20.0) & (abs(events.Electron.eta) < 2.1)
    electron_sel = ak.sum(electron_mask, axis=1) == 1

    # build and return selection results
    # "objects" maps source columns to new columns and selections to be applied on the old columns
    # to create them, e.g. {"Electron": {"MySelectedElectron": indices_applied_to_Electron}}
    return events, SelectionResult(
        steps={
            "electron": electron_sel,
        },
        objects={
            "Electron": {
                "Electron": electron_mask,
            },
        },
    )

@selector(
    uses={f"Jet.{var}" 
          for var in ("pt", "eta", )},
)
def jet_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    # example jet selection: at least one jet
    jet_mask = (events.Jet.pt >= 25.0) & (abs(events.Jet.eta) < 2.4)
    jet_sel = ak.sum(jet_mask, axis=1) >= 1

    # build and return selection results
    # "objects" maps source columns to new columns and selections to be applied on the old columns
    # to create them, e.g. {"Jet": {"MyCustomJetCollection": indices_applied_to_Jet}}
    return events, SelectionResult(
        steps={
            "jet": jet_sel,
        },
        objects={
            "Jet": {
                "Jet": sorted_indices_from_mask(jet_mask, events.Jet.pt, ascending=False),
            },
        },
        aux={
            "n_jets": ak.sum(jet_mask, axis=1),
        },
    )


#
# exposed selectors
# (those that can be invoked from the command line)
#

@selector(
    uses={
        # selectors / producers called within _this_ selector
        mc_weight, cutflow_features, process_ids, muon_selection, electron_selection,
        jet_selection, increment_stats,
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        mc_weight, cutflow_features, process_ids,
    },
    exposed=True,
)
def example(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # muon selection
    events, muon_results = self[muon_selection](events, **kwargs)
    results += muon_results

    # electron selection
    events, electron_results = self[electron_selection](events, **kwargs)
    results += electron_results

    # jet selection
    events, jet_results = self[jet_selection](events, **kwargs)
    results += jet_results

    # combined event selection after all steps
    results.event = results.steps.muon & results.steps.jet

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

    return events, results


@selector(
    uses={
        # selectors / producers called within _this_ selector
        mc_weight, cutflow_features, process_ids,
        "Jet.pt",
        increment_stats,
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        mc_weight, cutflow_features, process_ids,
    },
    exposed=True,
)
def empty(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    # prepare the selection results that are updated at every step
    results = SelectionResult(event=ak.Array(np.ones(len(events), dtype=np.bool_)))

    # jet selection
    events, jet_results = self[jet_selection](events, **kwargs)
    results += jet_results


    # create process ids
    if self.process_ids_dy is not None:
        events = self[self.process_ids_dy](events, **kwargs)
    else:
        events = self[process_ids](events, **kwargs)
   

    # add the mc weight
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)

    #

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
        }
        """
            # per jet multiplicity
              "njet": {
                "values": results.x.n_jets,
                "mask_fn": (lambda v: results.x.n_jets == v),
            }, """
        
    events, results = self[increment_stats](
        events,
        results,
        stats,
        weight_map=weight_map,
        group_map=group_map,
        **kwargs,
    )

    # sort jets
    results.objects.update(
        {
            "Jet": {
                "Jet": sorted_indices_from_mask(ak.ones_like(events.Jet.pt, dtype=bool), events.Jet.pt, ascending=False),
            },
        },
    )
    return events, results
