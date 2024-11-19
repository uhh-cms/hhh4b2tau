from collections import defaultdict

from operator import or_
from functools import reduce

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.cms.jets import jet_veto_map
from columnflow.selection.stats import increment_stats
from columnflow.selection.util import sorted_indices_from_mask
from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column
from hhh4b2tau.util import IF_RUN_2, IF_RUN_3

from hhh4b2tau.production.example import cutflow_features

np = maybe_import("numpy")
ak = maybe_import("awkward")

@selector(
    uses={ IF_RUN_3(jet_veto_map), IF_RUN_2("Jet.puId"),
        #   "trigger_ids",
          "Jet.{pt,eta,phi,mass,jetId,btagDeepFlavB}",
          "FatJet.{pt,eta,phi,mass,msoftdrop,jetId,subJetIdx1,subJetIdx2}",
          "SubJet.{pt,eta,phi,mass,btagDeepB}",
          },
)
def jet_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    from IPython import embed; embed()

    # local jet index
    # li = ak.local_index(events.Jet)

    # common ak4 jet mask for normal and vbf jets
    ak4_mask = (
        (events.Jet.jetId == 6) # &  # tight plus lepton veto
        # ak.all(events.Jet.metric_table(lepton_results.x.lepton_pair) > 0.5, axis=2)
    )

    jet_wp = self.config_inst.x.btag_working_points.deepjet.medium
    jets_btagged = events.Jet.btagDeepFlavB > jet_wp

    # default jets
    default_mask = (
        ak4_mask &
        (events.Jet.pt > 20.0) &
        (abs(events.Jet.eta) < 2.4) &
        jets_btagged
    )

    fatjet_mask = (
        (events.FatJet.jetId == 6) &  # tight plus lepton veto
        (events.FatJet.msoftdrop > 30.0) &
        (abs(events.FatJet.eta) < 2.4) &
        (events.FatJet.subJetIdx1 >= 0) &
        (events.FatJet.subJetIdx2 >= 0)
    )

    # store fatjet and subjet indices
    fatjet_indices = ak.local_index(events.FatJet.pt)[fatjet_mask]
    subjet_indices = ak.concatenate(
        [
            events.FatJet[fatjet_mask].subJetIdx1[..., None],
            events.FatJet[fatjet_mask].subJetIdx2[..., None],
        ],
        axis=2,
    )



    # discard the event in case the (first) fatjet with matching subjets is found
    # but they are not b-tagged (TODO: move to deepjet when available for subjets)
    if self.config_inst.campaign.x.run == 3:
        wp = self.config_inst.x.btag_working_points.particleNet.loose
    else:
        wp = self.config_inst.x.btag_working_points.deepcsv.loose

    subjets_btagged = ak.all(events.SubJet[ak.firsts(subjet_indices)].btagDeepB > wp, axis=1)

    # final event selection
    jet_sel2 = (
        (ak.sum(default_mask, axis=1) >= 2) &
        ak.fill_none(subjets_btagged, True) # was none for events with no matched fatjet
    )

    jet_sel3 = (
        (ak.sum(default_mask, axis=1) >= 3) &
        ak.fill_none(subjets_btagged, True)
    )

    jet_sel4 = (
        (ak.sum(default_mask, axis=1) >= 3) &
        ak.fill_none(subjets_btagged, True)
    )

    # pt sorted indices to convert mask
    jet_indices = sorted_indices_from_mask(default_mask, events.Jet.pt, ascending=False)

    # some final type conversions
    jet_indices = ak.values_astype(ak.fill_none(jet_indices, 0), np.int32)
    fatjet_indices = ak.values_astype(fatjet_indices, np.int32)

    # additional jet veto map, vetoing entire events
    # if self.has_dep(jet_veto_map):
    #     events, veto_result = self[jet_veto_map](events, **kwargs)
    #     result += veto_result

    # build and return selection results
    # "objects" maps source columns to new columns and selections to be applied on the old columns
    # to create them, e.g. {"Jet": {"MyCustomJetCollection": indices_applied_to_Jet}}
    return events, SelectionResult(
        steps={
            "two_bjet" : jet_sel2,
            "three_bjet" : jet_sel3,
            "four_bjet": jet_sel4,
        },
        objects={
            "Jet": {
                "Jet": jet_indices,
            },
            "FatJet": {
                "FatJet": fatjet_indices,
            },
            "SubJet": {
                "SubJet1": subjet_indices[..., 0],
                "SubJet2": subjet_indices[..., 1],
            },
        },
        aux={
            # jet mask that lead to the jet_indices
            "jet_mask": default_mask,
        },
    )