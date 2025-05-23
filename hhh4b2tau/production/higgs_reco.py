from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from hhh4b2tau.production.util import (table_combo, min_func_pair, dhh, chi2, mds)
import functools

ak = maybe_import("awkward")
np = maybe_import("numpy")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
class _HiggsReconstructor(Producer):
    
    def init_func(self):
        self.uses = {'Jet.{pt,eta,mass,phi,btagDeepFlavB}', attach_coffea_behavior}
        self.produces = {
            'BB{1,2}_idx',
            attach_coffea_behavior,
            }

    def load_jet_combinations(self, events: ak.Array):

        # b-tagged jet idx
        btag_wp = self.config_inst.x.btag_working_points.deepjet.medium
        btag_mask = events.Jet.btagDeepFlavB >= btag_wp
        self.n_bjet =ak.sum(btag_mask, axis =1)
        # btag_idx = ak.local_index(btag_mask)[btag_mask]
        # non b-tagged jet idx with highest pt
        non_btag_idx = ak.local_index(btag_mask)[~btag_mask][:, :1]
        # select jets with highest btag score over btag_wp
        # case == 3 btags add the remaining jet with highest pt
        sorted_btag_idx = ak.argsort(events.Jet.btagDeepFlavB, axis=1, ascending = False)
        sorted_btag_mask = (events.Jet.btagDeepFlavB[sorted_btag_idx] >= btag_wp)
        sorted_btag_idx = sorted_btag_idx[sorted_btag_mask]
        self.reco_idx = ak.concatenate([sorted_btag_idx, non_btag_idx], axis=1)[:, :4]

        self.jet_table_combo = table_combo(array = events.Jet, idx = self.reco_idx)
        self.jet_optimal_chi2_1 = ak.min(
            ak.min(self.jet_table_combo.chi2, axis=-1),
            axis=-1
        )
        self.jet_min_chi2_mask1 = self.jet_optimal_chi2_1 == self.jet_table_combo.chi2
        self.masked_table_combo = ak.mask(
            self.jet_table_combo, self.jet_table_combo.chi2==self.jet_optimal_chi2_1
        )

        self.lone_pair = ak.flatten(
                ak.drop_none(
                    ak.mask(self.jet_table_combo, self.jet_min_chi2_mask1),
                    axis=1,
                ),
                axis=2,
            )

        self.lone_pair_idx = ak.concatenate([self.lone_pair.idx1, self.lone_pair.idx2], axis=1)

        self.jet_num_mask = (ak.num(self.reco_idx) < 4) & (ak.num(self.reco_idx) > 1)

@_HiggsReconstructor.producer()     
def higgs_reco_mds(self, events: ak.Array, **kwargs):
    # get indicies of jet pair mass closest to 125 
    events = self[attach_coffea_behavior](events, **kwargs,)
    self.load_jet_combinations(events)

    bb1_idx, bb2_idx , min_func_val = min_func_pair(
        mds, 
        events.Jet, 
        self.jet_table_combo, 
        self.reco_idx, 
        ordered=True
        )

    # insert pairs for case < 4 jets
    bb1_idx = ak.where(self.jet_num_mask, self.lone_pair_idx, bb1_idx)

    # special treatment because of two seperate values that are minimized
    chi2_1 = chi2_12(events.Jet[bb1_idx]) 
    chi2_2 = chi2_12(events.Jet[bb2_idx])

    events = set_ak_column(events, 'BB1_idx', bb1_idx)
    events = set_ak_column(events, 'BB2_idx', bb2_idx)
    events = set_ak_column_f32(events, 'chi2_1', chi2_1)
    events = set_ak_column_f32(events, 'chi2_2', chi2_2)
    return events

@higgs_reco_mds.init
def higgs_reco_mds_init(self: Producer) -> None:
    super(higgs_reco_mds, self).init_func()
    self.produces |= {'chi2_{1,2}'}


def chi2_12(vectors):
    sum = ak.firsts(vectors[:, :1]) + ak.firsts(vectors[:, 1:2])
    chi2 = ((sum.mass-125)/125)**2
    return chi2


@_HiggsReconstructor.producer()
def higgs_reco_chi2(self, events: ak.Array, **kwargs,):
    events = self[attach_coffea_behavior](events, **kwargs)
    self.load_jet_combinations(events)
    bb1_idx, bb2_idx , min_func_val = min_func_pair(
        chi2, 
        events.Jet, 
        self.jet_table_combo, 
        self.reco_idx, 
        ordered=False
        )

    # insert pairs for case < 4 jets
    bb1_idx = ak.where(self.jet_num_mask, self.lone_pair_idx, bb1_idx)

    events = set_ak_column(events, 'BB1_idx', bb1_idx)
    events = set_ak_column(events, 'BB2_idx', bb2_idx)
    events = set_ak_column_f32(events, "chi2", min_func_val)

    return events

@higgs_reco_chi2.init
def higgs_reco_chi2_init(self: Producer) -> None:
    super(higgs_reco_chi2, self).init_func()
    self.produces |= {'chi2'}



@_HiggsReconstructor.producer()
def higgs_reco_dhh(self, events: ak.Array, **kwargs,):
    events = self[attach_coffea_behavior](events, **kwargs)
    self.load_jet_combinations(events)
    # bb1, bb2 , min_chisq = min_chi_sqr_pair(events.Jet, self.jet_table_combo, self.reco_idx)
    bb1_idx, bb2_idx , min_func_val = min_func_pair(
        dhh, 
        events.Jet, 
        self.jet_table_combo, 
        self.reco_idx, 
        ordered=True
        )

    # # insert pairs for case < 4 jets
    # bb1_idx = ak.where(self.jet_num_mask, self.lone_pair_idx, bb1_idx)

    events = set_ak_column(events, 'BB1_idx', bb1_idx)
    events = set_ak_column(events, 'BB2_idx', bb2_idx)
    events = set_ak_column_f32(events, 'dhh', min_func_val)

    return events

@higgs_reco_dhh.init
def higgs_reco_dhh_init(self: Producer) -> None:
    super(higgs_reco_dhh, self).init_func()
    self.produces |= {'dhh'}