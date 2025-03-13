from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column

from hhh4b2tau.production.util import table_combo, min_chi_sqr_pair

ak = maybe_import("awkward")


class _HiggsReconstructor(Producer):
    
    def init_func(self):
        self.uses = {"Jet.{pt,eta,mass,phi}"}
        self.produces = {"BB1_idx", "BB2_idx"}

    def load_jet_combinations(self, events: ak.Array):
        self.jet_table_combo = table_combo(events.Jet)
        self.jet_optimal_mass_diff1 = ak.min(
            ak.min(self.jet_table_combo.mass_diff, axis=-1),
            axis=-1
        )
        self.jet_min_diff_mask1 = self.jet_optimal_mass_diff1 == self.jet_table_combo.mass_diff
        self.masked_table_combo = ak.mask(
            self.jet_table_combo, self.jet_table_combo.mass_diff==self.jet_optimal_mass_diff1
        )
        self.lone_pair = ak.firsts(
        ak.flatten(
                ak.drop_none(
                    ak.mask(self.jet_table_combo, self.jet_min_diff_mask1),
                    axis=1,
                ),
                axis=2,
            )
        )
        self.jet_num_mask = (ak.num(events.Jet) < 4) & (ak.num(events.Jet) > 1)

@_HiggsReconstructor.producer()
def higgs_reco_mass_diff(self, events: ak.Array, **kwargs):
    # get indicies of jet pair mass closest to 125 
    
    
    self.load_jet_combinations(events)
    jet_min_diff_idx1 = self.masked_table_combo.idx1
    jet_min_diff_idx2 = self.masked_table_combo.idx2
    #reshape into single entry arrays
    jet_min_diff_idx1 = ak.sum(ak.sum(ak.sum(
        ak.singletons(jet_min_diff_idx1,axis=-1),axis=-1),axis=-1),axis=-1)
    jet_min_diff_idx2 = ak.sum(ak.sum(ak.sum(
        ak.singletons(jet_min_diff_idx2,axis=-1),axis=-1),axis=-1),axis=-1)
    # create mask to remove already used jets
    jet_idx_mask = ((self.jet_table_combo.idx1 != jet_min_diff_idx1) &
                    (self.jet_table_combo.idx2 != jet_min_diff_idx2) &
                    (self.jet_table_combo.idx1 != jet_min_diff_idx2) &
                    (self.jet_table_combo.idx2 != jet_min_diff_idx1) )
    
    jet_massdiff_table2 = ak.mask(self.jet_table_combo.mass_diff, jet_idx_mask)
    jet_optimal_mass_diff2 = ak.min(ak.min(jet_massdiff_table2, axis=-1),axis=-1)
    jet_min_diff_mask2 = jet_optimal_mass_diff2 == jet_massdiff_table2
    jet_min_diff_mask2 = ak.fill_none(jet_min_diff_mask2, False, axis=-1)
    # mask that only contains "optimal" combinations 
    final_jet_mask = (self.jet_min_diff_mask1 | jet_min_diff_mask2)
    final_jet_table = ak.mask(self.jet_table_combo,final_jet_mask)
    final_jet_table = ak.flatten(ak.drop_none(final_jet_table, axis=-1),axis=-1)

    # for now use ascending in mass_diff
    sorted_jet_idx = ak.argsort(final_jet_table.mass_diff, axis=-1, ascending=True)

    final_jet_table = final_jet_table[sorted_jet_idx]

    # final bb pairings
    bb1 = final_jet_table[:,0]
    bb2 = final_jet_table[:,1]
    # for <4 jets insert only jet pair afterwards
    
    bb1 = ak.where(self.jet_num_mask, self.lone_pair, bb1)

    bb1_idx = ak.concatenate([ak.unflatten(bb1.idx1,1),ak.unflatten(bb1.idx2,1)],axis=1)
    bb2_idx = ak.concatenate([ak.unflatten(bb2.idx1,1),ak.unflatten(bb2.idx2,1)],axis=1)

    events = set_ak_column(events, 'BB1_idx', bb1_idx)
    events = set_ak_column(events, 'BB2_idx', bb2_idx)

    return events

@_HiggsReconstructor.producer()
def higgs_reco_chi2(self, events: ak.Array, **kwargs):

    self.load_jet_combinations(events)
    jet_chi_table = min_chi_sqr_pair(events.Jet, self.jet_table_combo)
    bb1_chi = jet_chi_table[:,0]
    bb2_chi = jet_chi_table[:,1]
    # insert pairs for case < 4 jets
    bb1_chi = ak.where(self.jet_num_mask, self.lone_pair, bb1_chi)

    bb1_chi_idx = ak.concatenate([ak.unflatten(bb1_chi.idx1,1),ak.unflatten(bb1_chi.idx2,1)],axis=1)
    bb2_chi_idx = ak.concatenate([ak.unflatten(bb2_chi.idx1,1),ak.unflatten(bb2_chi.idx2,1)],axis=1)

    events = set_ak_column(events, 'BB1_idx', bb1_chi_idx)
    events = set_ak_column(events, 'BB2_idx', bb2_chi_idx)