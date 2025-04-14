from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from hhh4b2tau.production.util import table_combo, min_chi_sqr_pair
import functools

ak = maybe_import("awkward")
np = maybe_import("numpy")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
class _HiggsReconstructor(Producer):
    
    def init_func(self):
        self.uses = {"Jet.{pt,eta,mass,phi}", attach_coffea_behavior}
        self.produces = {
            "BB{1,2}_idx",
            optional("min_chisq"), optional("rando_mask"), 
            attach_coffea_behavior,
            }

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
    
    # from IPython import embed; embed(header="load_jet_comb")
    events = self[attach_coffea_behavior](events, **kwargs,)
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
    bb1_idx = ak.drop_none(bb1_idx)
    bb2_idx = ak.drop_none(bb2_idx)

    events = set_ak_column(events, 'BB1_idx', bb1_idx)
    events = set_ak_column(events, 'BB2_idx', bb2_idx)

    return events

@_HiggsReconstructor.producer()
def higgs_reco_chi2(self, events: ak.Array, **kwargs,):
    events = self[attach_coffea_behavior](events, **kwargs)
    self.load_jet_combinations(events)
    bb1, bb2 , min_chisq, rando_mask = min_chi_sqr_pair(events.Jet, self.jet_table_combo)

    # insert pairs for case < 4 jets
    bb1 = ak.where(self.jet_num_mask, self.lone_pair, bb1)
    # extract indicies
    bb1_idx = ak.concatenate([ak.unflatten(bb1.idx1,1),ak.unflatten(bb1.idx2,1)],axis=1)
    bb2_idx = ak.concatenate([ak.unflatten(bb2.idx1,1),ak.unflatten(bb2.idx2,1)],axis=1)
    bb1_idx = ak.drop_none(bb1_idx)
    bb2_idx = ak.drop_none(bb2_idx)

    events = set_ak_column(events, "BB1_idx", bb1_idx)
    events = set_ak_column(events, "BB2_idx", bb2_idx)
    events = set_ak_column_f32(events, "min_chisq", min_chisq)
    events = set_ak_column(events, "rando_mask", rando_mask)

    return events