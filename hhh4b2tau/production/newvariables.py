import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, set_ak_column
from columnflow.production.util import attach_coffea_behavior
from columnflow.columnar_util import (
    set_ak_column, remove_ak_column, attach_behavior, EMPTY_FLOAT, get_ak_routes, remove_ak_column,
    optional_column as optional
)
# coffea cannot be imported????
# from coffea.nanoevents.methods import vector 
# # creates a lorentz vector with all zeros
# zero_lorentz = ak.zip(
#         {
#             "x": 0.,
#             "y": 0.,
#             "z": 0.,
#             "t": 0.,
#         },
#         with_name="LorentzVector",
#         behavior=vector.behavior,
#     )

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses=(
        {
            f"Jet.{var}"
            for var in [
                'pt', 'eta', 'phi', 'mass'
            ]
        } | {attach_coffea_behavior}
    ),
    produces={
        "jet_delta_phi",
        "jet_delta_r_12",
        "jet_delta_r_13",
    },
)
def jet_angle_difference(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Calculate the difference of the azimutal angle between the two strongest jets
    """
    # attach coffea behavior for four-vector arithmetic
    events = self[attach_coffea_behavior](
        events,
        **kwargs,
    )
    n_jets = (
        ak.num(events.Jet.phi , axis=1) 
    )
    dijet_mask = n_jets >=2
    trijet_mask = n_jets >=3
    # calculate delta phi
    jets = ak.where(
                dijet_mask,
                events.Jet.phi,
                [[EMPTY_FLOAT]*2]*len(n_jets)
                )
    
    delta = jets[:,0]-jets[:,1]

    events = set_ak_column_f32(
        events,
        "jet_delta_phi",
        delta,
    )
    # calculate delta r between all jets
    all_delta_r = events.Jet.metric_table(events.Jet)
    # get delta r values of first (hardest) jet to all other jets
    hardest_delta_r = ak.firsts(all_delta_r)
    # get delta r values between hardest and 2nd hardest jet
    events = set_ak_column_f32(
        events,
        "jet_delta_r_12",
        ak.where(dijet_mask, ak.mask(hardest_delta_r, dijet_mask)[:, 1], EMPTY_FLOAT),
    )

    # get delta r values between hardest and 3rd hardest jet
    events = set_ak_column_f32(
        events,
        "jet_delta_r_13",
        ak.where(trijet_mask, ak.mask(hardest_delta_r, trijet_mask)[:, 2], EMPTY_FLOAT),
    )

    
    # return the events
    return events


@producer(
    uses=(
        {
            
        optional(f"gen_{mother}_to_{child}.{var}")
        for mother in ('h', 'tau', )
        for child in ('b', 'tau', 'taunu', 'electron', 'enu', 'muon', 'munu', )
        for var in ('pt', 'eta', 'phi', 'mass', 'pdgId', )
    } |
    {   
        optional(f"gen_{child}.{var}")
        for child in ('b', 'tau', 'taunu', 'electron', 'enu', 'muon', 'munu', )
        for var in ('pt', 'eta', 'phi', 'mass', 'pdgId', )} | 
        {
            attach_coffea_behavior,
        }
    ),
    produces={
        "mtautau_gen", "mbb_gen", "mhhh_gen", "mlnu_gen", "hpt_gen", "h1bpt_gen",
        "h2bpt_gen", "htaupt_gen",
        "delta_r_h12_gen", "delta_r_h13_gen", "delta_r_h23_gen", "delta_r_bb1_gen", 
        "delta_r_bb2_gen",
        "delta_r_tautau_gen",
        "cos_h12_gen", "cos_h13_gen", "cos_h23_gen", "cos_bb1_gen", "cos_bb2_gen",
        "cos_tautau_gen",
    },
)
def hhh_decay_invariant_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Construct b and tau invariant mass. 
    Also added pt, delta R and cos(delta) of H and their decay products.
    All on Gen-level
    """


    # attach coffea behavior for four-vector arithmetic
    events = self[attach_coffea_behavior](
        events,
        collections={ x : {
                "type_name": "GenParticle",
            } for x in [
            "gen_h_to_b", 
            "gen_b", 
            "gen_h_to_tau", 
            "gen_tau",
            "gen_taunu",
            "gen_electron",
            "gen_enu",
            "gen_munu",
            "gen_muon",
            ]},
        **kwargs,
    )

    # from IPython import embed; embed()

    # total number of partons per event
    n_taus = ak.num(events.gen_tau, axis=-1)
    n_bs = ak.num(events.gen_b, axis=-1)
    n_hs = (ak.num(events.gen_h_to_tau, axis=-1) +
    ak.num(events.gen_h_to_b, axis=-1))
  
    # total number of leptons, neutrinos per event
    n_taunu = ak.num(events.gen_taunu, axis=1)
    n_electron = ak.num(events.gen_electron, axis=1)
    n_muon = ak.num(events.gen_muon, axis=1)
    n_lep = n_electron + n_muon   

    # four-vector sum of first four elements of each
    # tau collection (possibly fewer)
    
    ditau = events.gen_tau.sum(axis=-1)
    dib = events.gen_b.sum(axis=-1)
    trih = events.gen_h_to_tau.sum(axis=-1) + events.gen_h_to_b.sum(axis=-1)
    # create 0 lorentz-vector
    zero_lorentz = ditau[0]-ditau[0]
    # lorentz vector sum of tau decay products
    ditaunu = ak.where(n_taunu>=1, ak.flatten(events.gen_taunu.sum(axis=1)), zero_lorentz)
    dielectron = ak.fill_none(ak.pad_none(events.gen_electron.sum(axis=1), 1), zero_lorentz[0])
    dienu = ak.fill_none(
        ak.pad_none(events.gen_enu.sum(axis=1), 1),
          zero_lorentz[0])
    dimuon = ak.fill_none(
        ak.pad_none(events.gen_muon.sum(axis=1), 1),
          zero_lorentz[0])
    dimunu = ak.fill_none(
        ak.pad_none(events.gen_munu.sum(axis=1), 1),
          zero_lorentz[0])

    # lep_sum = dielectron + dimuon
    # nu_sum = ditaunu + dienu + dimunu
    lep_nu_sum = ditaunu + dienu + dimunu + dielectron + dimuon

    # pt of all Higgs
    h_pt = ak.concatenate(
        (events.gen_h_to_b.pt, events.gen_h_to_tau.pt), axis=1)
    h1b_pt = events.gen_h_to_b.pt[:,0]
    h2b_pt = events.gen_h_to_b.pt[:,1]
    htau_pt = events.gen_h_to_tau.pt

    # delta_r between all Higgs
    delta_r_h12 = events.gen_h_to_b[:,0].delta_r(events.gen_h_to_b[:,1])
    delta_r_h13 = events.gen_h_to_b[:,0].delta_r(events.gen_h_to_tau)
    delta_r_h23 = events.gen_h_to_b[:,1].delta_r(events.gen_h_to_tau)

    # delta_r between H decay products
    delta_r_bb1 = events.gen_b[:,0,0].delta_r(events.gen_b[:,0,1])
    delta_r_bb2 = events.gen_b[:,1,0].delta_r(events.gen_b[:,1,1])
    delta_r_tautau = events.gen_tau[:,0,0].delta_r(events.gen_tau[:,0,1])
    # from IPython import embed; embed()

    # cosine of opening angle little delta between h and between their decay products 
    cos_h12 = events.gen_h_to_b[:,0].unit.dot(events.gen_h_to_b[:,1].unit)
    cos_h13 = events.gen_h_to_b[:,0].unit.dot(events.gen_h_to_tau.unit)
    cos_h23 = events.gen_h_to_b[:,1].unit.dot(events.gen_h_to_tau.unit)
    cos_bb1 = events.gen_b[:,0,0].unit.dot(events.gen_b[:,0,1].unit)
    cos_bb2 = events.gen_b[:,1,0].unit.dot(events.gen_b[:,1,1].unit)
    cos_tautau = events.gen_tau[:,0,0].unit.dot(events.gen_tau[:,0,1].unit)

    # four-lepton mass, taking into account only events with at least four leptons,
    # and otherwise substituting a predefined EMPTY_FLOAT value
    tautau_mass = ak.where(
        n_taus >= 2,
        ditau.mass,
        EMPTY_FLOAT,
    )

    b_mass = ak.where(
        n_bs >= 2,
        dib.mass,
        EMPTY_FLOAT,
    )

    h_mass = ak.where(
        n_hs >= 3,
        trih.mass,
        EMPTY_FLOAT,
    )

    lep_nu_mass = ak.where(
        n_lep >= 2,
        lep_nu_sum.mass,
        [[EMPTY_FLOAT]],
    )

    # write out the resulting mass to the `events` array,
    events = set_ak_column_f32(
        events,
        "mtautau_gen",
        tautau_mass,
    )

    events = set_ak_column_f32(
        events,
        "mbb_gen",
        b_mass,
    )

    events = set_ak_column_f32(
        events,
        "mhhh_gen",
        h_mass,
    )

    events = set_ak_column_f32(
        events,
        "mlnu_gen",
        lep_nu_mass,
    )
    # pt of h and their decay products 
    events = set_ak_column_f32(
        events,
        "hpt_gen",
        h_pt,
    )

    events = set_ak_column_f32(
        events,
        "htaupt_gen",
        htau_pt,
    )

    events = set_ak_column_f32(
        events,
        "h1bpt_gen",
        h1b_pt,
    )

    events = set_ak_column_f32(
        events,
        "h2bpt_gen",
        h2b_pt,
    )
    # delta r between h and between their decay products

    events = set_ak_column_f32(
        events,
        "delta_r_h12_gen",
        delta_r_h12,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h13_gen",
        delta_r_h13,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h23_gen",
        delta_r_h23,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_bb1_gen",
        delta_r_bb1,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_bb2_gen",
        delta_r_bb2,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_tautau_gen",
        delta_r_tautau,
    )

    # cos(delta) between h and between their decay products

    events = set_ak_column_f32(
        events,
        "cos_h12_gen",
        cos_h12,
    )

    events = set_ak_column_f32(
        events,
        "cos_h13_gen",
        cos_h13,
    )

    events = set_ak_column_f32(
        events,
        "cos_h23_gen",
        cos_h23,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb1_gen",
        cos_bb1,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb2_gen",
        cos_bb2,
    )

    events = set_ak_column_f32(
        events,
        "cos_tautau_gen",
        cos_tautau,
    )
    # return the events
    # from IPython import embed; embed()
    return events


@producer(
    uses=(
    {   
        optional(f"gen_tth_{child}.{var}")
        for child in ('b1', 'b2', 'tau', 'taunu',)
        for var in ('pt', 'eta', 'phi', 'mass', 'pdgId', )} | 
        {
            attach_coffea_behavior,
        }
    ),
    produces={
        "delta_r_bb1_gen", "delta_r_bb2_gen",
        "delta_r_tautau_gen",
        "cos_bb1_gen", "cos_bb2_gen",
        "cos_tautau_gen",
        "mhhh_gen",
    },
)
def tth_variables(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    '''
    migrating variables previously found in hhh_decay_invariant_mass 
    that are also compatible with background data from tth
    '''
    # attach lorentz vector behavior
    events = self[attach_coffea_behavior](
        events,
        collections={ x : {
                "type_name": "GenParticle",
            } for x in [
            "gen_tth_b1",
            "gen_tth_b2",
            "gen_tth_tau",
            "gen_tth_taunu",
            ]},
        **kwargs,
    )
    '''
    still work in progress since array dimensions differ from hhh processes
    ''' 
    # from IPython import embed; embed(header='tth_variables')

    # masks to require each come in pairs at least
    # b1_mask = ak.flatten(ak.num(events.gen_tth_b1, axis=-1) >= 2)
    # b2_mask = ak.num(events.gen_tth_b2, axis=-2) >= 2
    tau_mask = ak.num(events.gen_tth_tau, axis=-2) >= 2

    # (optionally) apply masks beforehand to avoid errors
    # for b1,b2 not necessary since the columns are always full 
    b1 = ak.flatten(events.gen_tth_b1)
    b2 = events.gen_tth_b2
    tau = ak.mask(events.gen_tth_tau, tau_mask)


    # delta_r between bb tautau
    delta_r_bb1 = b1[:,0].delta_r(b1[:,1])
    delta_r_bb2 = b2[:,0].delta_r(b2[:,1])
    delta_r_tautau = tau[:,0].delta_r(tau[:,1])

    # work still in progress since array dimensions differ from hhh
    cos_bb1 = b1[:,0].unit.dot(b1[:,1].unit)
    cos_bb2 = b2[:,0].unit.dot(b2[:,1].unit)
    cos_tautau = tau[:,0].unit.dot(tau[:,1].unit)

    # invariant mass of all 4b2tau but still call it mhhh for comparison
    mhhh = (b1.sum(axis=-1) + b2.sum(axis=-2) + tau.sum(axis=-2)).mass

    events = set_ak_column_f32(
        events,
        "delta_r_bb1_gen",
        delta_r_bb1,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_bb2_gen",
        delta_r_bb2,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_tautau_gen",
        delta_r_tautau,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb1_gen",
        cos_bb1,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb2_gen",
        cos_bb2,
    )

    events = set_ak_column_f32(
        events,
        "cos_tautau_gen",
        cos_tautau,
    )

    events = set_ak_column_f32(
        events,
        "mhhh_gen",
        mhhh,
    )

    return events

# variables on hadronic level for ggf hhh processes
@producer(
    uses=({"gen_b_jet.*", "GenVisTau.*", # "GenMET.*", "gen_taunu.*",
            }
          | {attach_coffea_behavior}
          ),
    produces={"mhhh_hadron", 
              "delta_r_bb1_hadron", "delta_r_bb2_hadron","delta_r_tautau_hadron",
              "cos_bb1_hadron", "cos_bb2_hadron","cos_tautau_hadron",
              "delta_r_h12_hadron", "delta_r_h13_hadron", "delta_r_h23_hadron",
              "cos_h12_hadron", "cos_h13_hadron", "cos_h23_hadron",
        },
        
)
def genHadron_variables(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events = self[attach_coffea_behavior](
        events,
        collections={ "GenVisTau" : {
                "type_name": "GenVisTau",

            }},
        **kwargs,
    )

    events = self[attach_coffea_behavior](
        events,
        collections={ "gen_b_jet" : {
                "type_name": "Jet",
        }},
        **kwargs,
    )

    # events = self[attach_coffea_behavior](
    #     events,
    #     collections={ "gen_taunu" : {
    #             "type_name": "GenParticle",
    #         }},
    #     **kwargs,
    # )

    # b_jet_mask = ak.num(events.gen_b_jet, axis=-1) >= 4
    # b_jet = ak.mask(events.gen_b_jet, b_jet_mask)
    # tau_mask = ak.num(events.GenVisTau, axis=-1) >= 2
    # tau = ak.mask(events.GenVisTau, tau_mask)

    b_jet = events.gen_b_jet
    tau = events.GenVisTau

    # taunu_mask = ak.num(events.gen_taunu, axis=-2) == 2 # this is always true for signal samples
    # taunu = ak.mask(events.gen_taunu, taunu_mask)
    # taunu_sum = ak.flatten(taunu.sum(axis=-2))



    # look out for all permutations of b quarks
    # there are up to 8 b jets and 5 tau

    # metric table of all possible lorentz vector sums of two jets
    b_jet_table = b_jet.metric_table(b_jet, metric=lambda a, b: (a+b))
    b_jet_massdiff_table1 = b_jet.metric_table(b_jet, 
                            metric=lambda a, b: abs((a+b).mass-125))
    # metric_table for delta_r and cos(delta)
    b_jet_delta_r_table = b_jet.metric_table(b_jet)
    b_jet_cos_table = b_jet.metric_table(b_jet, metric=lambda a, b: a.pvec.dot(b.pvec)/(a.pvec.absolute()*b.pvec.absolute()))

    # create dictionary for b jets with indicies
    b_idx1 = ak.local_index(b_jet_table, axis=-2)
    b_idx2 = ak.local_index(b_jet_table, axis=-1)
    b_table_combo = ak.zip({"pair_sum": b_jet_table, 
                            "mass_diff": b_jet_massdiff_table1,
                            "delta_r": b_jet_delta_r_table,
                            "cos": b_jet_cos_table, 
                            "idx1": b_idx1, "idx2": b_idx2})
    # remove duplicates and self sums
    b_table_combo = ak.mask(b_table_combo,b_table_combo.idx1<b_table_combo.idx2)
    # get indicies of jet pair mass closest to 125 
    b_optimal_mass_diff1 = ak.min(ak.min(b_table_combo.mass_diff, axis=-1),axis=-1)
    b_min_diff_mask1 = b_optimal_mass_diff1 == b_table_combo.mass_diff
    b_min_diff_idx1 = ak.mask(b_table_combo,b_table_combo.mass_diff==b_optimal_mass_diff1).idx1
    b_min_diff_idx2 = ak.mask(b_table_combo,b_table_combo.mass_diff==b_optimal_mass_diff1).idx2
    #reshape into single entry arrays
    b_min_diff_idx1 = ak.sum(ak.sum(ak.sum(
        ak.singletons(b_min_diff_idx1,axis=-1),axis=-1),axis=-1),axis=-1)
    b_min_diff_idx2 = ak.sum(ak.sum(ak.sum(
        ak.singletons(b_min_diff_idx2,axis=-1),axis=-1),axis=-1),axis=-1)
    # create mask to remove already used jets
    b_idx_mask = ((b_table_combo.idx1 != b_min_diff_idx1) &
                  (b_table_combo.idx2 != b_min_diff_idx2) &
                  (b_table_combo.idx1 != b_min_diff_idx2) &
                  (b_table_combo.idx2 != b_min_diff_idx1))
    
    b_jet_massdiff_table2 = ak.mask(b_jet_massdiff_table1, b_idx_mask)
    b_optimal_mass_diff2 = ak.min(ak.min(b_jet_massdiff_table2, axis=-1),axis=-1)
    b_min_diff_mask2 = b_optimal_mass_diff2 == b_jet_massdiff_table2
    b_min_diff_mask2 = ak.fill_none(b_min_diff_mask2, False, axis=-1)
    # mask that only contains "optimal" combinations 
    final_b_jet_mask = (b_min_diff_mask1 | b_min_diff_mask2)
    final_b_jet_table = ak.mask(b_table_combo,final_b_jet_mask)
    final_b_jet_table = ak.flatten(ak.drop_none(final_b_jet_table, axis=-1),axis=-1)
    # sort bb pairs by pt
    sorted_b_jet_idx = ak.argsort(final_b_jet_table.pair_sum.pt, axis=-1, ascending=False)
    final_b_jet_table = final_b_jet_table[sorted_b_jet_idx]

    # now do the same for taus but simpler since only two are required

    # tau_taunu_massdiff_table = tau.metric_table((tau+taunu_sum), 
    #                      metric=lambda a, b: abs((a+b).mass -125))
    tau_massdiff_table = tau.metric_table(tau, 
                         metric=lambda a, b: abs((a+b).mass -125))
    tau_table = tau.metric_table((tau), metric=lambda a, b: (a+b))
    tau_delta_r_table = tau.metric_table(tau)
    tau_cos_table = tau.metric_table(tau, metric=lambda a, b: a.pvec.dot(b.pvec)/(a.pvec.absolute()*b.pvec.absolute()))
    tau_idx1 = ak.local_index(tau_table, axis=-2)
    tau_idx2 = ak.local_index(tau_table, axis=-1)
    tau_table_combo = ak.zip({"pair_sum": tau_table, 
                              "mass_diff": tau_massdiff_table, 
                              "delta_r": tau_delta_r_table,
                              "cos": tau_cos_table,
                              "idx1": tau_idx1, "idx2": tau_idx2})

    tau_table_combo = ak.mask(tau_table_combo,tau_table_combo.idx1<tau_table_combo.idx2)
    tau_optimal_mass_diff = ak.min(ak.min(tau_table_combo.mass_diff, axis=-1),axis=-1)
    tau_min_diff_mask = tau_optimal_mass_diff == tau_table_combo.mass_diff
    final_tau_table = ak.drop_none(ak.mask(tau_table_combo, tau_min_diff_mask),axis=-1)



    # final tautau pairing
    tautau = ak.flatten(final_tau_table,axis=-1)[:,0]

    # final bb pairings
    bb1 = final_b_jet_table[:,0]
    bb2 = final_b_jet_table[:,1]

    # reconstructed higgs h1 and h2 are pt sorted from b-jets and h3 from taus
    h1 = bb1.pair_sum * 1
    h2 = bb2.pair_sum * 1
    h3 = tautau.pair_sum * 1


    events = set_ak_column_f32(
        events,
        "mhhh_hadron",
        (h1 + h2 + h3).mass,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_bb1_hadron",
        bb1.delta_r,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb1_hadron",
        bb1.cos,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_bb2_hadron",
        bb2.delta_r,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb2_hadron",
        bb2.cos,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_tautau_hadron",
        tautau.delta_r,
    )

    events = set_ak_column_f32(
        events,
        "cos_tautau_hadron",
        tautau.cos,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h12_hadron",
        h1.delta_r(h2),
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h13_hadron",
        h1.delta_r(h3),
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h23_hadron",
        h2.delta_r(h3),
    )

    events = set_ak_column_f32(
        events,
        "cos_h12_hadron",
        h1.pvec.dot(h2.pvec)/(h1.pvec.absolute()*h2.pvec.absolute()),
    )

    events = set_ak_column_f32(
        events,
        "cos_h13_hadron",
        h1.pvec.dot(h3.pvec)/(h1.pvec.absolute()*h3.pvec.absolute()),
    )

    events = set_ak_column_f32(
        events,
        "cos_h23_hadron",
        h2.pvec.dot(h3.pvec)/(h2.pvec.absolute()*h3.pvec.absolute()),
    )

    return events


# producer for analysis on detector level


@producer(
    uses=(
    {   
        optional(f"{field}.{var}")
        for field in ["Jet", "Tau", "Electron", "Muon", "FatJet"]
        for var in ["pt", "eta", "phi", "mass", "hadronFlavour", "charge",]
        } | 
        {
            attach_coffea_behavior,
        }
    ),
    produces={
        "delta_r_bb1", "delta_r_bb2", "delta_r_tautau", 
        "delta_r_h12", "delta_r_h13", "delta_r_h23",
        "cos_bb1", "cos_bb2", "cos_tautau",
        "cos_h12", "cos_h13", "cos_h23",
        "mhhh", "h1_mass", "h2_mass", "h3_mass",
        "n_b_jet", "n_fatjet",
        "h1_unsort_mass", "h2_unsort_mass"
    },
)
def dectector_variables(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    events = self[attach_coffea_behavior](
        events,
        **kwargs,
    )



    b_jet = events.Jet

    tau = events.Tau
    # from IPython import embed; embed(header="inside detector variables")
    # now copy-paste from hadron analysis 
    # final result will automatically filter all cases < 4 b-jets and < 2 tau out

    # metric table of all possible lorentz vector sums of two jets
    b_jet_table = b_jet.metric_table(b_jet, metric=lambda a, b: (a+b))
    b_jet_massdiff_table1 = b_jet.metric_table(b_jet, 
                            metric=lambda a, b: abs((a+b).mass-125))
    # metric_table for delta_r and cos(delta)
    b_jet_delta_r_table = b_jet.metric_table(b_jet)
    # b_jet_cos_table = b_jet.unit.metric_table(b_jet.unit, metric=lambda a, b: a.dot(b)) # no longer valid for newest version of coffea
    b_jet_cos_table = b_jet.metric_table(b_jet, metric=lambda a, b: a.pvec.dot(b.pvec)/(a.pvec.absolute()*b.pvec.absolute()))
    # create dictionary for b jets with indicies
    b_idx1 = ak.local_index(b_jet_table, axis=-2)
    b_idx2 = ak.local_index(b_jet_table, axis=-1)
    b_table_combo = ak.zip({"pair_sum": b_jet_table, 
                            "mass_diff": b_jet_massdiff_table1,
                            "delta_r": b_jet_delta_r_table,
                            "cos": b_jet_cos_table, 
                            "idx1": b_idx1, "idx2": b_idx2})
    # remove duplicates and self sums
    b_table_combo = ak.mask(b_table_combo,b_table_combo.idx1<b_table_combo.idx2)
    # get indicies of jet pair mass closest to 125 
    b_optimal_mass_diff1 = ak.min(ak.min(b_table_combo.mass_diff, axis=-1),axis=-1)
    b_min_diff_mask1 = b_optimal_mass_diff1 == b_table_combo.mass_diff
    b_min_diff_idx1 = ak.mask(b_table_combo,b_table_combo.mass_diff==b_optimal_mass_diff1).idx1
    b_min_diff_idx2 = ak.mask(b_table_combo,b_table_combo.mass_diff==b_optimal_mass_diff1).idx2
    #reshape into single entry arrays
    b_min_diff_idx1 = ak.sum(ak.sum(ak.sum(
        ak.singletons(b_min_diff_idx1,axis=-1),axis=-1),axis=-1),axis=-1)
    b_min_diff_idx2 = ak.sum(ak.sum(ak.sum(
        ak.singletons(b_min_diff_idx2,axis=-1),axis=-1),axis=-1),axis=-1)
    # create mask to remove already used jets
    b_idx_mask = ((b_table_combo.idx1 != b_min_diff_idx1) &
                  (b_table_combo.idx2 != b_min_diff_idx2) &
                  (b_table_combo.idx1 != b_min_diff_idx2) &
                  (b_table_combo.idx2 != b_min_diff_idx1) )
    
    b_jet_massdiff_table2 = ak.mask(b_jet_massdiff_table1, b_idx_mask)
    b_optimal_mass_diff2 = ak.min(ak.min(b_jet_massdiff_table2, axis=-1),axis=-1)
    b_min_diff_mask2 = b_optimal_mass_diff2 == b_jet_massdiff_table2
    b_min_diff_mask2 = ak.fill_none(b_min_diff_mask2, False, axis=-1)
    # mask that only contains "optimal" combinations 
    final_b_jet_mask = (b_min_diff_mask1 | b_min_diff_mask2)
    final_b_jet_table = ak.mask(b_table_combo,final_b_jet_mask)
    final_b_jet_table = ak.flatten(ak.drop_none(final_b_jet_table, axis=-1),axis=-1)
    # sort bb pairs by pt
    sorted_b_jet_idx = ak.argsort(final_b_jet_table.pair_sum.pt, axis=-1, ascending=False)
    final_b_jet_table = final_b_jet_table[sorted_b_jet_idx]

    tau_massdiff_table = tau.metric_table(tau, 
                         metric=lambda a, b: abs((a+b).mass -125))
    tau_table = tau.metric_table((tau), metric=lambda a, b: (a+b))
    tau_delta_r_table = tau.metric_table(tau)
    # tau_cos_table = tau.unit.metric_table(tau.unit, metric=lambda a, b: a.dot(b))
    tau_cos_table = tau.metric_table(tau, metric=lambda a, b: a.pvec.dot(b.pvec)/(a.pvec.absolute()*b.pvec.absolute()))
    tau_idx1 = ak.local_index(tau_table, axis=-2)
    tau_idx2 = ak.local_index(tau_table, axis=-1)
    tau_table_combo = ak.zip({"pair_sum": tau_table, 
                              "mass_diff": tau_massdiff_table, 
                              "delta_r": tau_delta_r_table,
                              "cos": tau_cos_table,
                              "idx1": tau_idx1, "idx2": tau_idx2})

    tau_table_combo = ak.mask(tau_table_combo,tau_table_combo.idx1<tau_table_combo.idx2)
    # ensure particle anti-particle pair (charge=0)
    tau_table_combo = ak.mask(tau_table_combo,tau_table_combo.pair_sum.charge==0)
    tau_optimal_mass_diff = ak.min(ak.min(tau_table_combo.mass_diff, axis=-1),axis=-1)
    tau_min_diff_mask = tau_optimal_mass_diff == tau_table_combo.mass_diff
    final_tau_table = ak.drop_none(ak.mask(tau_table_combo, tau_min_diff_mask),axis=-1)

    # final tautau pairing
    tautau = ak.flatten(final_tau_table,axis=-1)[:,0]

    # final bb pairings
    bb1 = final_b_jet_table[:,0]
    bb2 = final_b_jet_table[:,1]

    # reconstructed higgs h1 and h2 are pt sorted from b-jets and h3 from taus
    h1 = bb1.pair_sum *1
    h2 = bb2.pair_sum *1
    h3 = tautau.pair_sum *1

    # unsorted h into bb, with h1_unsort with closest mass to 125
    h1_unsort = final_b_jet_table[final_b_jet_table.mass_diff == b_optimal_mass_diff1].pair_sum
    h2_unsort = final_b_jet_table[final_b_jet_table.mass_diff == b_optimal_mass_diff2].pair_sum
    
    # from IPython import embed; embed()

    events = set_ak_column_f32(
        events,
        "mhhh",
        (h1 + h2 + h3).mass,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_bb1",
        bb1.delta_r,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb1",
        bb1.cos,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_bb2",
        bb2.delta_r,
    )

    events = set_ak_column_f32(
        events,
        "cos_bb2",
        bb2.cos,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_tautau",
        tautau.delta_r,
    )

    events = set_ak_column_f32(
        events,
        "cos_tautau",
        tautau.cos,
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h12",
        h1.delta_r(h2),
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h13",
        h1.delta_r(h3),
    )

    events = set_ak_column_f32(
        events,
        "delta_r_h23",
        h2.delta_r(h3),
    )

    events = set_ak_column_f32(
        events,
        "cos_h12",
        h1.pvec.dot(h2.pvec)/(h1.pvec.absolute()*h2.pvec.absolute()),
    )

    events = set_ak_column_f32(
        events,
        "cos_h13",
        h1.pvec.dot(h3.pvec)/(h1.pvec.absolute()*h3.pvec.absolute()),
    )

    events = set_ak_column_f32(
        events,
        "cos_h23",
        h2.pvec.dot(h3.pvec)/(h2.pvec.absolute()*h3.pvec.absolute()),
    )

    events = set_ak_column_f32(
        events,
        "h1_mass",
        h1.mass,
    )

    events = set_ak_column_f32(
        events,
        "h2_mass",
        h2.mass,
    )

    events = set_ak_column_f32(
        events,
        "h3_mass",
        h3.mass,
    )

    events = set_ak_column(
        events,
        "n_b_jet",
        ak.num(b_jet,axis=-1),
        value_type=np.int32
    )

    events = set_ak_column(
        events, 
        "n_fatjet",
        ak.num(events.FatJet.pt,axis=-1),
        value_type=np.int32
        )
    
    events = set_ak_column_f32(
        events,
        "h1_unsort_mass",
        h1_unsort.mass,
    )

    events = set_ak_column_f32(
        events,
        "h2_unsort_mass",
        h2_unsort.mass,
    )

    return events