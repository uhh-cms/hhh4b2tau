import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, set_ak_column
from columnflow.production.util import attach_coffea_behavior
from columnflow.columnar_util import (
    set_ak_column, remove_ak_column, attach_behavior, EMPTY_FLOAT, get_ak_routes, remove_ak_column,
    optional_column as optional
)


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
        "jet_delta_r",
        "jet_delta_r13",
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
        "jet_delta_r",
        ak.where(dijet_mask, ak.mask(hardest_delta_r, dijet_mask)[:, 1], EMPTY_FLOAT),
    )

    # get delta r values between hardest and 3rd hardest jet
    events = set_ak_column_f32(
        events,
        "jet_delta_r13",
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
        optional(f"gen_{child}_from_{mother}.{var}")
        for mother in ('h', 'tau', )
        for child in ('b', 'tau', 'taunu', 'electron', 'enu', 'muon', 'munu', )
        for var in ('pt', 'eta', 'phi', 'mass', 'pdgId', )} | 
        {
            attach_coffea_behavior,
        }
    ),
    produces={
        "mtautau", "mbb", "mhhh", "mlnu", "hpt", "h1bpt", "h2bpt", "htaupt",
    },
)
def h_decay_invariant_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Construct b and tau invariant mass.
    """


    # attach coffea behavior for four-vector arithmetic
    events = self[attach_coffea_behavior](
        events,
        collections={ x : {
                "type_name": "GenParticle",
            } for x in [
            "gen_h_to_b", 
            "gen_b_from_h", 
            "gen_h_to_tau", 
            "gen_tau_from_h",
            "gen_nu_from_tau",
            "gen_taunu_from_tau",
            "gen_electron_from_tau",
            "gen_enu_from_tau",
            "gen_munu_from_tau",
            "gen_muon_from_tau"
            ]},
        **kwargs,
    )

    # from IPython import embed; embed()

    # total number of taus per event
    n_taus = ak.num(events.gen_tau_from_h, axis=-1)
    n_bs = ak.num(events.gen_b_from_h, axis=-1)
    n_hs = (ak.num(events.gen_h_to_tau, axis=-1) +
    ak.num(events.gen_h_to_b, axis=-1))
  
    # total number of leptons, neutrinos per event
    n_taunu = ak.num(events.gen_taunu_from_tau, axis=1)
    n_electron = ak.num(events.gen_electron_from_tau, axis=1)
    n_muon = ak.num(events.gen_muon_from_tau, axis=1)
    n_munu = ak.num(events.gen_munu_from_tau, axis=1)
    n_enu = ak.num(events.gen_enu_from_tau, axis=1)
    n_lep = n_electron + n_muon   

    # four-vector sum of first four elements of each
    # tau collection (possibly fewer)
    
    ditau = events.gen_tau_from_h.sum(axis=-1)
    dib = events.gen_b_from_h.sum(axis=-1)
    trih = events.gen_h_to_tau.sum(axis=-1) + events.gen_h_to_b.sum(axis=-1)
    # create 0 lorenz-vector
    zero_lorenz = ditau[0]-ditau[0]
    # lorenz vector sum of tau decay products
    ditaunu = ak.where(n_taunu>=1, ak.flatten(events.gen_taunu_from_tau.sum(axis=1)), zero_lorenz)
    dielectron = ak.fill_none(ak.pad_none(events.gen_electron_from_tau.sum(axis=1), 1), zero_lorenz[0])
    dienu = ak.fill_none(ak.pad_none(events.gen_enu_from_tau.sum(axis=1), 1), zero_lorenz[0])
    dimuon = ak.fill_none(ak.pad_none(events.gen_muon_from_tau.sum(axis=1), 1), zero_lorenz[0])
    dimunu = ak.fill_none(ak.pad_none(events.gen_munu_from_tau.sum(axis=1), 1), zero_lorenz[0])

    # lep_sum = dielectron + dimuon
    # nu_sum = ditaunu + dienu + dimunu
    lep_nu_sum = ditaunu + dienu + dimunu + dielectron + dimuon

    # pt of all Higgs
    h_pt = ak.concatenate(
        (events.gen_h_to_b.pt, events.gen_h_to_tau.pt), axis=1)
    h1b_pt = events.gen_h_to_b.pt[:,0]
    h2b_pt = events.gen_h_to_b.pt[:,1]
    htau_pt = events.gen_h_to_tau.pt

    # from IPython import embed; embed()
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

    lep_nu_from_tau_mass = ak.where(
        n_lep >= 2,
        lep_nu_sum.mass,
        [[EMPTY_FLOAT]],
    )

    # write out the resulting mass to the `events` array,
    events = set_ak_column_f32(
        events,
        "mtautau",
        tautau_mass,
    )

    events = set_ak_column_f32(
        events,
        "mbb",
        b_mass,
    )

    events = set_ak_column_f32(
        events,
        "mhhh",
        h_mass,
    )

    events = set_ak_column_f32(
        events,
        "mlnu",
        lep_nu_from_tau_mass,
    )
    events = set_ak_column_f32(
        events,
        "hpt",
        h_pt,
    )

    events = set_ak_column_f32(
        events,
        "htaupt",
        htau_pt,
    )

    events = set_ak_column_f32(
        events,
        "h1bpt",
        h1b_pt,
    )

    events = set_ak_column_f32(
        events,
        "h2bpt",
        h2b_pt,
    )
    # return the events
    # from IPython import embed; embed()
    return events
