# coding: utf-8

"""
Definition of variables.
"""
from functools import partial

import order as od

from columnflow.columnar_util import EMPTY_FLOAT, attach_coffea_behavior, default_coffea_collections
from columnflow.util import maybe_import

ak = maybe_import("awkward")
np = maybe_import("numpy")

def delta_r12(*vectors):
        # delta r between first two elements
        if len(vectors[-1]) == 2:
            vec0 = vectors[0][0]
            vec1 = vectors[0][1]
        else:
            vec0 = ak.firsts(vectors[0][:, :1])
            vec1 = ak.firsts(vectors[0][:, 1:2])
        dr = vec0.delta_r(vec1)
        return dr

def cos(*vectors) -> ak.Array:
    if len(vectors[-1]) == 2:
        vec0 = vectors[0][0]
        vec1 = vectors[0][1]
    else:
        vec0 = ak.firsts(vectors[0][:, :1])
        vec1 = ak.firsts(vectors[0][:, 1:2])
    cos = vec0.pvec.dot(vec1.pvec)/(vec0.pvec.absolute()*vec1.pvec.absolute())
    return cos

def lvec_sum(*vectors) -> ak.Array:
    if len(vectors) == 1:
        sum = ak.firsts(vectors[0][:, :1]) + ak.firsts(vectors[0][:, 1:2])
    else: 
        sum = vectors[0]
        for vec in vectors[1:]:
            sum = sum + vec
    return sum


def output_vectors(
        obj=None, 
        bb1=None, 
        bb2=None, 
        leps=None, 
        h1=None, 
        h2=None, 
        h3=None, 
        jjj=None
    ):

    if obj == "bb1":
        vectors = bb1
    if obj == "bb2":
        vectors = bb2
    if obj == "leps":
        vectors = leps

    if obj == "h1":
        vectors = h1
    if obj == "h2":
        vectors = h2
    if obj == "h3":
        vectors = h3

    if obj == "h12":
        vectors = h1, h2
    if obj == "h13":
        vectors = h1, h3
    if obj == "h23":
        vectors = h2, h3

    if obj == "3j2l":
        vectors = (jjj + h3)
    if obj == "hhh":
        vectors = lvec_sum(h1, h2, h3)
        
    return vectors

# generic function to put out variable values
def output_value(vectors=None, var=None):

    if  var == "cos":
        value = cos(vectors)
    if var == "dr":
        value = delta_r12(vectors)
    if var == "mass":
        value = vectors.mass
    if var == "pt":
        value = vectors.pt
    if var == "eta":
        value = vectors.eta
    if var == "abs_eta":
        value = abs(vectors.eta)
    if var == "phi":
        value = vectors.phi
    if var == "energy":
        value = vectors.energy

    # if necessary reduce dimensions for categorizer
    try:
        value = ak.firsts(value)
    except:
        pass
    
    return ak.fill_none(value, EMPTY_FLOAT)


def build_higgs_reco(events, obj = None, var = None):
    events = attach_coffea_behavior(events)
    try:
        bb1 = events.Jet[events.BB1_idx] *1
        bb2 = events.Jet[events.BB2_idx] *1

        h1 = lvec_sum(bb1)
        h2 = lvec_sum(bb2)
    except:
        pass

    leps = ak.concatenate([events.Electron * 1, events.Muon * 1, events.Tau * 1], axis=1)[:, :2]

    h3 = lvec_sum(leps)

    btag_idx = ak.argsort(events.Jet.btagDeepFlavB, ascending=False)
    jjj = events.Jet[btag_idx][:, :3].sum(axis=1) *1

    vectors = output_vectors(
        obj=obj,
        bb1=bb1, 
        bb2=bb2,
        leps=leps,
        h1=h1,
        h2=h2,
        h3=h3,
        jjj=jjj,
    )

    return output_value(vectors=vectors,var=var)


build_higgs_reco.inputs = [
    "{Electron,Muon,Tau,Jet}.{pt,eta,phi,mass}", 
    "BB{1,2}_idx", 
    "Jet.btagDeepFlavB"
]

def build_jet_gen_matched(events, obj=None, var = None, mds=False):
    
    events = attach_coffea_behavior(events)
    bb1 = events.Jet[events.Gen_Matched_H1_idx] *1
    bb2 = events.Jet[events.Gen_Matched_H2_idx] *1
    if mds == False:
        from hhh4b2tau.production.util import swap_random
        bb1, bb2 = swap_random(bb1, bb2)

    h1 = lvec_sum(bb1)
    h2 = lvec_sum(bb2)
    # mass difference sorting: mass of H1 is closer to 125 GeV 
    if mds == True: 
        from hhh4b2tau.production.util import order_pairs
        mds_mask = ((h1.mass - 125)**2 <= (h2.mass - 125)**2)
        # if only one H is reconstructed it will be h1
        mds_mask = ak.where(ak.is_none(h1), False, mds_mask)
        mds_mask = ak.fill_none(mds_mask, True)
        h1, h2 = order_pairs(h1, h2, mds_mask)


    vectors = output_vectors(
        obj=obj,
        # bb1=bb1, 
        # bb2=bb2,
        # leps=leps,
        h1=h1,
        h2=h2,
        # h3=h3,
        # jjj=jjj,
    )

    return output_value(vectors=vectors,var=var)

build_jet_gen_matched.inputs = [
    "Jet.{pt,eta,phi,mass,btagDeepFlavB}",
    "Gen_Matched_H{1,2}_idx",
]

# gen level objects
def build_gen_higgs(events, obj=None, var=None):
    events = attach_coffea_behavior(
        events,
        collections={ x : {"type_name": "GenParticle"} for x in build_gen_higgs.objects}
    )

    h1 = events.gen_h_to_b[:,0] *1
    h2 = events.gen_h_to_b[:,1] *1
    h3 = ak.firsts(events.gen_h_to_tau) *1

    bb1 = events.gen_b[:, 0] * 1
    bb2 = events.gen_b[:, 1] * 1
    leps = ak.firsts(events.gen_tau) *1

    bbbb = ak.flatten(events.gen_b,axis=2)
    jjj_idx = ak.argsort(bbbb.pt, ascending=False)[:,:3]
    jjj = bbbb[jjj_idx].sum(axis=1)*1

    vectors = output_vectors(
        obj=obj,
        bb1=bb1, 
        bb2=bb2,
        leps=leps,
        h1=h1,
        h2=h2,
        h3=h3,
        jjj=jjj,
    )
    
    return output_value(vectors=vectors,var=var)

build_gen_higgs.objects = [
    "gen_b", 
    "gen_tau", 
    "gen_h_to_b", 
    "gen_h_to_tau",
    # "gen_taunu",
    # "gen_electron",
    # "gen_enu",
    # "gen_munu",
    # "gen_muon",
]

build_gen_higgs.variables = [
    'pt', 
    'eta', 
    'phi', 
    'mass', 
    'pdgId', 
]

build_gen_higgs.inputs = [
    f"{obj}.{var}" 
    for obj in build_gen_higgs.objects 
    for var in build_gen_higgs.variables
]  

# variables that use gen and detector level objects
# set op for a specific operation
def build_higgs_reco_gen(events, obj=None, var=None, op=None):
    events = attach_coffea_behavior(events)
    events = attach_coffea_behavior(
        events,
        collections={ x : {"type_name": "GenParticle"} for x in build_gen_higgs.objects}
    )
    reco_var = build_higgs_reco(events, obj=obj, var=var)
    gen_var = build_gen_higgs(events, obj=obj, var=var)
    # consider EMPTY_FLOAT values 
    ef_mask = (reco_var != EMPTY_FLOAT) & (gen_var != EMPTY_FLOAT)

    if op == "diff":
        value = gen_var - reco_var

    value = ak.where(ef_mask, value, EMPTY_FLOAT)

    return value

build_higgs_reco_gen.inputs = build_higgs_reco.inputs + build_gen_higgs.inputs



def add_variables(config: od.Config) -> None:
    # add variables
    # (the "event", "run" and "lumi" variables are required for some cutflow plotting task,
    # and also correspond to the minimal set of columns that coffea's nano scheme requires)
    add_variable(
        config,
        name="event",
        expression="event",
        binning=(1, 0.0, 1e9),
        x_title="Event number",
    )
    add_variable(
        config,
        name="run",
        expression="run",
        binning=(1, 100000.0, 500000.0),
        x_title="Run number",
        discrete_x=True,
    )
    add_variable(
        config,
        name="lumi",
        expression="luminosityBlock",
        binning=(1, 0.0, 5000.0),
        x_title="Luminosity block",
        discrete_x=True,
    )
    add_variable(
        config,
        name="n_jet",
        expression="n_jet",
        binning=(15, 0, 15),
        x_title="Number of jets",
        discrete_x=True,
    )
    # pt of all jets in every event
    add_variable(
        config,
        name="jets_pt",
        expression="Jet.pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$p_{T}$ of all jets",
    )
    # pt of the first jet in every event
    add_variable(
        config,
        name="jet1_pt",  # variable name, to be given to the "--variables" argument for the plotting task
        expression="Jet.pt[:,0]",  # content of the variable
          # value to be given if content not available for event
        binning=(50, 0.0, 500.0),  # (bins, lower edge, upper edge)
        unit="GeV",  # unit of the variable, if any
        x_title=r"Leading Jet $p_{T}$",  # x title of histogram when plotted
    )
    # eta of the first jet in every event
    add_variable(
        config,
        name="jet1_eta",
        expression="Jet.eta[:,0]",
        binning=(30, -3.0, 3.0),
        x_title=r"Leading Jet $\eta$",
    )

    def build_ht(events):
        objects = ak.concatenate([events.Electron * 1, events.Muon * 1, events.Tau * 1, events.Jet * 1], axis=1)[:, :]
        objects_sum = objects.sum(axis=1)
        return objects_sum.pt
    build_ht.inputs = ["{Electron,Muon,Tau,Jet}.{pt,eta,phi,mass}"]

    add_variable(
        config,
        name="ht",
        expression=partial(build_ht),
        aux={"inputs": build_ht.inputs},
        binning=[0, 80, 120, 160, 200, 240, 280, 320, 400, 500, 600, 800],
        unit="GeV",
        x_title="HT",
    )
    # weights
    add_variable(
        config,
        name="mc_weight",
        expression="mc_weight",
        binning=(200, -10, 10),
        x_title="MC weight",
    )
    # cutflow variables
    add_variable(
        config,
        name="cf_jet1_pt",
        expression="cutflow.jet1_pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"Leading Jet $p_{T}$",
    )
 # add new variables

    add_variable(
        config,
        name="jet_phi",
        expression="Jet.phi",
        binning=(30, -np.pi, np.pi),
        x_title=r"Jets $\phi$",
    )
    add_variable(
        config,
        name="jet_delta_phi",
        binning=(30, -np.pi, np.pi),
        x_title=r"Jets $\Delta\phi_{1,2}$",
    )

    add_variable(
        config,
        name="jet_delta_r_12",
        binning=(30, 0, 6),
        x_title=r"Jets $\Delta R_{1,2}$",
    )

    add_variable(
        config,
        name="jet_delta_r_13",
        binning=(30, 0, 6),
        x_title=r"Jets $\Delta R_{1,3}$",
    )


    add_variable(
        config,
        name="ele_eta",
        expression="Electron.eta[:,0]",
        binning=(30, -3.0, 3.0),
        x_title=r"Leading Electron $\eta$",
    )

    add_variable(
        config,
        name="electron_pt",
        expression="Electron.pt[:,0]",
        binning=(40, 0.0, 300.0),
        unit="GeV",
        x_title=r"Leading Electron $p_{T}$",
    )

    add_variable(
        config,
        name="muon_pt",
        expression="Muon.pt[:,0]",
        binning=(40, 0.0, 200.0),
        unit="GeV",
        x_title=r"Leading Muon $p_{T}$",
    )



    add_variable(
        config,
        name="mhhh_hadron",
        binning=(60, 0.0, 1500.0),
        unit="GeV",
        x_title=r"$m_{HHH}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="delta_r_bb1_hadron",
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $\Delta R^{gen,hadron}$",
    )

    add_variable(
        config,
        name="delta_r_bb2_hadron",
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $\Delta R^{gen,hadron}$",
    )

    add_variable(
        config,
        name="delta_r_tautau_hadron",
        binning=(35, 0, 7),
        x_title=r"$\tau\tau$ $\Delta R^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_bb1_hadron",
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $\cos(\delta)^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_bb2_hadron",
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $\cos(\delta)^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_tautau_hadron",
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $\cos(\delta)^{gen,hadron}$",
    )

    add_variable(
        config,
        name="delta_r_h12_hadron",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,2}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="delta_r_h13_hadron",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,3}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="delta_r_h23_hadron",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{2,3}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_h12_hadron",
        binning=(24, -1, +1),
        x_title=r"H $\cos(\delta)_{1,2}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_h13_hadron",
        binning=(24, -1, +1),
        x_title=r"H $\cos(\delta)_{1,3}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_h23_hadron",
        binning=(24, -1, +1),
        x_title=r"H $\cos(\delta)_{2,3}^{gen,hadron}$",
    )

    # gen-level variables
    add_variable(
        config,
        name="mtautau_gen",
        binning=(40, 123.0, 126.0),
        unit="GeV",
        x_title=r"$m_{H\rightarrow\tau\tau}^{gen}$",
    )

    add_variable(
        config,
        name="mbb_gen",
        binning=(40, 121.5, 126.0),
        unit="GeV",
        x_title=r"$m_{H \rightarrow bb}^{gen}$",
    )

    add_variable(
        config,
        name="mhhh_gen",
        binning=(60, 350.0, 1350.0),
        unit="GeV",
        x_title=r"$m_{HHH}^{gen}$",
    )

    add_variable(
        config,
        name="mlnu_gen",
        binning=(40, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{l\nu}^{gen}$",
    )

    add_variable(
        config,
        name="hpt_gen",
        binning=(60, 0.0, 800.0),
        unit="GeV",
        x_title=r"$p_{T,H}^{gen}$",
    )

    add_variable(
        config,
        name="h1bpt_gen",
        binning=(60, 0.0, 800.0),
        unit="GeV",
        x_title=r"$p_{T,H_1\rightarrow bb}^{gen}$",
    )

    add_variable(
        config,
        name="h2bpt_gen",
        binning=(60, 0.0, 800.0),
        unit="GeV",
        x_title=r"$p_{T,H_2\rightarrow bb}^{gen}$",
    )

    add_variable(
        config,
        name="htaupt_gen",
        binning=(60, 0.0, 800.0),
        unit="GeV",
        x_title=r"$p_{T,H\rightarrow\tau\tau}^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_h12_gen",
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H1,2}^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_h13_gen",
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H1,3}^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_h23_gen",
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H2,3}^{gen}$",
    )
    
    add_variable(
        config,
        name="delta_r_bb1_gen",
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H1\rightarrow bb}^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_bb2_gen",
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H2\rightarrow bb}^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_tautau_gen",
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H3\rightarrow \tau\tau}^{gen}$",
    )

    add_variable(
        config,
        name="cos_h12_gen",
        binning=(24, -1, +1),
        x_title=r"$cos(\delta_{H1,2})^{gen}$",
    )

    add_variable(
        config,
        name="cos_h13_gen",
        binning=(24, -1, +1),
        x_title=r"$cos(\delta_{H1,3})^{gen}$",
    )

    add_variable(
        config,
        name="cos_h23_gen",
        binning=(24, -1, +1),
        x_title=r"$cos(\delta_{H2,3})^{gen}$",
    )
    
    add_variable(
        config,
        name="cos_bb1_gen",
        binning=(24, -1, +1),
        x_title=r"$cos(\delta_{H1\rightarrow bb})^{gen}$",
    )

    add_variable(
        config,
        name="cos_bb2_gen",
        binning=(24, -1, +1),
        x_title=r"$cos(\delta_{H2\rightarrow bb})^{gen}$",
    )

    add_variable(
        config,
        name="cos_tautau_gen",
        binning=(24, -1, +1),
        x_title=r"$cos(\delta_{H3\rightarrow \tau\tau})^{gen}$",
    )

    add_variable(
        config,
        name="chi2",
        # binning=(50, 0, 0.4),
        # for flat-s binning
        binning=(5000, 0, 0.4),
        x_title=r"$\chi^2$",
    )

    add_variable(
        config,
        name="chi21",
        # binning=(50, 0, 0.15),
        # for flat-s binning
        binning=(5000, 0, 0.15),
        x_title=r"$\chi^2_{H1}$",
    )

    add_variable(
        config,
        name="chi22",
        # binning=(50, 0, 0.4),
        # for flat-s binning
        binning=(5000, 0, 0.4),
        x_title=r"$\chi^2_{H2}$",
    )

    add_variable(
        config,
        name="dhh",
        # binning=(50, 0, 20),
        # for flat-s binning
        binning=(5000, 0, 20),
        x_title=r"$D_{HH}$",
    )

    add_variable(
        config,
        name="m3j2l",
        expression=partial(build_higgs_reco, obj="3j2l", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        # binning=(60, 150.0, 1300.0),
        # for flat-s binning
        binning=(5000, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{3j2l}$",
    )

    add_variable(
        config,
        name="delta_m3j2l",
        expression=partial(
            build_higgs_reco_gen, 
            obj="3j2l", 
            var="mass", 
            op="diff"
            ),
        aux={"inputs": build_higgs_reco_gen.inputs},
        binning=(60, -100.0, 700.0),
        unit="GeV",
        x_title=r"$\Delta m_{3j2l}^{gen,det}$",
    )

    add_variable(
        config,
        name="delta_mhhh",
        expression=partial(
            build_higgs_reco_gen, 
            obj="hhh", 
            var="mass", 
            op="diff"
            ),
        aux={"inputs": build_higgs_reco_gen.inputs},
        binning=(60, -400.0, 700.0),
        unit="GeV",
        x_title=r"$\Delta m_{4j2l}^{gen,det}$",
    )

    add_variable(
        config,
        name="mhhh",
        expression=partial(build_higgs_reco, obj="hhh", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{4j2l}$",
    )

    add_variable(
        config,
        name="delta_r_jj1",
        expression=partial(build_higgs_reco, obj="bb1", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{jj1}$",
    )

    add_variable(
        config,
        name="delta_r_jj2",
        expression=partial(build_higgs_reco, obj="bb2", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{jj2}$",
    )

    add_variable(
        config,
        name="delta_r_ll",
        expression=partial(build_higgs_reco, obj="leps", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{ll}$",
    )

    add_variable(
        config,
        name="cos_jj1",
        expression=partial(build_higgs_reco, obj="bb1", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$\cos(\delta_{jj1})$",
    )

    add_variable(
        config,
        name="cos_jj2",
        expression=partial(build_higgs_reco, obj="bb2", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$\cos(\delta_{jj2})$",
    )

    add_variable(
        config,
        name="cos_ll",
        expression=partial(build_higgs_reco, obj="leps", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$\cos(\delta_{ll})$",
    )

    add_variable(
        config,
        name="delta_r_h12",
        expression=partial(build_higgs_reco, obj="h12", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H1,2}$",
    )

    add_variable(
        config,
        name="delta_r_h13",
        expression=partial(build_higgs_reco, obj="h13", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H1,3}$",
    )

    add_variable(
        config,
        name="delta_r_h23",
        expression=partial(build_higgs_reco, obj="h23", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{H2,3}$",
    )

    add_variable(
        config,
        name="cos_h12",
        expression=partial(build_higgs_reco, obj="h12", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$\cos(\delta)_{H1,2}$",
    )

    add_variable(
        config,
        name="cos_h13",
        expression=partial(build_higgs_reco, obj="h13", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$\cos(\delta)_{H1,3}$",
    )

    add_variable(
        config,
        name="cos_h23",
        expression=partial(build_higgs_reco, obj="h23", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$\cos(\delta)_{H2,3}$",
    )

    add_variable(
        config,
        name="mh1",
        expression=partial(build_higgs_reco, obj="h1", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}$",
    )

    add_variable(
        config,
        name="mh2",
        expression=partial(build_higgs_reco, obj="h2", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}$",
    )

    add_variable(
        config,
        name="mh3",
        expression=partial(build_higgs_reco, obj="h3", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H3}$",
    )

    add_variable(
        config,
        name="h1_energy",
        expression=partial(build_higgs_reco, obj="h1", var="energy"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$H_1$ energy",
    )

    add_variable(
        config,
        name="h2_energy",
        expression=partial(build_higgs_reco, obj="h2", var="energy"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$H_2$ energy",
    )

    add_variable(
        config,
        name="h3_energy",
        expression=partial(build_higgs_reco, obj="h3", var="energy"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$H_3$ energy",
    )

    add_variable(
        config,
        name="h1_pt",
        expression=partial(build_higgs_reco, obj="h1", var="pt"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$H_1$ $p_{T}$",
    )

    add_variable(
        config,
        name="h2_pt",
        expression=partial(build_higgs_reco, obj="h2", var="pt"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$H_2$ $p_{T}$",
    )

    add_variable(
        config,
        name="h3_pt",
        expression=partial(build_higgs_reco, obj="h3", var="pt"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$H_3$ $p_{T}$",
    )

    add_variable(
        config,
        name="h1_phi",
        expression=partial(build_higgs_reco, obj="h1", var="phi"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, -np.pi, np.pi),
        x_title=r"$H_1$ $\phi$",
    )

    add_variable(
        config,
        name="h2_phi",
        expression=partial(build_higgs_reco, obj="h2", var="phi"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, -np.pi, np.pi),
        x_title=r"$H_2$ $\phi$",
    )

    add_variable(
        config,
        name="h3_phi",
        expression=partial(build_higgs_reco, obj="h3", var="phi"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, -np.pi, np.pi),
        x_title=r"$H_3$ $\phi$",
    )

    add_variable(
        config,
        name="h1_eta",
        expression=partial(build_higgs_reco, obj="h1", var="eta"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, -3.0, 3.0),
        x_title=r"$H_1$ $\eta$",
    )

    add_variable(
        config,
        name="h2_eta",
        expression=partial(build_higgs_reco, obj="h2", var="eta"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, -3.0, 3.0),
        x_title=r"$H_2$ $\eta$",
    )

    add_variable(
        config,
        name="h3_eta",
        expression=partial(build_higgs_reco, obj="h3", var="eta"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, -3.0, 3.0),
        x_title=r"$H_3$ $\eta$",
    )

    add_variable(
        config,
        name="h1_abs_eta",
        expression=partial(build_higgs_reco, obj="h1", var="abs_eta"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, 0.0, 3.0),
        x_title=r"$H_1$ $|\eta|$",
    )

    add_variable(
        config,
        name="h2_abs_eta",
        expression=partial(build_higgs_reco, obj="h2", var="abs_eta"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, 0.0, 3.0),
        x_title=r"$H_2$ $|\eta|$",
    )

    add_variable(
        config,
        name="h3_abs_eta",
        expression=partial(build_higgs_reco, obj="h3", var="abs_eta"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(30, 0.0, 3.0),
        x_title=r"$H_3$ $|\eta|$",
    )

    add_variable(
        config,
        name="h1_mass_gm",
        expression=partial(build_jet_gen_matched, obj="h1", var="mass"),
        aux={"inputs": build_jet_gen_matched.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}^{gm}$",
    )

    add_variable(
        config,
        name="h2_mass_gm",
        expression=partial(build_jet_gen_matched, obj="h2", var="mass"),
        aux={"inputs": build_jet_gen_matched.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}^{gm}$",
    )

    add_variable(
        config,
        name="h1_mass_gm_mds",
        expression=partial(build_jet_gen_matched, obj="h1", var="mass", mds=True),
        aux={"inputs": build_jet_gen_matched.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}^{gm,mds}$",
    )

    add_variable(
        config,
        name="h2_mass_gm_mds",
        expression=partial(build_jet_gen_matched, obj="h2", var="mass", mds=True),
        aux={"inputs": build_jet_gen_matched.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}^{gm,mds}$",
    )


    # build variables for dilepton, and jet_lepton

    
    def build_dilep(events, which=None):
        events = attach_coffea_behavior(events)
        leps = ak.concatenate([events.Electron * 1, events.Muon * 1, events.Tau * 1], axis=1)[:, :2]
        if which == "dr":
            return delta_r12(leps)
        dilep = leps.sum(axis=1)
        if which is None:
            return dilep * 1
        if which == "mass":
            return dilep.mass
        if which == "pt":
            return dilep.pt
        if which == "eta":
            return dilep.eta
        if which == "abs_eta":
            return abs(dilep.eta)
        if which == "phi":
            return dilep.phi
        if which == "energy":
            return dilep.energy
        min_jet = events.Jet[ak.unflatten(ak.argmin(leps[:,0].delta_r(events.Jet),axis=1),1)] * 1
        if which == "jet_dr":
            return leps[:,0].delta_r(min_jet)
        if which == "lep_jet_mass":
            return (leps[:,0] + min_jet).mass
        raise ValueError(f"Unknown which: {which}")

    build_dilep.inputs = ["{Electron,Muon,Tau,Jet}.{pt,eta,phi,mass}"]

    add_variable(
        config,
        name="lep_jet_dr",
        expression=partial(build_dilep, which="jet_dr"),
        aux={"inputs": build_dilep.inputs},
        binning=(35, 0, 7),
        x_title=r"$\Delta R_{l,j}$",
    )

    add_variable(
        config,
        name="lep_jet_mass",
        expression=partial(build_dilep, which="lep_jet_mass"),
        aux={"inputs": build_dilep.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{l,j}$",
    )



    def build_dibjet(events, which=None):
        events = attach_coffea_behavior(events, {"HHBJet": default_coffea_collections["Jet"]})
        # from IPython import embed; embed(header="dibjet")
        hhbjets = events.HHBJet[:, :2]
        if which == "dr":
            return delta_r12(hhbjets)
        dijet = hhbjets.sum(axis=1)
        if which is None:
            return dijet * 1
        if which == "mass":
            return dijet.mass
        if which == "pt":
            return dijet.pt
        if which == "eta":
            return dijet.eta
        if which == "abs_eta":
            return abs(dijet.eta)
        if which == "phi":
            return dijet.phi
        if which == "energy":
            return dijet.energy
        raise ValueError(f"Unknown which: {which}")

    build_dibjet.inputs = ["HHBJet.{pt,eta,phi,mass}"]

    add_variable(
        config,
        name="hhbjet_mass",
        expression=partial(build_dibjet, which="mass"),
        aux={"inputs": build_dibjet.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{hhbjet}$",
    )


    def build_bjets(events, which=None, algo=None):
        wp = "medium"
        if algo == "btagPNetB":
            wp_value = config.x.btag_working_points["particleNet"][wp]
        elif algo == "btagDeepFlavB":
            wp_value = config.x.btag_working_points["deepjet"][wp]
        else:
            raise ValueError(f"Unknown which: {algo}")
        bjet_mask = events.Jet[algo] >= wp_value
        objects = events.Jet[bjet_mask]
        if which == "energy":
            return objects.energy
        raise ValueError(f"Unknown which: {which}")

    build_bjets.inputs = ["Jet.{btagPNetB,btagDeepFlavB}"]



        






# helper to add a variable to the config with some defaults
def add_variable(config: od.Config, *args, **kwargs) -> od.Variable:
    kwargs.setdefault("null_value", EMPTY_FLOAT)

    # create the variable
    variable = config.add_variable(*args, **kwargs)

    # defaults
    # if not variable.has_aux("underflow"):
    #     variable.x.underflow = True
    # if not variable.has_aux("overflow"):
    #     variable.x.overflow = True

    return variable

