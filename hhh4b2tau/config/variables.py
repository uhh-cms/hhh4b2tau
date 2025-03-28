# coding: utf-8

"""
Definition of variables.
"""
from functools import partial

import order as od

from columnflow.columnar_util import EMPTY_FLOAT, attach_coffea_behavior, default_coffea_collections
from columnflow.util import maybe_import

ak = maybe_import("awkward")

def add_variables(config: od.Config) -> None:
    # add variables
    # (the "event", "run" and "lumi" variables are required for some cutflow plotting task,
    # and also correspond to the minimal set of columns that coffea's nano scheme requires)
    add_variable(
        config,
        name="event",
        expression="event",
        binning=(1, 0.0, 1.0e6),
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
        binning=(30, -3.14, 3.14),
        x_title=r"Jets $\phi$",
    )
    add_variable(
        config,
        name="jet_delta_phi",
        binning=(30, -3.14, 3.14),
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
        x_title=r"$bb_1$ $cos(\delta)^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_bb2_hadron",
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $cos(\delta)^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_tautau_hadron",
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $cos(\delta)^{gen,hadron}$",
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
        x_title=r"H $cos(\delta)_{1,2}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_h13_hadron",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,3}^{gen,hadron}$",
    )

    add_variable(
        config,
        name="cos_h23_hadron",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{2,3}^{gen,hadron}$",
    )


    # detector level

    add_variable(
        config,
        name="mhhh",
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{4b2\tau}$",
    )

    add_variable(
        config,
        name="cos_taulep",
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="delta_r_taulep",
        binning=(35, 0, 7),
        x_title=r"$\tau\tau$ $\Delta R$",
    )

    add_variable(
        config,
        name="h3_mass",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H3}$",
    )

    add_variable(
        config,
        name="n_fatjet",
        expression="n_fatjet",
        binning=(5, 0, 5),
        x_title="Number of fat jets",
        discrete_x=True,
    )

    add_variable(
        config,
        name="h1_unsort_mass",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}^{unsorted}$",
    )

    add_variable(
        config,
        name="h2_unsort_mass",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}^{unsorted}$",
    )

    add_variable(
        config,
        name="delta_r_h12",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,2}$",
    )

    add_variable(
        config,
        name="delta_r_h13",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,3}$",
    )

    add_variable(
        config,
        name="delta_r_h23",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{2,3}$",
    )
    
    add_variable(
        config,
        name="delta_r_bb1",
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $\Delta R$",
    )

    add_variable(
        config,
        name="delta_r_bb2",
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $\Delta R$",
    )

    add_variable(
        config,
        name="cos_h12",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,2}$",
    )

    add_variable(
        config,
        name="cos_h13",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,3}$",
    )

    add_variable(
        config,
        name="cos_h23",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{2,3}$",
    )
    
    add_variable(
        config,
        name="cos_bb1",
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="cos_bb2",
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="h1_mass",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}$",
    )

    add_variable(
        config,
        name="h2_mass",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}$",
    )

    add_variable(
        config,
        name="m_3btaulep",
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{3b2\tau}, (b_3,hhbtag)$",
    )

    add_variable(
        config,
        name="m_3btaulep_pt",
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{3b2\tau}, (b_3,pt)$",
    )


    # gen-level variables
    add_variable(
        config,
        name="mtautau_gen",
        binning=(40, 100.0, 150.0),
        unit="GeV",
        x_title=r"$m_{\tau\tau}^{gen}$",
    )

    add_variable(
        config,
        name="mbb_gen",
        binning=(40, 100.0, 150.0),
        unit="GeV",
        x_title=r"$m_{bb}^{gen}$",
    )

    add_variable(
        config,
        name="mhhh_gen",
        binning=(60, 150.0, 1300.0),
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
        x_title=r"$p_{TH}^{gen}$",
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
        x_title=r"H $\Delta R_{1,2}^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_h13_gen",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,3}^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_h23_gen",
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{2,3}^{gen}$",
    )
    
    add_variable(
        config,
        name="delta_r_bb1_gen",
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $\Delta R^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_bb2_gen",
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $\Delta R^{gen}$",
    )

    add_variable(
        config,
        name="delta_r_tautau_gen",
        binning=(35, 0, 7),
        x_title=r"$\tau\tau$ $\Delta R^{gen}$",
    )

    add_variable(
        config,
        name="cos_h12_gen",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,2}^{gen}$",
    )

    add_variable(
        config,
        name="cos_h13_gen",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,3}^{gen}$",
    )

    add_variable(
        config,
        name="cos_h23_gen",
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{2,3}^{gen}$",
    )
    
    add_variable(
        config,
        name="cos_bb1_gen",
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $cos(\delta)^{gen}$",
    )

    add_variable(
        config,
        name="cos_bb2_gen",
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $cos(\delta)^{gen}$",
    )

    add_variable(
        config,
        name="cos_tautau_gen",
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $cos(\delta)^{gen}$",
    )


    ### detector level but with experimental Delta chi**2 minimization to H mass for jet pairing

    add_variable(
        config,
        name="min_chi",
        binning=(50, 0, 0.2),
        x_title=r"minimal $\Delta\chi^2$",
    )

    add_variable(
        config,
        name="mds_h1_mass_chi",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}$ $(\Delta\chi^2)$ (mds)",
    )

    add_variable(
        config,
        name="mds_h2_mass_chi",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}$ $(\Delta\chi^2)$ (mds)",
    )

    add_variable(
        config,
        name="delta_r_h12_chi",
        binning=(35, 0, 7),
        x_title=r"H $(\Delta\chi^2)$ $\Delta R_{1,2}$",
    )

    add_variable(
        config,
        name="delta_r_h13_chi",
        binning=(35, 0, 7),
        x_title=r"H $(\Delta\chi^2)$ $\Delta R_{1,3}$",
    )

    add_variable(
        config,
        name="delta_r_h23_chi",
        binning=(35, 0, 7),
        x_title=r"H $(\Delta\chi^2)$ $\Delta R_{2,3}$",
    )
    
    add_variable(
        config,
        name="delta_r_bb1_chi",
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $(\Delta\chi^2)$ $\Delta R$",
    )

    add_variable(
        config,
        name="delta_r_bb2_chi",
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $(\Delta\chi^2)$ $\Delta R$",
    )

    add_variable(
        config,
        name="cos_h12_chi",
        binning=(24, -1, +1),
        x_title=r"H $(\Delta\chi^2)$ $cos(\delta)_{1,2}$",
    )

    add_variable(
        config,
        name="cos_h13_chi",
        binning=(24, -1, +1),
        x_title=r"H $(\Delta\chi^2)$ $cos(\delta)_{1,3}$",
    )

    add_variable(
        config,
        name="cos_h23_chi",
        binning=(24, -1, +1),
        x_title=r"H $(\Delta\chi^2)$ $cos(\delta)_{2,3}$",
    )
    
    add_variable(
        config,
        name="cos_bb1_chi",
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $(\Delta\chi^2)$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="cos_bb2_chi",
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $(\Delta\chi^2)$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="h1_mass_chi",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}$ $(\Delta\chi^2)$",
    )

    add_variable(
        config,
        name="h2_mass_chi",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}$ $(\Delta\chi^2)$",
    )

    add_variable(
        config,
        name="m_3btaulep_chi",
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{3b2\tau} (\Delta\chi^2), (b_3,hhbtag)$",
    )

    add_variable(
        config,
        name="m_3btaulep_pt_chi",
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{3b2\tau} (\Delta\chi^2), (b_3,pt)$",
    )


    add_variable(
        config,
        name="mds_h1_mass_gm",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}^{mds,gm}$",
    )

    add_variable(
        config,
        name="mds_h2_mass_gm",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}^{mds,gm}$",
    )

    # build variables for dilepton, and jet_lepton
    def delta_r12(vectors):
        # delta r between first two elements
        dr = ak.firsts(vectors[:, :1], axis=1).delta_r(ak.firsts(vectors[:, 1:2], axis=1))
        return ak.fill_none(dr, EMPTY_FLOAT)

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
        x_title=r"$\Delta R_{l,jet}$",
    )

    add_variable(
        config,
        name="lep_jet_mass",
        expression=partial(build_dilep, which="lep_jet_mass"),
        aux={"inputs": build_dilep.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{l,jet}$",
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
    if not variable.has_aux("underflow"):
        variable.x.underflow = True
    if not variable.has_aux("overflow"):
        variable.x.overflow = True

    return variable
