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

    # add_variable(
    #     config,
    #     name="mhhh",
    #     binning=(60, 150.0, 1300.0),
    #     unit="GeV",
    #     x_title=r"$m_{4b2\tau}$",
    # )

    # add_variable(
    #     config,
    #     name="cos_taulep",
    #     binning=(24, -1, +1),
    #     x_title=r"$\tau\tau$ $cos(\delta)$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_taulep",
    #     binning=(35, 0, 7),
    #     x_title=r"$\tau\tau$ $\Delta R$",
    # )

    # add_variable(
    #     config,
    #     name="h3_mass",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H3}$",
    # )

    # add_variable(
    #     config,
    #     name="n_fatjet",
    #     expression="n_fatjet",
    #     binning=(5, 0, 5),
    #     x_title="Number of fat jets",
    #     discrete_x=True,
    # )

    # add_variable(
    #     config,
    #     name="h1_unsort_mass",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H1}^{unsorted}$",
    # )

    # add_variable(
    #     config,
    #     name="h2_unsort_mass",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H2}^{unsorted}$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_h12",
    #     binning=(35, 0, 7),
    #     x_title=r"H $\Delta R_{1,2}$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_h13",
    #     binning=(35, 0, 7),
    #     x_title=r"H $\Delta R_{1,3}$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_h23",
    #     binning=(35, 0, 7),
    #     x_title=r"H $\Delta R_{2,3}$",
    # )
    
    # add_variable(
    #     config,
    #     name="delta_r_bb1",
    #     binning=(35, 0, 7),
    #     x_title=r"$bb_1$ $\Delta R$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_bb2",
    #     binning=(35, 0, 7),
    #     x_title=r"$bb_2$ $\Delta R$",
    # )

    # add_variable(
    #     config,
    #     name="cos_h12",
    #     binning=(24, -1, +1),
    #     x_title=r"H $cos(\delta)_{1,2}$",
    # )

    # add_variable(
    #     config,
    #     name="cos_h13",
    #     binning=(24, -1, +1),
    #     x_title=r"H $cos(\delta)_{1,3}$",
    # )

    # add_variable(
    #     config,
    #     name="cos_h23",
    #     binning=(24, -1, +1),
    #     x_title=r"H $cos(\delta)_{2,3}$",
    # )
    
    # add_variable(
    #     config,
    #     name="cos_bb1",
    #     binning=(24, -1, +1),
    #     x_title=r"$bb_1$ $cos(\delta)$",
    # )

    # add_variable(
    #     config,
    #     name="cos_bb2",
    #     binning=(24, -1, +1),
    #     x_title=r"$bb_2$ $cos(\delta)$",
    # )

    # add_variable(
    #     config,
    #     name="h1_mass",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H1}$",
    # )

    # add_variable(
    #     config,
    #     name="h2_mass",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H2}$",
    # )

    # add_variable(
    #     config,
    #     name="m_3btaulep",
    #     binning=(60, 150.0, 1300.0),
    #     unit="GeV",
    #     x_title=r"$m_{3b2\tau}, (b_3,hhbtag)$",
    # )

    # add_variable(
    #     config,
    #     name="m_3btaulep_pt",
    #     binning=(60, 150.0, 1300.0),
    #     unit="GeV",
    #     x_title=r"$m_{3b2\tau}, (b_3,pt)$",
    # )


    # gen-level variables
    add_variable(
        config,
        name="mtautau_gen",
        binning=(40, 120.0, 130.0),
        unit="GeV",
        x_title=r"$m_{\tau\tau}^{gen}$",
    )

    add_variable(
        config,
        name="mbb_gen",
        binning=(40, 120.0, 130.0),
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

    # add_variable(
    #     config,
    #     name="mds_h1_mass_chi",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H1}$ $(\Delta\chi^2)$ (mds)",
    # )

    # add_variable(
    #     config,
    #     name="mds_h2_mass_chi",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H2}$ $(\Delta\chi^2)$ (mds)",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_h12_chi",
    #     binning=(35, 0, 7),
    #     x_title=r"H $(\Delta\chi^2)$ $\Delta R_{1,2}$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_h13_chi",
    #     binning=(35, 0, 7),
    #     x_title=r"H $(\Delta\chi^2)$ $\Delta R_{1,3}$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_h23_chi",
    #     binning=(35, 0, 7),
    #     x_title=r"H $(\Delta\chi^2)$ $\Delta R_{2,3}$",
    # )
    
    # add_variable(
    #     config,
    #     name="delta_r_bb1_chi",
    #     binning=(35, 0, 7),
    #     x_title=r"$bb_1$ $(\Delta\chi^2)$ $\Delta R$",
    # )

    # add_variable(
    #     config,
    #     name="delta_r_bb2_chi",
    #     binning=(35, 0, 7),
    #     x_title=r"$bb_2$ $(\Delta\chi^2)$ $\Delta R$",
    # )

    # add_variable(
    #     config,
    #     name="cos_h12_chi",
    #     binning=(24, -1, +1),
    #     x_title=r"H $(\Delta\chi^2)$ $cos(\delta)_{1,2}$",
    # )

    # add_variable(
    #     config,
    #     name="cos_h13_chi",
    #     binning=(24, -1, +1),
    #     x_title=r"H $(\Delta\chi^2)$ $cos(\delta)_{1,3}$",
    # )

    # add_variable(
    #     config,
    #     name="cos_h23_chi",
    #     binning=(24, -1, +1),
    #     x_title=r"H $(\Delta\chi^2)$ $cos(\delta)_{2,3}$",
    # )
    
    # add_variable(
    #     config,
    #     name="cos_bb1_chi",
    #     binning=(24, -1, +1),
    #     x_title=r"$bb_1$ $(\Delta\chi^2)$ $cos(\delta)$",
    # )

    # add_variable(
    #     config,
    #     name="cos_bb2_chi",
    #     binning=(24, -1, +1),
    #     x_title=r"$bb_2$ $(\Delta\chi^2)$ $cos(\delta)$",
    # )

    # add_variable(
    #     config,
    #     name="h1_mass_chi",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H1}$ $(\Delta\chi^2)$",
    # )

    # add_variable(
    #     config,
    #     name="h2_mass_chi",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H2}$ $(\Delta\chi^2)$",
    # )

    # add_variable(
    #     config,
    #     name="m_3btaulep_chi",
    #     binning=(60, 150.0, 1300.0),
    #     unit="GeV",
    #     x_title=r"$m_{3b2\tau} (\Delta\chi^2), (b_3,hhbtag)$",
    # )

    # add_variable(
    #     config,
    #     name="m_3btaulep_pt_chi",
    #     binning=(60, 150.0, 1300.0),
    #     unit="GeV",
    #     x_title=r"$m_{3b2\tau} (\Delta\chi^2), (b_3,pt)$",
    # )


    # add_variable(
    #     config,
    #     name="mds_h1_mass_gm",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H1}^{mds,gm}$",
    # )

    # add_variable(
    #     config,
    #     name="mds_h2_mass_gm",
    #     binning=(40, 0.0, 400.0),
    #     unit="GeV",
    #     x_title=r"$m_{H2}^{mds,gm}$",
    # )

    def delta_r12(vectors):
        # delta r between first two elements
        # from IPython import embed; embed(header="delta_r12")
        dr = ak.firsts(vectors[:, :1], axis=1).delta_r(ak.firsts(vectors[:, 1:2], axis=1))
        return ak.fill_none(dr, EMPTY_FLOAT)
    
    def delta_r12_alt(array0, array1):
        dr = ak.firsts(array0, axis=1).delta_r(ak.firsts(array1, axis=1))
        return ak.fill_none(dr, EMPTY_FLOAT)
    

    # cosine of 3D angle between two particles
    def cos(vectors) -> ak.Array:
        cos = vectors[:, :1].pvec.dot(vectors[:, 1:2].pvec)/(vectors[:, :1].pvec.absolute()*vectors[:, 1:2].pvec.absolute())
        cos = ak.firsts(cos)
        return ak.fill_none(cos, EMPTY_FLOAT)
    
    def cos_alt(array0, array1) -> ak.Array:
        array0 = ak.firsts(array0)
        array1 = ak.firsts(array1)
        cos = array0.pvec.dot(array1.pvec)/(array0.pvec.absolute()*array1.pvec.absolute())
        return ak.fill_none(cos, EMPTY_FLOAT)
    
    def lvec_sum(vectors) -> ak.Array:
        return vectors[:, :1] + vectors[:, 1:2]

    def build_higgs_reco(events, obj = None, var = None):
        events = attach_coffea_behavior(events)

        bb1 = events.Jet[events.BB1_idx] *1
        bb2 = events.Jet[events.BB2_idx] *1
        leps = ak.concatenate([events.Electron * 1, events.Muon * 1, events.Tau * 1], axis=1)[:, :2]
        
        h1 = lvec_sum(bb1)
        h2 = lvec_sum(bb2)
        h3 = lvec_sum(leps)
        # from IPython import embed; embed(header="build_higgs_reco")
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

        # special treatment for these combinations because of emerging of union types by concatenating
        if obj == "h12":
            if  var == "cos":
                return cos_alt(h1,h2)
            if var == "dr":
                return delta_r12_alt(h1,h2)

        if obj == "h13":
            if  var == "cos":
                return cos_alt(h1,h3)
            if var == "dr":
                return delta_r12_alt(h1,h3)
            
        if obj == "h23":
            if  var == "cos":
                return cos_alt(h2,h3)
            if var == "dr":
                return delta_r12_alt(h2,h3)
        

        if obj == "3b2tau":
            mds_mask = (h1.mass-125)**2 <= (h2.mass-125)**2
            h_best = ak.where(mds_mask, h1, h2) *1
            h_best_idx = ak.where(mds_mask, events.BB1_idx, events.BB2_idx)
            rem_jet_mask = ((ak.local_index(events.Jet, axis=-1) != h_best_idx[:,0]) &
                            (ak.local_index(events.Jet, axis=-1) != h_best_idx[:,1]) )
            j3 = ak.unflatten(events.Jet[rem_jet_mask][:,0],1) *1
            vectors = h_best + h3 + j3

        if obj == "hhh":
            vectors = h1 + h2 + h3

        if  var == "cos":
            return cos(vectors)
        if var == "dr":
            return delta_r12(vectors)
        if var == "mass":
            return vectors.mass
        if var == "pt":
            return vectors.pt
        if var == "eta":
            return vectors.eta
        if var == "abs_eta":
            return abs(vectors.eta)
        if var == "phi":
            return vectors.phi
        if var == "energy":
            return vectors.energy
        
        
        raise ValueError(f"Unknown obj or var: {obj, var}")
    
    build_higgs_reco.inputs = ["{Electron,Muon,Tau,Jet}.{pt,eta,phi,mass}", "BB{1,2}_idx"]

    add_variable(
        config,
        name="min_chisq",
        binning=(50, 0, 0.4),
        x_title=r"minimal $\Delta\chi^2$",
    )

    add_variable(
        config,
        name="m_3btaulep_pt",
        expression=partial(build_higgs_reco, obj="3b2tau", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{3b2\tau}, (b_3,pt)$",
    )

    add_variable(
        config,
        name="mhhh",
        expression=partial(build_higgs_reco, obj="hhh", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(60, 150.0, 1300.0),
        unit="GeV",
        x_title=r"$m_{4b2\tau}$",
    )

    add_variable(
        config,
        name="delta_r_bb1",
        expression=partial(build_higgs_reco, obj="bb1", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $\Delta R$",
    )

    add_variable(
        config,
        name="delta_r_bb2",
        expression=partial(build_higgs_reco, obj="bb2", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $\Delta R$",
    )

    add_variable(
        config,
        name="delta_r_taulep",
        expression=partial(build_higgs_reco, obj="leps", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$\tau\tau$ $\Delta R$",
    )

    add_variable(
        config,
        name="cos_bb1",
        expression=partial(build_higgs_reco, obj="bb1", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="cos_bb2",
        expression=partial(build_higgs_reco, obj="bb2", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="cos_taulep",
        expression=partial(build_higgs_reco, obj="leps", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $cos(\delta)$",
    )

    add_variable(
        config,
        name="delta_r_h12",
        expression=partial(build_higgs_reco, obj="h12", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$H$ $\Delta R_{1,2}$",
    )

    add_variable(
        config,
        name="delta_r_h13",
        expression=partial(build_higgs_reco, obj="h13", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$H$ $\Delta R_{1,3}$",
    )

    add_variable(
        config,
        name="delta_r_h23",
        expression=partial(build_higgs_reco, obj="h23", var="dr"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(35, 0, 7),
        x_title=r"$H$ $\Delta R_{2,3}$",
    )

    add_variable(
        config,
        name="cos_h12",
        expression=partial(build_higgs_reco, obj="h12", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$H$ $cos(\delta)_{1,2}$",
    )

    add_variable(
        config,
        name="cos_h13",
        expression=partial(build_higgs_reco, obj="h13", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$H$ $cos(\delta)_{1,3}$",
    )

    add_variable(
        config,
        name="cos_h23",
        expression=partial(build_higgs_reco, obj="h23", var="cos"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(24, -1, +1),
        x_title=r"$H$ $cos(\delta)_{2,3}$",
    )

    add_variable(
        config,
        name="h1_mass",
        expression=partial(build_higgs_reco, obj="h1", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H1}$",
    )

    add_variable(
        config,
        name="h2_mass",
        expression=partial(build_higgs_reco, obj="h2", var="mass"),
        aux={"inputs": build_higgs_reco.inputs},
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$m_{H2}$",
    )

    add_variable(
        config,
        name="h3_mass",
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
        binning=(30, -3.0, 3.0),
        x_title=r"$H_3$ $|\eta|$",
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

