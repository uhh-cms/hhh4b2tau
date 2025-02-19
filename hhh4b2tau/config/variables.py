from columnflow.columnar_util import EMPTY_FLOAT
from columnflow.util import maybe_import

ak = maybe_import("awkward")

def add_variables(
    cfg
):
    # add variables
    # (the "event", "run" and "lumi" variables are required for some cutflow plotting task,
    # and also correspond to the minimal set of columns that coffea's nano scheme requires)
    cfg.add_variable(
        name="event",
        expression="event",
        binning=(1, 0.0, 1.0e9),
        x_title="Event number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="run",
        expression="run",
        binning=(1, 100000.0, 500000.0),
        x_title="Run number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="lumi",
        expression="luminosityBlock",
        binning=(1, 0.0, 5000.0),
        x_title="Luminosity block",
        discrete_x=True,
    )
    cfg.add_variable(
        name="n_jet",
        expression="n_jet",
        binning=(15, 0, 15),
        x_title="Number of jets",
        discrete_x=True,
    )
    # pt of all jets in every event
    cfg.add_variable(
        name="jets_pt",
        expression="Jet.pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$p_{T}$ of all jets",
    )
    # pt of the first jet in every event
    cfg.add_variable(
        name="jet1_pt",  # variable name, to be given to the "--variables" argument for the plotting task
        expression="Jet.pt[:,0]",  # content of the variable
        null_value=EMPTY_FLOAT,  # value to be given if content not available for event
        binning=(50, 0.0, 500.0),  # (bins, lower edge, upper edge)
        unit="GeV",  # unit of the variable, if any
        x_title=r"Jet 1 $p_{T}$",  # x title of histogram when plotted
    )
    # eta of the first jet in every event
    cfg.add_variable(
        name="jet1_eta",
        expression="Jet.eta[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(30, -3.0, 3.0),
        x_title=r"Jet 1 $\eta$",
    )
    cfg.add_variable(
        name="ht",
        expression=lambda events: ak.sum(events.Jet.pt, axis=1),
        binning=(40, 0.0, 800.0),
        unit="GeV",
        x_title="HT",
    )
    # weights
    cfg.add_variable(
        name="mc_weight",
        expression="mc_weight",
        binning=(200, -10, 10),
        x_title="MC weight",
    )
    # cutflow variables
    cfg.add_variable(
        name="cf_jet1_pt",
        expression="cutflow.jet1_pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"Jet 1 $p_{T}$",
    )
 # add new variables

    cfg.add_variable(
        name="jet_phi",
        expression="Jet.phi",
        null_value=EMPTY_FLOAT,
        binning=(30, -3.14, 3.14),
        x_title=r"Jets $\phi$",
    )
    cfg.add_variable(
        name="jet_delta_phi",
        null_value=EMPTY_FLOAT,
        binning=(30, -3.14, 3.14),
        x_title=r"Jets $\Delta\phi_{1,2}$",
    )

    cfg.add_variable(
        name="jet_delta_r_12",
        null_value=EMPTY_FLOAT,
        binning=(30, 0, 6),
        x_title=r"Jets $\Delta R_{1,2}$",
    )

    cfg.add_variable(
        name="jet_delta_r_13",
        null_value=EMPTY_FLOAT,
        binning=(30, 0, 6),
        x_title=r"Jets $\Delta R_{1,3}$",
    )


    cfg.add_variable(
        name="ele_eta",
        expression="Electron.eta",
        null_value=EMPTY_FLOAT,
        binning=(30, -3.0, 3.0),
        x_title=r"Electrons $\eta$",
    )

    cfg.add_variable(
        name="electron_pt",
        expression="Electron.pt",
        binning=(40, 0.0, 300.0),
        unit="GeV",
        x_title=r"$p_{T}$ of all electrons",
    )

    cfg.add_variable(
        name="muon_pt",
        expression="Muon.pt",
        binning=(40, 0.0, 200.0),
        unit="GeV",
        x_title=r"$p_{T}$ of all muons",
    )



    cfg.add_variable(
        name="mhhh_hadron",
        binning=(60, 0.0, 1500.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{HHH}^{gen,hadron}$",
    )

    cfg.add_variable(
        name="delta_r_bb1_hadron",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $\Delta R^{gen,hadron}$",
    )

    cfg.add_variable(
        name="delta_r_bb2_hadron",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $\Delta R^{gen,hadron}$",
    )

    cfg.add_variable(
        name="delta_r_tautau_hadron",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$\tau\tau$ $\Delta R^{gen,hadron}$",
    )

    cfg.add_variable(
        name="cos_bb1_hadron",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $cos(\delta)^{gen,hadron}$",
    )

    cfg.add_variable(
        name="cos_bb2_hadron",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $cos(\delta)^{gen,hadron}$",
    )

    cfg.add_variable(
        name="cos_tautau_hadron",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $cos(\delta)^{gen,hadron}$",
    )

    cfg.add_variable(
        name="delta_r_h12_hadron",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,2}^{gen,hadron}$",
    )

    cfg.add_variable(
        name="delta_r_h13_hadron",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,3}^{gen,hadron}$",
    )

    cfg.add_variable(
        name="delta_r_h23_hadron",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{2,3}^{gen,hadron}$",
    )

    cfg.add_variable(
        name="cos_h12_hadron",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,2}^{gen,hadron}$",
    )

    cfg.add_variable(
        name="cos_h13_hadron",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,3}^{gen,hadron}$",
    )

    cfg.add_variable(
        name="cos_h23_hadron",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{2,3}^{gen,hadron}$",
    )


    # detector level

    cfg.add_variable(
        name="mhhh",
        binning=(60, 150.0, 1300.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{4b2\tau}$",
    )

    cfg.add_variable(
        name="cos_taulep",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $cos(\delta)$",
    )

    cfg.add_variable(
        name="delta_r_taulep",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$\tau\tau$ $\Delta R$",
    )

    cfg.add_variable(
        name="h3_mass",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H3}$",
    )

    cfg.add_variable(
        name="n_fatjet",
        expression="n_fatjet",
        binning=(5, 0, 5),
        x_title="Number of fat jets",
        discrete_x=True,
    )

    cfg.add_variable(
        name="h1_unsort_mass",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H1}^{unsorted}$",
    )

    cfg.add_variable(
        name="h2_unsort_mass",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H2}^{unsorted}$",
    )

    cfg.add_variable(
        name="delta_r_h12",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,2}$",
    )

    cfg.add_variable(
        name="delta_r_h13",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,3}$",
    )

    cfg.add_variable(
        name="delta_r_h23",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{2,3}$",
    )
    
    cfg.add_variable(
        name="delta_r_bb1",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $\Delta R$",
    )

    cfg.add_variable(
        name="delta_r_bb2",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $\Delta R$",
    )

    cfg.add_variable(
        name="cos_h12",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,2}$",
    )

    cfg.add_variable(
        name="cos_h13",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,3}$",
    )

    cfg.add_variable(
        name="cos_h23",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{2,3}$",
    )
    
    cfg.add_variable(
        name="cos_bb1",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $cos(\delta)$",
    )

    cfg.add_variable(
        name="cos_bb2",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $cos(\delta)$",
    )

    cfg.add_variable(
        name="h1_mass",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H1}$",
    )

    cfg.add_variable(
        name="h2_mass",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H2}$",
    )

    cfg.add_variable(
        name="m_3btaulep",
        binning=(60, 150.0, 1300.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{3b2\tau}, (b_3,hhbtag)$",
    )

    cfg.add_variable(
        name="m_3btaulep_pt",
        binning=(60, 150.0, 1300.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{3b2\tau}, (b_3,pt)$",
    )


    # gen-level variables
    cfg.add_variable(
        name="mtautau_gen",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{\tau\tau}^{gen}$",
    )

    cfg.add_variable(
        name="mbb_gen",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{bb}^{gen}$",
    )

    cfg.add_variable(
        name="mhhh_gen",
        binning=(60, 150.0, 1300.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{HHH}^{gen}$",
    )

    cfg.add_variable(
        name="mlnu_gen",
        binning=(40, 0.0, 200.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{l\nu}^{gen}$",
    )

    cfg.add_variable(
        name="hpt_gen",
        binning=(60, 0.0, 800.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$p_{TH}^{gen}$",
    )

    cfg.add_variable(
        name="h1bpt_gen",
        binning=(60, 0.0, 800.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$p_T_{H_1\rightarrow bb}^{gen}$",
    )

    cfg.add_variable(
        name="h2bpt_gen",
        binning=(60, 0.0, 800.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$p_T_{H_2\rightarrow bb}^{gen}$",
    )

    cfg.add_variable(
        name="htaupt_gen",
        binning=(60, 0.0, 800.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$p_T_{H\rightarrow\tau\tau}^{gen}$",
    )

    cfg.add_variable(
        name="delta_r_h12_gen",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,2}^{gen}$",
    )

    cfg.add_variable(
        name="delta_r_h13_gen",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{1,3}^{gen}$",
    )

    cfg.add_variable(
        name="delta_r_h23_gen",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $\Delta R_{2,3}^{gen}$",
    )
    
    cfg.add_variable(
        name="delta_r_bb1_gen",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $\Delta R^{gen}$",
    )

    cfg.add_variable(
        name="delta_r_bb2_gen",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $\Delta R^{gen}$",
    )

    cfg.add_variable(
        name="delta_r_tautau_gen",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$\tau\tau$ $\Delta R^{gen}$",
    )

    cfg.add_variable(
        name="cos_h12_gen",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,2}^{gen}$",
    )

    cfg.add_variable(
        name="cos_h13_gen",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{1,3}^{gen}$",
    )

    cfg.add_variable(
        name="cos_h23_gen",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $cos(\delta)_{2,3}^{gen}$",
    )
    
    cfg.add_variable(
        name="cos_bb1_gen",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $cos(\delta)^{gen}$",
    )

    cfg.add_variable(
        name="cos_bb2_gen",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $cos(\delta)^{gen}$",
    )

    cfg.add_variable(
        name="cos_tautau_gen",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$\tau\tau$ $cos(\delta)^{gen}$",
    )


    ### detector level but with experimental chi**2 minimization to H mass for jet pairing

    cfg.add_variable(
        name="min_chi",
        null_value=EMPTY_FLOAT,
        binning=(50, 0, 0.2),
        x_title=r"minimal $\chi^2$",
    )

    cfg.add_variable(
        name="mds_h1_mass_chi",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H1}$ $(\chi^2)$ (mds)",
    )

    cfg.add_variable(
        name="mds_h2_mass_chi",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H2}$ $(\chi^2)$ (mds)",
    )

    cfg.add_variable(
        name="delta_r_h12_chi",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $(\chi^2)$ $\Delta R_{1,2}$",
    )

    cfg.add_variable(
        name="delta_r_h13_chi",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $(\chi^2)$ $\Delta R_{1,3}$",
    )

    cfg.add_variable(
        name="delta_r_h23_chi",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"H $(\chi^2)$ $\Delta R_{2,3}$",
    )
    
    cfg.add_variable(
        name="delta_r_bb1_chi",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_1$ $(\chi^2)$ $\Delta R$",
    )

    cfg.add_variable(
        name="delta_r_bb2_chi",
        null_value=EMPTY_FLOAT,
        binning=(35, 0, 7),
        x_title=r"$bb_2$ $(\chi^2)$ $\Delta R$",
    )

    cfg.add_variable(
        name="cos_h12_chi",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $(\chi^2)$ $cos(\delta)_{1,2}$",
    )

    cfg.add_variable(
        name="cos_h13_chi",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $(\chi^2)$ $cos(\delta)_{1,3}$",
    )

    cfg.add_variable(
        name="cos_h23_chi",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"H $(\chi^2)$ $cos(\delta)_{2,3}$",
    )
    
    cfg.add_variable(
        name="cos_bb1_chi",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_1$ $(\chi^2)$ $cos(\delta)$",
    )

    cfg.add_variable(
        name="cos_bb2_chi",
        null_value=EMPTY_FLOAT,
        binning=(24, -1, +1),
        x_title=r"$bb_2$ $(\chi^2)$ $cos(\delta)$",
    )

    cfg.add_variable(
        name="h1_mass_chi",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H1}$ $(\chi^2)$",
    )

    cfg.add_variable(
        name="h2_mass_chi",
        binning=(40, 0.0, 400.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{H2}$ $(\chi^2)$",
    )

    cfg.add_variable(
        name="m_3btaulep_chi",
        binning=(60, 150.0, 1300.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{3bl\tau} (\chi^2), (b_3,hhbtag)$",
    )

    cfg.add_variable(
        name="m_3btaulep_pt_chi",
        binning=(60, 150.0, 1300.0),
        null_value=EMPTY_FLOAT,
        unit="GeV",
        x_title=r"$m_{3bl\tau} (\chi^2), (b_3,pt)$",
    )