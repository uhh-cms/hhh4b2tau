import functools

import order as od
from scinum import Number

from columnflow.util import DotDict
from columnflow.columnar_util import ColumnCollection
from columnflow.config_util import (
    get_root_processes_from_campaign, add_shift_aliases, get_shifts_from_sources,
    verify_config_processes,
)


def add_config(
    analysis: od.Analysis,
    campaign: od.Campaign,
    config_name: str or None = None,
    config_id: int or None = None,
    limit_dataset_files: int or None = None,
    **kwargs,
) -> od.Config:

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    # from IPython import embed; embed(header="After getting procs")
    # create a config by passing the campaign, so id and name will be identical
    cfg = od.Config(name=config_name, id=config_id, campaign=campaign)

    # gather campaign data
    year = campaign.x.year

    # add processes we are interested in
    from cmsdb.processes.hhh import __all__ as all_hhh_processes
    process_names = (
        [
            "data",
            "tt",
        ] + [
            x for x in all_hhh_processes
            if all(s in x for s in ["c3", "d4", "4b2tau"])
        ]
    )

    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))

        # configuration of colors, labels, etc. can happen here
        if proc.is_mc:         
            coupling_with_colors = (
                # (c3, d4, color) with colors recommended by cms
                  (0, 0, "#000000"),
                  (0, 99, "#3f90da"),
                  (0, 'm1', "#ffa90e"),
                  (19, 19, "#bd1f01"),
                  (1, 0, "#94a4a2"),
                  (1, 2, "#832db6"),
                  (2, 'm1', "#a96b59"),
                  (4, 9, "#e76300"),
                  ('m1', 0, "#b9ac70"),
                  ('m1', 'm1', "#717581"),
                  ('m1p5', 'm0p5', "#92dadd"),
            )

            for c3,d4,color in coupling_with_colors:
                if proc.name == f"hhh_4b2tau_c3{c3}_d4{d4}":
                    proc.color1 = color
            if proc.name == "tt":
                proc.color1 = "#e41a1c"   # else "#377eb8"

        # from hhh4b2tau.config.styles import stylize_processes
        # stylize_processes(cfg)

    # add datasets we need to study
    dataset_names = ([
        # data
        "data_mu_d",
        # backgrounds
        "tt_sl_powheg",        
    ] + [
            f"{x}_amcatnlo" for x in all_hhh_processes
            if all(s in x for s in ["c3", "d4", "4b2tau"])
        ]
    )
    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))

        # for testing purposes, limit the number of files
        if limit_dataset_files:
            for info in dataset.info.values():
                info.n_files = min(info.n_files, limit_dataset_files)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "example"
    cfg.x.default_selector = "example"
    cfg.x.default_producer = "example"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories = ("incl",)
    cfg.x.default_variables = ("n_jet", "jet1_pt")

    # set default weight_producer
    cfg.x.default_weight_producer = "all_weights"

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {}

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {
        "hhh_couplings": [
            f"{x}_amcatnlo" for x in all_hhh_processes
            if all(s in x for s in ["c3", "d4", "4b2tau"])
        ],
        "hhh_compare_1": [
            f"hhh_4b2tau_c3{x}_d4{y}_amcatnlo" for x,y in ((0, 0), (1, 0), ("m1", 0), (0, 99), (0, "m1"), (2, "m1"))
        ],

        "hhh_compare_2": [
            f"hhh_4b2tau_c3{x}_d4{y}_amcatnlo" for x,y in ((0, 0), (19, 19), (4, 9), ("m1p5", "m0p5"), ("m1", "m1"), (1, 2))
        ],
    }

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # general_settings groups for conveniently looping over different values for the general-settings parameter
    # (used during plotting)
    cfg.x.general_settings_groups = {
        "compare_shapes": {"skip_ratio": True, "shape_norm": True, "cms_label": "simpw"},
    }

    # process_settings groups for conveniently looping over different values for the process-settings parameter
    # (used during plotting)
    cfg.x.process_settings_groups = {
        "unstack_processes": {f"{x}": {"unstack": True} for x in all_hhh_processes
            if all(s in x for s in ["c3", "d4", "4b2tau"])},

        # "unstack_all": {"*": {"unstack": True}}, # does not work somehow????
    }

    # variable_settings groups for conveniently looping over different values for the variable-settings parameter
    # (used during plotting)
    cfg.x.variable_settings_groups = {}

    # custom_style_config groups for conveniently looping over certain style configs
    # (used during plotting)
    cfg.x.custom_style_config_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["muon", "jet"],
    }

    # calibrator groups for conveniently looping over certain calibrators
    # (used during calibration)
    cfg.x.calibrator_groups = {}

    # producer groups for conveniently looping over certain producers
    # (used during the ProduceColumns task)
    cfg.x.producer_groups = {}

    # ml_model groups for conveniently looping over certain ml_models
    # (used during the machine learning tasks)
    cfg.x.ml_model_groups = {}


    # custom method and sandbox for determining dataset lfns
    cfg.x.get_dataset_lfns = None
    cfg.x.get_dataset_lfns_sandbox = None

    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False

    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    if year == 2022:
        cfg.x.luminosity = Number(38010, {
            "total": 0.014j,
        })
    elif year == 2023:
        cfg.x.luminosity = Number(27208, {
            "lumi_13TeV_correlated": 0.0j,
        })
    else:  # 2024
        cfg.x.luminosity = Number(0, {
            "lumi_13TeV_correlated": 0.0j,
        })

    # names of muon correction sets and working points
    # (used in the muon producer)
    cfg.x.muon_sf_names = ("NUM_TightRelIso_DEN_TightIDandIPCut", f"{year}_UL")

    # register shifts
    cfg.add_shift(name="nominal", id=0)

    # tune shifts are covered by dedicated, varied datasets, so tag the shift as "disjoint_from_nominal"
    # (this is currently used to decide whether ML evaluations are done on the full shifted dataset)
    cfg.add_shift(name="tune_up", id=1, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="tune_down", id=2, type="shape", tags={"disjoint_from_nominal"})

    # fake jet energy correction shift, with aliases flaged as "selection_dependent", i.e. the aliases
    # affect columns that might change the output of the event selection
    cfg.add_shift(name="jec_up", id=20, type="shape")
    cfg.add_shift(name="jec_down", id=21, type="shape")
    add_shift_aliases(
        cfg,
        "jec",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            "MET.pt": "MET.pt_{name}",
            "MET.phi": "MET.phi_{name}",
        },
    )

    # event weights due to muon scale factors
    cfg.add_shift(name="mu_up", id=10, type="shape")
    cfg.add_shift(name="mu_down", id=11, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})

    # external files
    # json_mirror = "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-9ea86c4c"
    # cfg.x.external_files = DotDict.wrap({
    #     # lumi files
    #     "lumi": {
    #         "golden": ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt", "v1"),  # noqa
    #         "normtag": ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
    #     },

    #     # muon scale factors
    #     "muon_sf": (f"{json_mirror}/POG/MUO/{year}_UL/muon_Z.json.gz", "v1"),
    # })

# external files
    json_mirror = "/afs/cern.ch/user/a/anhaddad/public/jsonpog-integration"
    json_mirror_alt = "/afs/cern.ch/user/a/anhaddad/public/jsonpog_alt"
    # TODO later: add factors for other POGs when available
    year2="22"
    year_postfix = ""
    cfg.x.external_files = DotDict.wrap({
        # pileup weight corrections
        "pu_sf": (f"{json_mirror}/POG/LUM/{year}_Summer{year2}{year_postfix}/puWeights.json.gz", "v1"),

        # jet energy correction
        "jet_jerc": (f"{json_mirror}/POG/JME/{year}_Summer{year2}{year_postfix}/jet_jerc.json.gz", "v1"),

        # tau energy correction and scale factors
        # "tau_sf": (f"{json_mirror}/POG/TAU/{year_folder}/tau.json.gz", "v1"),

        # electron scale factors
        # "electron_sf": (f"{json_mirror}/POG/EGM/{year_folder}/electron.json.gz", "v1"),

        # muon scale factors
        # "muon_sf": (f"{json_mirror}/POG/MUO/{year_folder}/muon_Z.json.gz", "v1"),

        # btag scale factor
        "btag_sf_corr": (f"{json_mirror}/POG/BTV/{year}_Summer{year2}{year_postfix}/btagging.json.gz", "v1"),

        # met phi corrector
        # "met_phi_corr": (f"{json_mirror}/POG/JME/2018_UL/met.json.gz", "v1"),

        # hh-btag repository (lightweight) with TF saved model directories
        "hh_btag_repo": ("https://github.com/hh-italian-group/HHbtag/archive/df5220db5d4a32d05dc81d652083aece8c99ccab.tar.gz", "v2"),  # noqa
    })

    if year == 2022:
        cfg.x.external_files.update(DotDict.wrap({
            # Add Muon POG scale factors
            "muon_sf": (f"{json_mirror}/POG/MUO/{year}{year_postfix}_27Jun2023/muon_Z.json.gz", "v1"),

            # electron scale factors
            "electron_sf": (f"{json_mirror}/POG/EGM/{year}_Summer{year2}{year_postfix}/electron.json.gz", "v1"),

            # tau energy correction and scale factors
            #"tau_sf": (f"{json_mirror_alt}/POG/TAU/{year}_{postfixEE}/tau_DeepTau2018v2p5_2022_{postfixEE}.json.gz", "v1"),  # noqa

            # tau trigger
            #"tau_trigger_sf": (f"{json_mirror_alt}/POG/TAU/output/tau_trigger_DeepTau2018v2p5_{year}{postfixEE}.json", "v1"),  # noqa
        }))

    # external files with more complex year dependence # TODO: check this
    if year == 2022:
        cfg.x.external_files.update(DotDict.wrap({
            # lumi files
            "lumi": {
                "golden": ("/afs/cern.ch/user/a/anhaddad/public/Collisions22/Cert_Collisions2022_355100_362760_Golden.json", "v1"),  # noqa
                "normtag": ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json", "v1"),
            },

            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData?rev=45#Pileup_JSON_Files_For_Run_II
            "pu": {
                "json": ("/afs/cern.ch/user/a/anhaddad/public/Collisions22/pileup_JSON.txt", "v1"),  # noqa
                # Problem No file for 2022 --> using 2023 no matching shapes with root shape
                "mc_profile": ("https://raw.githubusercontent.com/cms-sw/cmssw/203834e3ae301f2564423dd1cc84bebf660519b9/SimGeneral/MixingModule/python/Run3_2022_LHC_Simulation_10h_2h_cfi.py", "v1"),  # noqa
                "data_profile": {
                    "nominal": (f"/afs/cern.ch/user/a/anhaddad/public/Collisions22/pileupHistogram-Cert_Collisions2022_355100_362760_GoldenJson-13p6TeV-69200ub-100bins.root", "v1"),  # noqa
                    "minbias_xs_up": (f"/afs/cern.ch/user/a/anhaddad/public/Collisions22/pileupHistogram-Cert_Collisions2022_355100_362760_GoldenJson-13p6TeV-72400ub-100bins.root", "v1"),  # noqa
                    "minbias_xs_down": (f"/afs/cern.ch/user/a/anhaddad/public/Collisions22/pileupHistogram-Cert_Collisions2022_355100_362760_GoldenJson-13p6TeV-66000ub-100bins.root", "v1"),  # noqa
                },
            },
        }))
    else:  # year 2023
        cfg.x.external_files.update(DotDict.wrap({
            # lumi files
            "lumi": {
                "golden": ("/afs/cern.ch/user/a/anhaddad/public/Collisions23/Cert_Collisions2023_366442_370790_Golden.json", "v1"),  # noqa
                "normtag": ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
            },

            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData?rev=45#Pileup_JSON_Files_For_Run_II
            "pu": {
                "json": ("/afs/cern.ch/user/a/anhaddad/public/Collisions23/pileup_JSON.txt", "v1"),  # noqa
                "mc_profile": ("https://raw.githubusercontent.com/cms-sw/cmssw/203834e3ae301f2564423dd1cc84bebf660519b9/SimGeneral/MixingModule/python/mix_2023_25ns_EraCD_PoissonOOTPU_cfi.py", "v1"),  # noqa
                "data_profile": {
                    "nominal": (f"/afs/cern.ch/user/a/anhaddad/public/Collisions23/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-69200ub-100bins.root", "v1"),  # noqa
                    "minbias_xs_up": (f"/afs/cern.ch/user/a/anhaddad/public/Collisions23/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-72400ub-100bins.root", "v1"),  # noqa
                    "minbias_xs_down": (f"/afs/cern.ch/user/a/anhaddad/public/Collisions23/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-66000ub-100bins.root", "v1"),  # noqa
                },
            },
        }))






    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0

    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # general event info, mandatory for reading files with coffea
            ColumnCollection.MANDATORY_COFFEA,  # additional columns can be added as strings, similar to object info
            # object info
            "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.btagDeepFlavB", "Jet.hadronFlavour",           
            "MET.pt", "MET.phi", "MET.significance", "MET.covXX", "MET.covXY", "MET.covYY",
            "Muon.*",
            "Electron.*",
            "Tau.pt", "Tau.eta", "Tau.phi", "Tau.mass", "Tau.idDeepTau2017v2p1VSe", "Tau.charge",
            "Tau.idDeepTau2017v2p1VSmu", "Tau.idDeepTau2017v2p1VSjet", "Tau.genPartFlav",
            "Tau.decayMode",
            "PV.npvs",
            # all columns added during selection using a ColumnCollection flag
            ColumnCollection.ALL_FROM_SELECTOR,
        },
        "cf.MergeSelectionMasks": {
            "cutflow.*",
        },
        "cf.UniteColumns": {
            "*",
        },
    })

    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        # "muon_weight": get_shifts("mu"),
    })

    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    cfg.x.versions = {
        # "cf.CalibrateEvents": "prod1",
        # "cf.SelectEvents": (lambda cls, inst, params: "prod1" if params.get("selector") == "default" else "dev1"),
        # ...
    }

    # channels
    # (just one for now)
    cfg.add_channel(name="mutau", id=1)

    # add categories using the "add_category" tool which adds auto-generated ids
    # the "selection" entries refer to names of categorizers, e.g. in categorization/example.py
    # note: it is recommended to always add an inclusive category with id=1 or name="incl" which is used
    #       in various places, e.g. for the inclusive cutflow plots and the "empty" selector
    
    from hhh4b2tau.config.categories import add_categories
    add_categories(cfg)

    from hhh4b2tau.config.variables import add_variables
    add_variables(cfg)

    return cfg