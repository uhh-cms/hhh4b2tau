from __future__ import annotations
import functools

import order as od
import law
import re
from scinum import Number

from columnflow.util import DotDict, dev_sandbox
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

    # gather campaign data
    run = campaign.x.run
    year = campaign.x.year
    year2 = year % 100

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    # from IPython import embed; embed(header="After getting procs")
    # create a config by passing the campaign, so id and name will be identical
    cfg = od.Config(name=config_name, id=config_id, campaign=campaign)

    # helper to enable processes / datasets only for a specific era
    def if_era(
        *,
        run: int | None = None,
        year: int | None = None,
        postfix: str | None = None,
        tag: str | None = None,
        values: list[str] | None = None,
    ) -> list[str]:
        match = (
            (run is None or campaign.x.run == run) and
            (year is None or campaign.x.year == year) and
            (postfix is None or campaign.x.postfix == postfix) and
            (tag is None or campaign.has_tag(tag))
        )
        return (values or []) if match else []

    # add processes we are interested in
    from cmsdb.processes.hhh import __all__ as all_hhh_processes
    process_names = (
        "data",
        "tt",
        "dy",
        # hhh signal processes
        *if_era(run=3, year=2022, tag="preEE", values=[x for x in all_hhh_processes
            if all(s in x for s in ["c3", "d4", "4b2tau"])
        ]),
        *if_era(run=3, year=2022, tag="preEE", values=[
            # hh background
            "hh_ggf_hbb_htt_kl1_kt1",
            "hh_ggf_hbb_htt_kl0_kt1",
            "hh_ggf_hbb_htt_kl2p45_kt1",
            "hh_ggf_hbb_htt_kl5_kt1",
            "hh_ggf_hbb_htt_kl0_kt1_c21",
            "hh_ggf_hbb_htt_kl1_kt1_c23",
            "hh_vbf_hbb_htt_kv1_k2v1_kl1",
            "hh_vbf_hbb_htt_kv1_k2v0_kl1",
            "hh_vbf_hbb_htt_kv1_k2v1_kl2",
            "hh_vbf_hbb_htt_kv1_k2v2_kl1",
            "hh_vbf_hbb_htt_kv1p74_k2v1p37_kl14p4",
            "hh_vbf_hbb_htt_kvm0p012_k2v0p03_kl10p2",
            "hh_vbf_hbb_htt_kvm0p758_k2v1p44_klm19p3",
            "hh_vbf_hbb_htt_kvm0p962_k2v0p959_klm1p43",
            "hh_vbf_hbb_htt_kvm1p21_k2v1p94_klm0p94",
            "hh_vbf_hbb_htt_kvm1p6_k2v2p72_klm1p36",
            "hh_vbf_hbb_htt_kvm1p83_k2v3p57_klm3p39",
            "hh_vbf_hbb_htt_kvm2p12_k2v3p87_klm5p96",
            # tth background
            "tth",
        ]),
    )

    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))

        

    from hhh4b2tau.config.styles import stylize_processes
    stylize_processes(cfg)

    # add datasets we need to study
    dataset_names = (
        # data
        "data_mu_d",
        # backgrounds
        "tt_sl_powheg",
        "tt_dl_powheg",
        "tt_fh_powheg",
        *if_era(run=3, year=2022, values=[
            "hh_ggf_hbb_htt_kl1_kt1_powheg",
            "hh_ggf_hbb_htt_kl0_kt1_powheg",
            "hh_ggf_hbb_htt_kl2p45_kt1_powheg",
            "hh_ggf_hbb_htt_kl5_kt1_powheg",
            # vbf
            "hh_vbf_hbb_htt_kv1_k2v1_kl1_madgraph",
            "hh_vbf_hbb_htt_kv1_k2v1_kl2_madgraph",
            "hh_vbf_hbb_htt_kv1_k2v0_kl1_madgraph",
            "hh_vbf_hbb_htt_kv1_k2v2_kl1_madgraph",
            "hh_vbf_hbb_htt_kv1p74_k2v1p37_kl14p4_madgraph",
            "hh_vbf_hbb_htt_kvm0p012_k2v0p03_kl10p2_madgraph",
            "hh_vbf_hbb_htt_kvm0p758_k2v1p44_klm19p3_madgraph",
            "hh_vbf_hbb_htt_kvm0p962_k2v0p959_klm1p43_madgraph",
            "hh_vbf_hbb_htt_kvm1p21_k2v1p94_klm0p94_madgraph",
            "hh_vbf_hbb_htt_kvm1p6_k2v2p72_klm1p36_madgraph",
            "hh_vbf_hbb_htt_kvm1p83_k2v3p57_klm3p39_madgraph",
            "hh_vbf_hbb_htt_kvm2p12_k2v3p87_klm5p96_madgraph",
            # dy
            "dy_m4to10_amcatnlo",
            "dy_m10to50_amcatnlo",
            "dy_m50toinf_amcatnlo",
            "dy_m50toinf_0j_amcatnlo",
            "dy_m50toinf_1j_amcatnlo",
            "dy_m50toinf_2j_amcatnlo",
            "dy_m50toinf_1j_pt40to100_amcatnlo",
            "dy_m50toinf_1j_pt100to200_amcatnlo",
            "dy_m50toinf_1j_pt200to400_amcatnlo",
            "dy_m50toinf_1j_pt400to600_amcatnlo",
            "dy_m50toinf_1j_pt600toinf_amcatnlo",
            "dy_m50toinf_2j_pt40to100_amcatnlo",
            "dy_m50toinf_2j_pt100to200_amcatnlo",
            "dy_m50toinf_2j_pt200to400_amcatnlo",
            "dy_m50toinf_2j_pt400to600_amcatnlo",
            "dy_m50toinf_2j_pt600toinf_amcatnlo",
            # ttH
            "tth_hbb_powheg",
            "tth_hnonbb_powheg",
        ]),
        *if_era(run=3, year=2022, values=[
            f"{x}_amcatnlo" for x in all_hhh_processes
            if all(s in x for s in ["c3", "d4", "4b2tau"])
        ]),
    )
    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))

        # add tags to datasets
        if dataset.name.startswith("tt"):
            dataset.add_tag(("has_top", "is_ttbar"))
        elif dataset.name.startswith("st"):
            dataset.add_tag(("has_top", "is_single_top"))
        if dataset.name.startswith("dy"):
            dataset.add_tag("is_dy")
        if re.match(r"^(ww|wz|zz)_.*pythia$", dataset.name):
            dataset.add_tag("no_lhe_weights")
        if dataset_name.startswith("hh_"):
            dataset.add_tag("signal")
            dataset.add_tag("nonresonant_signal")
        if dataset_name.startswith(("graviton_hh_", "radion_hh_")):
            dataset.add_tag("signal")
            dataset.add_tag("resonant_signal")

        # apply an optional limit on the number of files
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

    # add a hist hook to work with histograms of data samples
    # i.e. morph data for hypothetical coupling that has not been generated yet
    import hhh4b2tau.plotting.morphing as morphing

    cfg.x.hist_hooks = {
    "morphing": morphing.morphing_hook,
    }

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "backgrounds": (backgrounds := [
            # "h",
            "tt",
            "dy",
            # "qcd",
            # "st",
            # "v",
            # "multiboson",
            # "tt_multiboson",
            # "ewk",
            # "tth",
        ]),
         "split_backgrounds": (split_backgrounds := [
            # "h",
            "tt_sl",
            "tt_dl",
            "tt_fh",
            "dy",
            # "qcd",
            # "st",
            # "v",
            # "multiboson",
            # "tt_multiboson",
            # "ewk",
            # "tth",
        ]),
        "hhh_couplings": [
            f"{x}" for x in all_hhh_processes
            if all(s in x for s in ["c3", "d4", "4b2tau"])
        ],
        "hhh_compare_1": [
            f"hhh_4b2tau_c3{x}_d4{y}" for x,y in ((0, 0), (1, 0), ("m1", 0), (0, 99), (0, "m1"), (2, "m1"))
        ],

        "hhh_compare_2": [
            f"hhh_4b2tau_c3{x}_d4{y}" for x,y in ((0, 0), (19, 19), (4, 9), ("m1p5", "m0p5"), ("m1", "m1"), (1, 2))
        ],
        "sm_higgs": (sm_higgs := [
            "tth",
            "hhh_4b2tau_c30_d40",
            "hh_ggf_hbb_htt_kl1_kt1",
        ]),
        "sm": sorted(list(set(split_backgrounds + sm_higgs))),
    }

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
        "hhh_compare_to_morph": [
            "hhh_4b2tau_c3{c3}_d4{d4}_amcatnlo".format(
                      c3=str(c3).replace("-", "m").replace(".", "p"),
                      d4=str(d4).replace("-", "m").replace(".", "p"),
                      ) for c3,d4 in morphing.morphing_coupling_combinations
        ],
        "sm_higgs": (sm_higgs := [
            "tth_hbb_powheg",
            "tth_hnonbb_powheg",
            "hhh_4b2tau_c30_d40_amcatnlo",
            "hh_ggf_hbb_htt_kl1_kt1_powheg",
        ]),
    }

    # define inclusive datasets for the dy process identification with corresponding leaf processes
    if run == 3:
        cfg.x.dy_stitching = {
            "m50toinf": {
                "inclusive_dataset": cfg.datasets.n.dy_m50toinf_amcatnlo,
                "leaf_processes": [
                    # the following processes cover the full njet and pt phasespace
                    procs.n.dy_m50toinf_0j,
                    *(
                        procs.get(f"dy_m50toinf_{nj}j_pt{pt}")
                        for nj in [1, 2]
                        for pt in ["0to40", "40to100", "100to200", "200to400", "400to600", "600toinf"]
                    ),
                    procs.n.dy_m50toinf_ge3j,
                ],
            },
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
        "unstack_morph": {"hhh_4b2tau_c3{c3}_d4{d4}".format(
                      c3=str(c3).replace("-", "m").replace(".", "p"),
                      d4=str(d4).replace("-", "m").replace(".", "p"),
                      ) : {"unstack": True}
                      for c3,d4 in morphing.all_cc}
    }

    # variable_settings groups for conveniently looping over different values for the variable-settings parameter
    # (used during plotting)
    cfg.x.variable_settings_groups = {}

    # custom_style_config groups for conveniently looping over certain style configs
    # (used during plotting)
    cfg.x.custom_style_config_groups = {
        "small_legend": {
            "legend_cfg": {"ncols": 2, "fontsize": 16, "columnspacing": 0.6},
        },
    }
    cfg.x.default_custom_style_config = "small_legend"

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

    ################################################################################################
    # external files
    ################################################################################################

    cfg.x.external_files = DotDict()

    # helper
    def add_external(name, value):
        if isinstance(value, dict):
            value = DotDict.wrap(value)
        cfg.x.external_files[name] = value

    if run == 2:
        json_postfix = ""
        if year == 2016:
            json_postfix = f"{'pre' if campaign.has_tag('preVFP') else 'post'}VFP"
        json_pog_era = f"{year}{json_postfix}_UL"
        json_mirror = "/afs/cern.ch/user/m/mrieger/public/mirrors/jsonpog-integration-6ce37404"
    elif run == 3:
        json_pog_era = f"{year}_Summer{year2}{campaign.x.postfix}"
        json_mirror = "/afs/cern.ch/user/m/mrieger/public/mirrors/jsonpog-integration-6ce37404"
    else:
        assert False

    # common files
    # lumi files
    add_external("lumi", {
        "golden": {
            2016: ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt", "v1"),  # noqa
            2017: ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt", "v1"),  # noqa
            2018: ("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt", "v1"),  # noqa,
            # TODO: document source
            2022: ("/afs/cern.ch/user/a/anhaddad/public/Collisions22/Cert_Collisions2022_355100_362760_Golden.json", "v1"),  # noqa
            2023: ("/afs/cern.ch/user/a/anhaddad/public/Collisions23/Cert_Collisions2023_366442_370790_Golden.json", "v1"),  # noqa
        }[year],
        "normtag": {
            2016: ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
            2017: ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
            2018: ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
            # TODO: check
            2022: ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json", "v1"),
            2023: ("/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json", "v1"),
        }[year],
    })
    # pileup weight corrections
    add_external("pu_sf", (f"{json_mirror}/POG/LUM/{json_pog_era}/puWeights.json.gz", "v1"))
    # jet energy correction
    add_external("jet_jerc", (f"{json_mirror}/POG/JME/{json_pog_era}/jet_jerc.json.gz", "v1"))
    # jet veto map
    add_external("jet_veto_map", (f"{json_mirror}/POG/JME/{json_pog_era}/jetvetomaps.json.gz", "v1"))
    # btag scale factor
    add_external("btag_sf_corr", (f"{json_mirror}/POG/BTV/{json_pog_era}/btagging.json.gz", "v1"))
    # hh-btag repository (lightweight) with TF saved model directories
    add_external("hh_btag_repo", ("https://github.com/hh-italian-group/HHbtag/archive/df5220db5d4a32d05dc81d652083aece8c99ccab.tar.gz", "v2"))  # noqa
    # Tobias' tautauNN (https://github.com/uhh-cms/tautauNN)
    add_external("res_pdnn", ("/afs/cern.ch/work/m/mrieger/public/hbt/models/res_prod3/model_fold0.tgz", "v1"))

    # run specific files
    if run == 2:
        # tau energy correction and scale factors
        add_external("tau_sf", (f"{json_mirror}/POG/TAU/{json_pog_era}/tau.json.gz", "v1"))
        # tau trigger scale factors
        add_external("tau_trigger_sf", (f"{json_mirror}/POG/TAU/{json_pog_era}/tau.json.gz", "v1"))
        # electron scale factors
        add_external("electron_sf", (f"{json_mirror}/POG/EGM/{json_pog_era}/electron.json.gz", "v1"))
        # muon scale factors
        add_external("muon_sf", (f"{json_mirror}/POG/MUO/{json_pog_era}/muon_Z.json.gz", "v1"))
        # met phi correction
        add_external("met_phi_corr", (f"{json_mirror}/POG/JME/{json_pog_era}/met.json.gz", "v1"))
    elif run == 3:
        if year == 2022 and campaign.has_tag("preEE"):
            # muon scale factors
            add_external("muon_sf", (f"{json_mirror}/POG/MUO/{json_pog_era}/muon_Z.json.gz", "v1"))
            # electron scale factors
            add_external("electron_sf", (f"{json_mirror}/POG/EGM/{json_pog_era}/electron.json.gz", "v1"))
            # tau energy correction and scale factors
            # TODO: remove tag pog mirror once integrated centrally
            json_mirror_tau_pog = "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-taupog"
            tau_pog_era = f"{year}_{'pre' if campaign.has_tag('preEE') else 'post'}EE"
            add_external("tau_sf", (f"{json_mirror_tau_pog}/POG/TAU/{tau_pog_era}/tau_DeepTau2018v2p5_{tau_pog_era}.json.gz", "v1"))  # noqa
    else:
        assert False

    ################################################################################################
    # tau settings
    ################################################################################################

    # tau tagger name
    # (needed by TECConfig below as well as tau selection)
    if run == 2:
        # TODO: still correct? what about 2p5?
        cfg.x.tau_tagger = "DeepTau2017v2p1"
    elif run == 3:
        # https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
        cfg.x.tau_tagger = "DeepTau2018v2p5"
    else:
        assert False

    # tec config
    from columnflow.calibration.cms.tau import TECConfig
    corrector_kwargs = {"wp": "Tight", "wp_VSe": "Tight"} if run == 3 else {}
    cfg.x.tec = TECConfig(tagger=cfg.x.tau_tagger, corrector_kwargs=corrector_kwargs)

    # tau ID working points
    if campaign.x.version < 10:
        cfg.x.tau_id_working_points = DotDict.wrap({
            "tau_vs_e": {"vvvloose": 1, "vvloose": 2, "vloose": 4, "loose": 8, "medium": 16, "tight": 32, "vtight": 64, "vvtight": 128},  # noqa
            "tau_vs_jet": {"vvvloose": 1, "vvloose": 2, "vloose": 4, "loose": 8, "medium": 16, "tight": 32, "vtight": 64, "vvtight": 128},  # noqa
            "tau_vs_mu": {"vloose": 1, "loose": 2, "medium": 4, "tight": 8},
        })
    else:
        cfg.x.tau_id_working_points = DotDict.wrap({
            "tau_vs_e": {"vvvloose": 1, "vvloose": 2, "vloose": 3, "loose": 4, "medium": 5, "tight": 6, "vtight": 7, "vvtight": 8},  # noqa
            "tau_vs_jet": {"vvvloose": 1, "vvloose": 2, "vloose": 3, "loose": 4, "medium": 5, "tight": 6, "vtight": 7, "vvtight": 8},  # noqa
            "tau_vs_mu": {"vloose": 1, "loose": 2, "medium": 3, "tight": 4},
        })

    # tau trigger working points
    cfg.x.tau_trigger_working_points = DotDict.wrap({
        "id_vs_jet_v0": "VVLoose",
        "id_vs_jet_gv0": ("Loose", "VVLoose"),
        "id_vs_mu_single": "Tight",
        "id_vs_mu_cross": "VLoose",
        "id_vs_e_single": "VVLoose",
        "id_vs_e_cross": "VVLoose",
        "trigger_corr": "VVLoose",
    })

    ################################################################################################
    # electron settings
    ################################################################################################

    # names of electron correction sets and working points
    # (used in the electron_sf producer)
    if run == 2:
        e_postfix = ""
        if year == 2016:
            e_postfix = "preVFP" if campaign.has_tag("preVFP") else "postVFP"
        cfg.x.electron_sf_names = (
            "UL-Electron-ID-SF",
            f"{year}{e_postfix}",
            "wp80iso",
        )
    elif run == 3 and year == 2022:
        cfg.x.electron_sf_names = (
            "Electron-ID-SF",
            "2022Re-recoBCD" if campaign.has_tag("preEE") else "2022Re-recoE+PromptFG",
            "wp80iso",
        )
    else:
        assert False

    ################################################################################################
    # muon settings
    ################################################################################################

    # names of muon correction sets and working points
    # (used in the muon producer)
    if run == 2:
        mu_postfix = ""
        if year == 2016:
            mu_postfix = "preVFP" if campaign.has_tag("preVFP") else "postVFP"
        cfg.x.muon_sf_names = (
            "NUM_TightRelIso_DEN_TightIDandIPCut",
            f"{year}{mu_postfix}_UL",
        )
    elif run == 3 and year == 2022:
        cfg.x.muon_sf_names = (
            "NUM_TightPFIso_DEN_TightID",
            "2022_preEE" if campaign.has_tag("preEE") else "2022_postEE",
        )
    else:
        assert False

    ################################################################################################
    # met settings
    ################################################################################################

    # name of the MET phi correction set
    # (used in the met_phi calibrator)
    cfg.x.met_phi_correction_set = r"{variable}_metphicorr_pfmet_{data_source}"

    ################################################################################################
    # b tagging
    ################################################################################################

    # b-tag working points
    btag_key = f"{year}{campaign.x.postfix}"
    if run == 2:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP?rev=6
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP?rev=8
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17?rev=15
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18?rev=18
        cfg.x.btag_working_points = DotDict.wrap({
            "deepjet": {
                "loose": {"2016APV": 0.0508, "2016": 0.0480, "2017": 0.0532, "2018": 0.0490}[btag_key],
                "medium": {"2016APV": 0.2598, "2016": 0.2489, "2017": 0.3040, "2018": 0.2783}[btag_key],
                "tight": {"2016APV": 0.6502, "2016": 0.6377, "2017": 0.7476, "2018": 0.7100}[btag_key],
            },
            "deepcsv": {
                "loose": {"2016APV": 0.2027, "2016": 0.1918, "2017": 0.1355, "2018": 0.1208}[btag_key],
                "medium": {"2016APV": 0.6001, "2016": 0.5847, "2017": 0.4506, "2018": 0.4168}[btag_key],
                "tight": {"2016APV": 0.8819, "2016": 0.8767, "2017": 0.7738, "2018": 0.7665}[btag_key],
            },
        })
    elif run == 3:
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
        # TODO later: complete WP when data becomes available
        cfg.x.btag_working_points = DotDict.wrap({
            "deepjet": {
                "loose": {"2022": 0.0583, "2022EE": 0.0614, "2023": 0.0, "2024": 0.0}[btag_key],
                "medium": {"2022": 0.3086, "2022EE": 0.3196, "2023": 0.0, "2024": 0.0}[btag_key],
                "tight": {"2022": 0.7183, "2022EE": 0.73, "2023": 0.0, "2024": 0.0}[btag_key],
                "xtight": {"2022": 0.8111, "2022EE": 0.8184, "2023": 0.0, "2024": 0.0}[btag_key],
                "xxtight": {"2022": 0.9512, "2022EE": 0.9542, "2023": 0.0, "2024": 0.0}[btag_key],
            },
            "robustParticleTransformer": {
                "loose": {"2022": 0.0849, "2022EE": 0.0897, "2023": 0.0, "2024": 0.0}[btag_key],
                "medium": {"2022": 0.4319, "2022EE": 0.451, "2023": 0.0, "2024": 0.0}[btag_key],
                "tight": {"2022": 0.8482, "2022EE": 0.8604, "2023": 0.0, "2024": 0.0}[btag_key],
                "xtight": {"2022": 0.9151, "2022EE": 0.9234, "2023": 0.0, "2024": 0.0}[btag_key],
                "xxtight": {"2022": 0.9874, "2022EE": 0.9893, "2023": 0.0, "2024": 0.0}[btag_key],
            },
            "particleNet": {
                "loose": {"2022": 0.047, "2022EE": 0.0499, "2023": 0.0, "2024": 0.0}[btag_key],
                "medium": {"2022": 0.245, "2022EE": 0.2605, "2023": 0.0, "2024": 0.0}[btag_key],
                "tight": {"2022": 0.6734, "2022EE": 0.6915, "2023": 0.0, "2024": 0.0}[btag_key],
                "xtight": {"2022": 0.7862, "2022EE": 0.8033, "2023": 0.0, "2024": 0.0}[btag_key],
                "xxtight": {"2022": 0.961, "2022EE": 0.9664, "2023": 0.0, "2024": 0.0}[btag_key],
            },
        })
    else:
        assert False

    # JEC uncertainty sources propagated to btag scale factors
    # (names derived from contents in BTV correctionlib file)
    cfg.x.btag_sf_jec_sources = [
        "",  # same as "Total"
        "Absolute",
        "AbsoluteMPFBias",
        "AbsoluteScale",
        "AbsoluteStat",
        f"Absolute_{year}",
        "BBEC1",
        f"BBEC1_{year}",
        "EC2",
        f"EC2_{year}",
        "FlavorQCD",
        "Fragmentation",
        "HF",
        f"HF_{year}",
        "PileUpDataMC",
        "PileUpPtBB",
        "PileUpPtEC1",
        "PileUpPtEC2",
        "PileUpPtHF",
        "PileUpPtRef",
        "RelativeBal",
        "RelativeFSR",
        "RelativeJEREC1",
        "RelativeJEREC2",
        "RelativeJERHF",
        "RelativePtBB",
        "RelativePtEC1",
        "RelativePtEC2",
        "RelativePtHF",
        "RelativeSample",
        f"RelativeSample_{year}",
        "RelativeStatEC",
        "RelativeStatFSR",
        "RelativeStatHF",
        "SinglePionECAL",
        "SinglePionHCAL",
        "TimePtEta",
    ]

    from columnflow.production.cms.btag import BTagSFConfig
    cfg.x.btag_sf = BTagSFConfig(
        correction_set="particleNet_shape",
        jec_sources=cfg.x.btag_sf_jec_sources,
        discriminator="btagPNetB",
    )

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
            "PFCandidate.*",
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
    ################################################################################################
    # weights
    ################################################################################################

    # configurations for all possible event weight columns as keys in an OrderedDict,
    # mapped to shift instances they depend on
    # (this info is used by weight producers)
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        # "pdf_weight": get_shifts("pdf"),
        # "murmuf_weight": get_shifts("murmuf"),
        # "normalized_pu_weight": get_shifts("minbias_xs"),
        # "normalized_njet_btag_weight": get_shifts(*(f"btag_{unc}" for unc in cfg.x.btag_unc_names)),
        # "electron_weight": get_shifts("e"),
        # "muon_weight": get_shifts("mu"),
        # "tau_weight": get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names)),
        # "tau_trigger_weight": get_shifts("etau_trigger", "mutau_trigger", "tautau_trigger"),
    })

    # define per-dataset event weights
    # for dataset in cfg.datasets:
    #     if dataset.has_tag("is_ttbar"):
    #         dataset.x.event_weights = {"top_pt_weight": get_shifts("top_pt")}

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

    from hhh4b2tau.config.met_filters import add_met_filters
    add_met_filters(cfg)

    ################################################################################################
    # LFN settings
    ################################################################################################

    # custom method and sandbox for determining dataset lfns
    cfg.x.get_dataset_lfns = None
    cfg.x.get_dataset_lfns_sandbox = None

    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    cfg.x.validate_dataset_lfns = limit_dataset_files is None

    # custom lfn retrieval method in case the underlying campaign is custom uhh
    if cfg.campaign.x("custom", {}).get("creator") == "uhh":
        def get_dataset_lfns(
            dataset_inst: od.Dataset,
            shift_inst: od.Shift,
            dataset_key: str,
        ) -> list[str]:
            # destructure dataset_key into parts and create the lfn base directory
            dataset_id, full_campaign, tier = dataset_key.split("/")[1:]
            main_campaign, sub_campaign = full_campaign.split("-", 1)
            lfn_base = law.wlcg.WLCGDirectoryTarget(
                f"/store/{dataset_inst.data_source}/{main_campaign}/{dataset_id}/{tier}/{sub_campaign}/0",
                fs=f"wlcg_fs_{cfg.campaign.x.custom['name']}",
            )

            # loop though files and interpret paths as lfns
            return [
                lfn_base.child(basename, type="f").path
                for basename in lfn_base.listdir(pattern="*.root")
            ]

        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns

        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")

        # define custom remote fs's to look at
        cfg.x.get_dataset_lfns_remote_fs = lambda dataset_inst: [
            f"local_fs_{cfg.campaign.x.custom['name']}",
            f"wlcg_fs_{cfg.campaign.x.custom['name']}",
        ]

    return cfg