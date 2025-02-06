# coding: utf-8

"""
Configuration of the HHH → bbbb𝜏𝜏 analysis.
"""

from __future__ import annotations

import os
import itertools
import functools

import yaml
import law
import order as od
from scinum import Number
import re


from columnflow.util import DotDict, dev_sandbox
from columnflow.config_util import (
    get_root_processes_from_campaign, add_shift_aliases, get_shifts_from_sources,
    verify_config_processes,
)
from columnflow.columnar_util import ColumnCollection, skip_column

thisdir = os.path.dirname(os.path.abspath(__file__))

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
        # "data_mu_d",
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

            # data
            *if_era(year=2022, tag="preEE", values=[
                f"data_{stream}_{period}" for stream in ["mu", "e", "tau"] for period in "cd"
            ]),
            *if_era(year=2022, tag="postEE", values=[
                f"data_{stream}_{period}" for stream in ["mu", "e", "tau"] for period in "efg"
            ]),
            *if_era(year=2023, tag="preBPix", values=[
                f"data_{stream}_c{v}" for stream in ["mu", "e", "tau"] for v in "1234"
            ]),
            *if_era(year=2023, tag="postBPix", values=[
                f"data_{stream}_d{v}" for stream in ["mu", "e", "tau"] for v in "12"
            ]),
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
        if dataset.name.startswith("data_e_"):
            dataset.add_tag({"etau", "emu_from_e", "ee"})
        if dataset.name.startswith("data_mu_"):
            dataset.add_tag({"mutau", "emu_from_mu", "mumu"})    
        if dataset.name.startswith("data_tau_"):
            dataset.add_tag({"tautau"})    
        if dataset.name.startswith("tt_"):
            dataset.add_tag({"has_top", "ttbar", "tt"})
        elif dataset.name.startswith("st"):
            dataset.add_tag({"has_top", "single_top", "st"})
        if dataset.name.startswith("dy"):
            dataset.add_tag("is_dy")
        if dataset.name.startswith("w_lnu_"):
            dataset.add_tag("w_lnu")
        # datasets that are known to have no lhe info at all
        if law.util.multi_match(dataset.name, [
            r"^(ww|wz|zz)_.*pythia$",
            r"^tt(w|z)_.*amcatnlo$",
        ]):
            dataset.add_tag("no_lhe_weights")
        # datasets that are allowed to contain some events with missing lhe infos
        # (known to happen for amcatnlo)
        if dataset.name.endswith("_amcatnlo"):
            dataset.add_tag("partial_lhe_weights")
        if dataset_name.startswith("hh_"):
            dataset.add_tag("signal")
            dataset.add_tag("nonresonant_signal")
        if dataset_name.startswith(("graviton_hh_", "radion_hh_")):
            dataset.add_tag("signal")
            dataset.add_tag("resonant_signal")
            if dataset_name.startswith(("graviton_hh_ggf_", "radion_hh_ggf")):
                dataset.add_tag("ggf")
            elif dataset_name.startswith(("graviton_hh_vbf_", "radion_hh_vbf")):
                dataset.add_tag("vbf")
        if re.match(r"^(ww|wz|zz)_.*pythia$", dataset.name):
            dataset.add_tag("no_lhe_weights")
        
        # apply an optional limit on the number of files
        if limit_dataset_files:
            for info in dataset.info.values():
                info.n_files = min(info.n_files, limit_dataset_files)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "default"
    cfg.x.default_selector = "default"
    cfg.x.default_producer = "default"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories = ("incl",)
    cfg.x.default_variables = ("n_jet", "jet1_pt")

    # set default weight_producer
    cfg.x.default_weight_producer = "default"

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
            # "dy",
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

        "hhh_compare_3": [
            f"hhh_4b2tau_c3{x}_d4{y}" for x,y in ((1, 0), (4, 9))
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
    cfg.x.variable_groups = {
        "all": ["delta_r_bb1", "delta_r_bb2", "delta_r_taulep", 
                "delta_r_h12", "delta_r_h13", "delta_r_h23",
                "cos_bb1", "cos_bb2", "cos_taulep",
                "cos_h12", "cos_h13", "cos_h23",
                "mhhh", "h1_mass", "h2_mass", "h3_mass",
                "n_b_jet", "n_fatjet",
                "m_3btaulep", "m_3btaulep_pt",
                # "h1_unsort_mass", "h2_unsort_mass",
                "delta_r_bb1_chi", "delta_r_bb2_chi",
                "delta_r_h12_chi", "delta_r_h13_chi", "delta_r_h23_chi",
                "cos_bb1_chi", "cos_bb2_chi",
                "cos_h12_chi", "cos_h13_chi", "cos_h23_chi",
                "h1_mass_chi", "h2_mass_chi",
                "m_3btaulep_chi", "m_3btaulep_pt_chi",
                "mds_h1_mass_chi", "mds_h2_mass_chi", "min_chi"
                ],
    }

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
        "default": ["json", "trigger", "met_filter", "jet_veto_map", "lepton", "jet2", "bjet"],
        "3b2tau": ["one_jet", "two_jet", "three_jet", "one_tau", "two_tau",],
        "var": ["delta_r_bb1", "delta_r_bb2", "delta_r_tautau", 
                "delta_r_h12", "delta_r_h13", "delta_r_h23",
                "cos_bb1", "cos_bb2", "cos_tautau",
                "cos_h12", "cos_h13", "cos_h23",
                "mhhh", "h3_mass", "m_3b2tau", "m_3b2tau_pt",],
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

    ################################################################################################
    # luminosity and normalization
    ################################################################################################

    # lumi values in 1/pb (= 1000/fb)
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=7
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3?rev=25
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    # difference pre-post VFP: https://cds.cern.ch/record/2854610/files/DP2023_006.pdf
    if year == 2016 and campaign.has_tag("preVFP"):
        cfg.x.luminosity = Number(19_500, {
            "lumi_13TeV_2016": 0.01j,
            "lumi_13TeV_correlated": 0.006j,
        })
    elif year == 2016 and campaign.has_tag("postVFP"):
        cfg.x.luminosity = Number(16_800, {
            "lumi_13TeV_2016": 0.01j,
            "lumi_13TeV_correlated": 0.006j,
        })
    elif year == 2017:
        cfg.x.luminosity = Number(41_480, {
            "lumi_13TeV_2017": 0.02j,
            "lumi_13TeV_1718": 0.006j,
            "lumi_13TeV_correlated": 0.009j,
        })
    elif year == 2018:
        cfg.x.luminosity = Number(59_830, {
            "lumi_13TeV_2017": 0.015j,
            "lumi_13TeV_1718": 0.002j,
            "lumi_13TeV_correlated": 0.02j,
        })
    elif year == 2022 and campaign.has_tag("preEE"):
        cfg.x.luminosity = Number(7_980.4, {
            "lumi_13p6TeV_correlated": 0.014j,
        })
    elif year == 2022 and campaign.has_tag("postEE"):
        cfg.x.luminosity = Number(26_671.7, {
            "lumi_13p6TeV_correlated": 0.014j,
        })
    elif year == 2023 and campaign.has_tag("preBPix"):
        cfg.x.luminosity = Number(17_794, {
            "lumi_13p6TeV_correlated": 0.013j,
        })
    elif year == 2023 and campaign.has_tag("postBPix"):
        cfg.x.luminosity = Number(9_451, {
            "lumi_13p6TeV_correlated": 0.013j,
        })
    else:
        assert False

    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=52#Recommended_cross_section
    cfg.x.minbias_xs = Number(69.2, 0.046j)

    ################################################################################################
    # met settings
    ################################################################################################

    if run == 2:
        cfg.x.met_name = "MET"
        cfg.x.raw_met_name = "RawMET"
    elif run == 3:
        cfg.x.met_name = "PuppiMET"
        cfg.x.raw_met_name = "RawPuppiMET"
    else:
        assert False

    # name of the MET phi correction set
    # (used in the met_phi calibrator)
    if run == 2:
        cfg.x.met_phi_correction_set = r"{variable}_metphicorr_pfmet_{data_source}"

    ################################################################################################
    # jet settings
    # TODO: keep a single table somewhere that configures all settings: btag correlation, year
    #       dependence, usage in calibrator, etc
    ################################################################################################

    # common jec/jer settings configuration
    if run == 2:
        # https://cms-jerc.web.cern.ch/Recommendations/#run-2
        # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=204
        # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=109
        jec_campaign = f"Summer19UL{year2}{campaign.x.postfix}"
        jec_version = {2016: "V7", 2017: "V5", 2018: "V5"}[year]
        jer_campaign = f"Summer{'20' if year == 2016 else '19'}UL{year2}{campaign.x.postfix}"
        jer_version = "JR" + {2016: "V3", 2017: "V2", 2018: "V2"}[year]
        jet_type = "AK4PFchs"
    elif run == 3:
        # https://cms-jerc.web.cern.ch/Recommendations/#2022
        jerc_postfix = {2022: "_22Sep2023", 2023: "Prompt23"}[year]
        jec_campaign = f"Summer{year2}{campaign.x.postfix}{jerc_postfix}"
        jec_version = {2022: "V2", 2023: "V1"}[year]
        jer_campaign = f"Summer{year2}{campaign.x.postfix}{jerc_postfix}"
        # special "Run" fragment in 2023 jer campaign
        if year == 2023:
            jer_campaign += f"_Run{'Cv1234' if campaign.has_tag('preBPix') else 'D'}"
        jer_version = "JR" + {2022: "V1", 2023: "V1"}[year]
        jet_type = "AK4PFPuppi"
    else:
        assert False

    cfg.x.jec = DotDict.wrap({
        "Jet": {
            "campaign": jec_campaign,
            "version": jec_version,
            "jet_type": jet_type,
            "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],
            "levels_for_type1_met": ["L1FastJet"],
            "uncertainty_sources": [
                # "AbsoluteStat",
                # "AbsoluteScale",
                # "AbsoluteSample",
                # "AbsoluteFlavMap",
                # "AbsoluteMPFBias",
                # "Fragmentation",
                # "SinglePionECAL",
                # "SinglePionHCAL",
                # "FlavorQCD",
                # "TimePtEta",
                # "RelativeJEREC1",
                # "RelativeJEREC2",
                # "RelativeJERHF",
                # "RelativePtBB",
                # "RelativePtEC1",
                # "RelativePtEC2",
                # "RelativePtHF",
                # "RelativeBal",
                # "RelativeSample",
                # "RelativeFSR",
                # "RelativeStatFSR",
                # "RelativeStatEC",
                # "RelativeStatHF",
                # "PileUpDataMC",
                # "PileUpPtRef",
                # "PileUpPtBB",
                # "PileUpPtEC1",
                # "PileUpPtEC2",
                # "PileUpPtHF",
                # "PileUpMuZero",
                # "PileUpEnvelope",
                # "SubTotalPileUp",
                # "SubTotalRelative",
                # "SubTotalPt",
                # "SubTotalScale",
                # "SubTotalAbsolute",
                # "SubTotalMC",
                "Total",
                # "TotalNoFlavor",
                # "TotalNoTime",
                # "TotalNoFlavorNoTime",
                # "FlavorZJet",
                # "FlavorPhotonJet",
                # "FlavorPureGluon",
                # "FlavorPureQuark",
                # "FlavorPureCharm",
                # "FlavorPureBottom",
                "CorrelationGroupMPFInSitu",
                "CorrelationGroupIntercalibration",
                "CorrelationGroupbJES",
                "CorrelationGroupFlavor",
                "CorrelationGroupUncorrelated",
            ],
        },
    })

    # JER
    cfg.x.jer = DotDict.wrap({
        "Jet": {
            "campaign": jer_campaign,
            "version": jer_version,
            "jet_type": jet_type,
        },
    })

    ################################################################################################
    # tau settings
    ################################################################################################

    # tau tagger name
    # (needed by TECConfig below as well as tau selection)
    if run == 2:
        # TODO: still correct? what about 2p5?
        cfg.x.tau_tagger = "DeepTau2017v2p1"
    elif run == 3:
        # https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3?rev=9
        cfg.x.tau_tagger = "DeepTau2018v2p5"
    else:
        assert False

    # tec config
    from columnflow.calibration.cms.tau import TECConfig
    corrector_kwargs = {"wp": "Medium", "wp_VSe": "VVLoose"} if run == 3 else {}
    cfg.x.tec = TECConfig(tagger=cfg.x.tau_tagger, corrector_kwargs=corrector_kwargs)

    # pec config
    from columnflow.calibration.cms.egamma import EGammaCorrectionConfig

    cfg.x.eec = EGammaCorrectionConfig(correction_set="Scale")
    cfg.x.eer = EGammaCorrectionConfig(correction_set="Smearing")

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
    from columnflow.production.cms.electron import ElectronSFConfig
    if run == 2:
        e_postfix = ""
        if year == 2016:
            e_postfix = "preVFP" if campaign.has_tag("preVFP") else "postVFP"
        cfg.x.electron_sf_names = ElectronSFConfig(
            correction="UL-Electron-ID-SF",
            campaign=f"{year}{e_postfix}",
            working_point="wp80iso",
        )
    elif run == 3:
        if year == 2022:
            cmpgn = "2022Re-recoBCD" if campaign.has_tag("preEE") else "2022Re-recoE+PromptFG"
        elif year == 2023:
            cmpgn = "2023PromptC" if campaign.has_tag("preBPix") else "2023PromptD"
        cfg.x.electron_sf_names = ElectronSFConfig(
            correction="Electron-ID-SF",
            campaign=cmpgn,
            working_point="wp80iso",
        )
    else:
        assert False

    ################################################################################################
    # muon settings
    ################################################################################################

    # names of muon correction sets and working points
    # (used in the muon producer)
    from columnflow.production.cms.muon import MuonSFConfig
    if run == 2:
        cfg.x.muon_sf_names = MuonSFConfig(
            correction="NUM_TightRelIso_DEN_TightIDandIPCut",
        )
    elif run == 3:
        cfg.x.muon_sf_names = MuonSFConfig(
            correction="NUM_TightPFIso_DEN_TightID",
        )
    else:
        assert False

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
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer23
        # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer23BPix
        cfg.x.btag_working_points = DotDict.wrap({
            "deepjet": {
                "loose": {"2022": 0.0583, "2022EE": 0.0614, "2023": 0.0479, "2023BPix": 0.048}[btag_key],
                "medium": {"2022": 0.3086, "2022EE": 0.3196, "2023": 0.2431, "2023BPix": 0.2435}[btag_key],
                "tight": {"2022": 0.7183, "2022EE": 0.73, "2023": 0.6553, "2023BPix": 0.6563}[btag_key],
                "xtight": {"2022": 0.8111, "2022EE": 0.8184, "2023": 0.7667, "2023BPix": 0.7671}[btag_key],
                "xxtight": {"2022": 0.9512, "2022EE": 0.9542, "2023": 0.9459, "2023BPix": 0.9483}[btag_key],
            },
            "particleNet": {
                "loose": {"2022": 0.047, "2022EE": 0.0499, "2023": 0.0358, "2023BPix": 0.0359}[btag_key],
                "medium": {"2022": 0.245, "2022EE": 0.2605, "2023": 0.1917, "2023BPix": 0.1919}[btag_key],
                "tight": {"2022": 0.6734, "2022EE": 0.6915, "2023": 0.6172, "2023BPix": 0.6133}[btag_key],
                "xtight": {"2022": 0.7862, "2022EE": 0.8033, "2023": 0.7515, "2023BPix": 0.7544}[btag_key],
                "xxtight": {"2022": 0.961, "2022EE": 0.9664, "2023": 0.9659, "2023BPix": 0.9688}[btag_key],
            },
            "robustParticleTransformer": {
                "loose": {"2022": 0.0849, "2022EE": 0.0897, "2023": 0.0681, "2023BPix": 0.0683}[btag_key],
                "medium": {"2022": 0.4319, "2022EE": 0.451, "2023": 0.3487, "2023BPix": 0.3494}[btag_key],
                "tight": {"2022": 0.8482, "2022EE": 0.8604, "2023": 0.7969, "2023BPix": 0.7994}[btag_key],
                "xtight": {"2022": 0.9151, "2022EE": 0.9234, "2023": 0.8882, "2023BPix": 0.8877}[btag_key],
                "xxtight": {"2022": 0.9874, "2022EE": 0.9893, "2023": 0.9883, "2023BPix": 0.9883}[btag_key],
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
    cfg.x.btag_sf_deepjet = BTagSFConfig(
        correction_set="deepJet_shape",
        jec_sources=cfg.x.btag_sf_jec_sources,
        discriminator="btagDeepFlavB",
    )
    if run == 3:
        cfg.x.btag_sf_pnet = BTagSFConfig(
            correction_set="particleNet_shape",
            jec_sources=cfg.x.btag_sf_jec_sources,
            discriminator="btagPNetB",
        )

    ################################################################################################
    # shifts
    ################################################################################################

    # load jec sources
    with open(os.path.join(thisdir, "jec_sources.yaml"), "r") as f:
        all_jec_sources = yaml.load(f, yaml.Loader)["names"]

    # register shifts
    cfg.add_shift(name="nominal", id=0)

    cfg.add_shift(name="tune_up", id=1, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="tune_down", id=2, type="shape", tags={"disjoint_from_nominal"})

    cfg.add_shift(name="hdamp_up", id=3, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="hdamp_down", id=4, type="shape", tags={"disjoint_from_nominal"})

    cfg.add_shift(name="mtop_up", id=5, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="mtop_down", id=6, type="shape", tags={"disjoint_from_nominal"})

    cfg.add_shift(name="minbias_xs_up", id=7, type="shape")
    cfg.add_shift(name="minbias_xs_down", id=8, type="shape")
    add_shift_aliases(
        cfg,
        "minbias_xs",
        {
            "pu_weight": "pu_weight_{name}",
            "normalized_pu_weight": "normalized_pu_weight_{name}",
        },
    )

    cfg.add_shift(name="top_pt_up", id=9, type="shape")
    cfg.add_shift(name="top_pt_down", id=10, type="shape")
    add_shift_aliases(cfg, "top_pt", {"top_pt_weight": "top_pt_weight_{direction}"})

    for jec_source in cfg.x.jec.Jet.uncertainty_sources:
        idx = all_jec_sources.index(jec_source)
        cfg.add_shift(
            name=f"jec_{jec_source}_up",
            id=5000 + 2 * idx,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        cfg.add_shift(
            name=f"jec_{jec_source}_down",
            id=5001 + 2 * idx,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        add_shift_aliases(
            cfg,
            f"jec_{jec_source}",
            {
                "Jet.pt": "Jet.pt_{name}",
                "Jet.mass": "Jet.mass_{name}",
                f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
                f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
            },
        )
        # TODO: check the JEC de/correlation across years and the interplay with btag weights
        if ("" if jec_source == "Total" else jec_source) in cfg.x.btag_sf_jec_sources:
            add_shift_aliases(
                cfg,
                f"jec_{jec_source}",
                {
                    "normalized_btag_deepjet_weight": "normalized_btag_deepjet_weight_{name}",
                    "normalized_njet_btag_deepjet_weight": "normalized_njet_btag_deepjet_weight_{name}",
                    "normalized_btag_pnet_weight": "normalized_btag_pnet_weight_{name}",
                    "normalized_njet_btag_pnet_weight": "normalized_njet_btag_pnet_weight_{name}",
                },
            )

    cfg.add_shift(name="jer_up", id=6000, type="shape", tags={"jer"})
    cfg.add_shift(name="jer_down", id=6001, type="shape", tags={"jer"})
    add_shift_aliases(
        cfg,
        "jer",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
            f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
        },
    )

    for i, (match, dm) in enumerate(itertools.product(["jet", "e"], [0, 1, 10, 11])):
        cfg.add_shift(name=f"tec_{match}_dm{dm}_up", id=20 + 2 * i, type="shape", tags={"tec"})
        cfg.add_shift(name=f"tec_{match}_dm{dm}_down", id=21 + 2 * i, type="shape", tags={"tec"})
        add_shift_aliases(
            cfg,
            f"tec_{match}_dm{dm}",
            {
                "Tau.pt": "Tau.pt_{name}",
                "Tau.mass": "Tau.mass_{name}",
                f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
                f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
            },
        )

    # start at id=50
    cfg.x.tau_unc_names = [
        "jet_dm0", "jet_dm1", "jet_dm10",
        "e_barrel", "e_endcap",
        "mu_0p0To0p4", "mu_0p4To0p8", "mu_0p8To1p2", "mu_1p2To1p7", "mu_1p7To2p3",
    ]
    for i, unc in enumerate(cfg.x.tau_unc_names):
        cfg.add_shift(name=f"tau_{unc}_up", id=50 + 2 * i, type="shape")
        cfg.add_shift(name=f"tau_{unc}_down", id=51 + 2 * i, type="shape")
        add_shift_aliases(cfg, f"tau_{unc}", {"tau_weight": f"tau_weight_{unc}_" + "{direction}"})

    cfg.add_shift(name="tautau_trigger_up", id=80, type="shape")
    cfg.add_shift(name="tautau_trigger_down", id=81, type="shape")
    add_shift_aliases(cfg, "tautau_trigger", {"tau_trigger_weight": "tau_trigger_weight_tautau_{direction}"})
    cfg.add_shift(name="etau_trigger_up", id=82, type="shape")
    cfg.add_shift(name="etau_trigger_down", id=83, type="shape")
    add_shift_aliases(cfg, "etau_trigger", {"tau_trigger_weight": "tau_trigger_weight_etau_{direction}"})
    cfg.add_shift(name="mutau_trigger_up", id=84, type="shape")
    cfg.add_shift(name="mutau_trigger_down", id=85, type="shape")
    add_shift_aliases(cfg, "mutau_trigger", {"tau_trigger_weight": "tau_trigger_weight_mutau_{direction}"})
    # no uncertainty for di-tau VBF trigger existing yet
    # cfg.add_shift(name="mutau_trigger_up", id=86, type="shape")
    # cfg.add_shift(name="tautauvbf_trigger_down", id=86, type="shape")
    # add_shift_aliases(cfg, "tautauvbf_trigger", {"tau_trigger_weight": "tau_trigger_weight_tautauvbf_{direction}"})

    cfg.add_shift(name="e_up", id=90, type="shape")
    cfg.add_shift(name="e_down", id=91, type="shape")
    add_shift_aliases(cfg, "e", {"electron_weight": "electron_weight_{direction}"})

    # electron shifts
    # TODO: energy corrections are currently only available for 2022 (Jan 2025)
    #       include them when available
    if run == 3 and year == 2022:
        cfg.add_shift(name="eec_up", id=92, type="shape", tags={"eec"})
        cfg.add_shift(name="eec_down", id=93, type="shape", tags={"eec"})
        add_shift_aliases(
            cfg,
            "eec",
            {
                "Electron.pt": "Electron.pt_scale_{direction}",
            },
        )

        cfg.add_shift(name="eer_up", id=94, type="shape", tags={"eer"})
        cfg.add_shift(name="eer_down", id=95, type="shape", tags={"eer"})
        add_shift_aliases(
            cfg,
            "eer",
            {
                "Electron.pt": "Electron.pt_res_{direction}",
            },
        )

    cfg.add_shift(name="mu_up", id=100, type="shape")
    cfg.add_shift(name="mu_down", id=101, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})

    cfg.x.btag_unc_names = [
        "hf", "lf",
        f"hfstats1_{year}", f"hfstats2_{year}",
        f"lfstats1_{year}", f"lfstats2_{year}",
        "cferr1", "cferr2",
    ]
    for i, unc in enumerate(cfg.x.btag_unc_names):
        cfg.add_shift(name=f"btag_{unc}_up", id=110 + 2 * i, type="shape")
        cfg.add_shift(name=f"btag_{unc}_down", id=111 + 2 * i, type="shape")
        add_shift_aliases(
            cfg,
            f"btag_{unc}",
            {
                "normalized_btag_deepjet_weight": f"normalized_btag_deepjet_weight_{unc}_" + "{direction}",
                "normalized_njet_btag_deepjet_weight": f"normalized_njet_btag_deepjet_weight_{unc}_" + "{direction}",
                # TODO: pnet here, or is this another shift? probably the latter
            },
        )

    cfg.add_shift(name="pdf_up", id=130, type="shape", tags={"lhe_weight"})
    cfg.add_shift(name="pdf_down", id=131, type="shape", tags={"lhe_weight"})
    add_shift_aliases(
        cfg,
        "pdf",
        {
            "pdf_weight": "pdf_weight_{direction}",
            "normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )

    cfg.add_shift(name="murmuf_up", id=140, type="shape", tags={"lhe_weight"})
    cfg.add_shift(name="murmuf_down", id=141, type="shape", tags={"lhe_weight"})
    add_shift_aliases(
        cfg,
        "murmuf",
        {
            "murmuf_weight": "murmuf_weight_{direction}",
            "normalized_murmuf_weight": "normalized_murmuf_weight_{direction}",
        },
    )


    # # register shifts
    # cfg.add_shift(name="nominal", id=0)

    # # tune shifts are covered by dedicated, varied datasets, so tag the shift as "disjoint_from_nominal"
    # # (this is currently used to decide whether ML evaluations are done on the full shifted dataset)
    # cfg.add_shift(name="tune_up", id=1, type="shape", tags={"disjoint_from_nominal"})
    # cfg.add_shift(name="tune_down", id=2, type="shape", tags={"disjoint_from_nominal"})

    # # fake jet energy correction shift, with aliases flaged as "selection_dependent", i.e. the aliases
    # # affect columns that might change the output of the event selection
    # cfg.add_shift(name="jec_up", id=20, type="shape")
    # cfg.add_shift(name="jec_down", id=21, type="shape")
    # add_shift_aliases(
    #     cfg,
    #     "jec",
    #     {
    #         "Jet.pt": "Jet.pt_{name}",
    #         "Jet.mass": "Jet.mass_{name}",
    #         "MET.pt": "MET.pt_{name}",
    #         "MET.phi": "MET.phi_{name}",
    #     },
    # )

    # # event weights due to muon scale factors
    # cfg.add_shift(name="mu_up", id=10, type="shape")
    # cfg.add_shift(name="mu_down", id=11, type="shape")
    # add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})

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
        json_mirror = "/afs/cern.ch/user/m/mrieger/public/mirrors/jsonpog-integration-7439b936"
    elif run == 3:
        json_pog_era = f"{year}_Summer{year2}{campaign.x.postfix}"
        json_mirror = "/afs/cern.ch/user/m/mrieger/public/mirrors/jsonpog-integration-7439b936"
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

                    # TODO: electron (and photon) energy corrections and smearing are only available for 2022
        #       include them when available
        if year == 2022:
            # electron energy correction and smearing
            add_external("electron_ss", (f"{json_mirror}/POG/EGM/{json_pog_era}/electronSS.json.gz", "v1"))
            
            # tau energy correction and scale factors
            # TODO: remove tag pog mirror once integrated centrally
            json_mirror_tau_pog = "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-taupog"
            tau_pog_era = f"{year}_{'pre' if campaign.has_tag('preEE') else 'post'}EE"
            add_external("tau_sf", (f"{json_mirror_tau_pog}/POG/TAU/{tau_pog_era}/tau_DeepTau2018v2p5_{tau_pog_era}.json.gz", "v1"))  # noqa
    else:
        assert False

    



    ################################################################################################
    # reductions
    ################################################################################################

    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0

    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # general event info, mandatory for reading files with coffea
            ColumnCollection.MANDATORY_COFFEA,  # additional columns can be added as strings, similar to object info
            # object info
            "Jet.*", 
            "FatJet.*",          
            "MET.pt", "MET.phi", "MET.significance", "MET.covXX", "MET.covXY", "MET.covYY",
            "Muon.*",
            "Electron.*",
            "Tau.*",
            "PV.npvs",
            "PFCandidate.*",
            "GenJet.*", "GenJetAK8.*", "GenVisTau.*",
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

    # # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    # get_shifts = functools.partial(get_shifts_from_sources, cfg)

    ################################################################################################
    # weights
    ################################################################################################

    # configurations for all possible event weight columns as keys in an OrderedDict,
    # mapped to shift instances they depend on
    # (this info is used by weight producers)
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        "normalization_weight_inclusive": [],
        "pdf_weight": get_shifts("pdf"),
        "murmuf_weight": get_shifts("murmuf"),
        "normalized_pu_weight": get_shifts("minbias_xs"),
        # TODO: enable again once we have btag cuts
        # "normalized_njet_btag_deepjet_weight": get_shifts(*(f"btag_{unc}" for unc in cfg.x.btag_unc_names)),
        "electron_weight": get_shifts("e"),
        "muon_weight": get_shifts("mu"),
        "tau_weight": get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names)),
        "tau_trigger_weight": get_shifts("etau_trigger", "mutau_trigger", "tautau_trigger"),
    })

    # define per-dataset event weights
    for dataset in cfg.datasets:
        if dataset.has_tag("ttbar"):
            dataset.x.event_weights = {"top_pt_weight": get_shifts("top_pt")}

    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    cfg.x.versions = {
        # "cf.CalibrateEvents": "prod1",
        # "cf.SelectEvents": (lambda cls, inst, params: "prod1" if params.get("selector") == "default" else "dev1"),
        # ...
    }

    ################################################################################################
    # external configs: channels, categories, met filters, triggers, variables, hist hooks
    ################################################################################################

    # channels
    cfg.add_channel(name="etau", id=1, label=r"$e\tau_{h}$")
    cfg.add_channel(name="mutau", id=2, label=r"$\mu\tau_{h}$")
    cfg.add_channel(name="tautau", id=3, label=r"$\tau_{h}\tau_{h}$")
    cfg.add_channel(name="ee", id=4, label=r"$ee$")
    cfg.add_channel(name="mumu", id=5, label=r"$\mu\mu$")
    cfg.add_channel(name="emu", id=6, label=r"$e\mu$")

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

    # add triggers
    if year == 2016:
        from hhh4b2tau.config.triggers import add_triggers_2016
        add_triggers_2016(cfg)
    elif year == 2017:
        from hhh4b2tau.config.triggers import add_triggers_2017
        add_triggers_2017(cfg)
    elif year == 2018:
        from hhh4b2tau.config.triggers import add_triggers_2018
        add_triggers_2018(cfg)
    elif year == 2022:
        from hhh4b2tau.config.triggers import add_triggers_2022
        add_triggers_2022(cfg)
    elif year == 2023:
        from hhh4b2tau.config.triggers import add_triggers_2023
        add_triggers_2023(cfg)
    else:
        raise False

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