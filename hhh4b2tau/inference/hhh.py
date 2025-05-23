# coding: utf-8

"""
Example inference model.
"""

from columnflow.inference import inference_model, ParameterType, ParameterTransformation


@inference_model
def hhh(self):

    #
    # categories
    #

    variable = "m3j2l"
    # variable = "chi2"
    kin_cat =["3b0j", "3b1j", "4b"]
    data = ["TT", "ggHH_kl1_kt1_13p6TeV_hbbhtt", "ttH"]

    self.add_category(
        "incl",
        config_category="incl",
        config_variable=variable,
        data_from_processes=data,
        mc_stats=True,
    )


    self.add_category(
        f"mutau__incl__os__iso" ,
        config_category=f"mutau__incl__os__iso" ,
        config_variable=variable,
        data_from_processes=data,
        mc_stats=True,
    )
    self.add_category(
        f"mutau__incl__os__iso__var_cuts" ,
        config_category=f"mutau__incl__os__iso__var_cuts" ,
        config_variable=variable,
        data_from_processes=data,
        mc_stats=True,
    )
    # self.add_category(
    #     f"mutau__incl__os__iso__var_cuts__3j2l" ,
    #     config_category=f"mutau__incl__os__iso__var_cuts__3j2l" ,
    #     config_variable=variable,
    #     data_from_processes=data,
    #     mc_stats=True,
    # )

    #
    # processes
    #

    self.add_process(
        "TT",
        config_process="tt",
        config_mc_datasets=["tt_sl_powheg", "tt_dl_powheg", "tt_fh_powheg"],
    )

    self.add_process(
        "ggHH_kl1_kt1_13p6TeV_hbbhtt",
        config_process="hh_ggf_hbb_htt_kl1_kt1",
    )
    self.add_process(
        "ttH",
        config_process="tth",
    )

    self.add_process(
        "ggHHH_c30_d40_13p6TeV_hbbhbbhtt",
        config_process="hhh_4b2tau_c30_d40",
        is_signal=True,
    )

    #
    # parameters
    #

    # groups
    self.add_parameter_group("experiment")
    self.add_parameter_group("theory")

    # lumi
    lumi = self.config_inst.x.luminosity
    for unc_name in lumi.uncertainties:
        self.add_parameter(
            unc_name,
            type=ParameterType.rate_gauss,
            effect=lumi.get(names=unc_name, direction=("down", "up"), factor=True),
            transformations=[ParameterTransformation.symmetrize],
            group=["experiment"],
        )

    self.add_parameter(
        "BR_hbb",
        type=ParameterType.rate_gauss,
        process=["*_hbb", "*_hbbhtt"],
        effect=(0.9874, 1.0124),
        group=["theory", "signal_norm_xsbr"],
    )
    self.add_parameter(
        "BR_htt",
        type=ParameterType.rate_gauss,
        process=["*_htt", "*_hbbhtt"],
        effect=(0.9837, 1.0165),
        group=["theory", "signal_norm_xsbr"],
    )
    self.add_parameter(
        "pdf_gg",  # contains alpha_s
        type=ParameterType.rate_gauss,
        process="TT",
        effect=1.042,
        group=["theory"],
    )
    self.add_parameter(
        "pdf_Higgs_ggHH",  # contains alpha_s
        type=ParameterType.rate_gauss,
        process="ggHH_*",
        effect=1.023,
        group=["theory", "signal_norm_xs", "signal_norm_xsbr"],
    )
    self.add_parameter(
        "QCDscale_ttbar",
        type=ParameterType.rate_gauss,
        process="TT",
        effect=(0.965, 1.024),
        group=["theory"],
    )

    self.add_parameter(
        "xs_hhh_total",
        type=ParameterType.rate_gauss,
        process="ggHHH_*",
        effect=1.2,
        group=["theory"],
    )
    
