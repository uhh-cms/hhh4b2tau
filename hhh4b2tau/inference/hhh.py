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
    kin_cat =[
        "incl",
        "3b0j", 
        "3b1j", 
        "4b",
        "1b",
        "2b",
        "3b",
    ]
    
    # data = ["TT", "ggHH_kl1_kt1_13p6TeV_hbbhtt", "ttH"]
    data = ["TT", "ggHH_kl1_kt1_13p6TeV_hbbhtt", "ttH", "DY"]
    categories = [
        'incl', 
        # "mutau__incl__os__iso" , 
        # "mutau__incl__os__iso__var_cuts" , 
        # "mutau__1b__os__iso" , 
        # "mutau__1b__os__iso__var_cuts" , 
        # "mutau__2b__os__iso" , 
        # "mutau__2b__os__iso__var_cuts" , 
        # "mutau__3b__os__iso" , 
        # "mutau__3b__os__iso__var_cuts" , 
        # "mutau__4b__os__iso" , 
        # "mutau__4b__os__iso__var_cuts" , 
        # "mutau__5b__os__iso" , 
        # "mutau__5b__os__iso__var_cuts" , 
        # "mutau__3b0j__os__iso" , 
        # "mutau__3b0j__os__iso__var_cuts" ,
        # "mutau__3b1j__os__iso" , 
        # "mutau__3b1j__os__iso__var_cuts" ,
        # "mutau__incl__os__iso__mh3" , 
        # "mutau__incl__os__iso__mh3__leps_dr" , 
        # "mutau__j1_pt__os__iso__var_cuts" , 
        # "mutau__incl__os__iso__var_cuts__bb1_dr" ,
        # "mutau__j1_pt__os__iso__var_cuts__bb1_dr" ,
    ]

    for cat in categories:

        self.add_category(
            cat,
            config_category=cat,
            config_variable=variable,
            data_from_processes=data,
            mc_stats=True,
        )

    # self.add_category(
    #     "incl",
    #     config_category="incl",
    #     config_variable=variable,
    #     data_from_processes=data,
    #     mc_stats=True,
    # )


    # self.add_category(
    #     "mutau__incl__os__iso" ,
    #     config_category="mutau__incl__os__iso" ,
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
        "DY",
        config_process="dy",
        config_mc_datasets=[
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
        ],
    )

    self.add_process(
        "ggHHH_c30_d40_13p6TeV_hbbhbbhtt",
        config_process="hhh_4b2tau_c30_d40",
        is_signal=True,
    )

    # from hhh4b2tau.hist_hooks.morphing import morphing_coupling_combinations
    # self.add_parameter_group("switches")
    # for c3,d4 in morphing_coupling_combinations:
    #     proc_name = "ggHHH_c3{c3}_d4{d4}_13p6TeV_hbbhbbhtt".format(
    #                   c3=str(c3).replace("-", "m").replace(".", "p"),
    #                   d4=str(d4).replace("-", "m").replace(".", "p"),
    #                   )
    #     self.add_process(
    #         proc_name,
    #         config_process="hhh_4b2tau_c3{c3}_d4{d4}".format(
    #                   c3=str(c3).replace("-", "m").replace(".", "p"),
    #                   d4=str(d4).replace("-", "m").replace(".", "p"),
    #                   ),
    #         is_signal=True,
    #     )

    #     self.add_parameter(
    #         f"switch_c3{c3}_d4{d4}",
    #         type=ParameterType.rate_unconstrained,
    #         process=proc_name,
    #         effect=1 if (c3 == 0 and d4 == 0) else 0,
    #         group=["switches"],
    #     )



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

    
    
