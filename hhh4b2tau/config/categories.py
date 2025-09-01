# coding: utf-8

"""
Definition of categories.
"""

import functools

import order as od

from columnflow.config_util import add_category, create_category_combinations

def name_fn(categories):
    return "__".join(cat.name for cat in categories.values() if cat)

def kwargs_fn(categories, add_qcd_group=True):
    # build auxiliary information
    aux = {}
    if add_qcd_group:
        aux["qcd_group"] = name_fn({
            name: cat for name, cat in categories.items()
            if name not in {"sign", "tau2"}
        })
    # return the desired kwargs
    return {
        # just increment the category id
        # NOTE: for this to be deterministic, the order of the categories must no change!
        "id": "+",
        # join all tags
        "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
        # auxiliary information
        "aux": aux,
        # label
        "label": "\n".join([
            cat.label or cat.name
            for cat in categories.values()
            if cat.name != "os"  # os is the default
            if cat.name != "iso"  # iso is the default
            if cat.name != "incl"  # incl is the default
        ]) or None,
    }

def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    # lepton channels
    add_category(config, name="etau", id=1, selection="cat_etau", label=config.channels.n.etau.label)
    add_category(config, name="mutau", id=2, selection="cat_mutau", label=config.channels.n.mutau.label)
    add_category(config, name="tautau", id=3, selection="cat_tautau", label=config.channels.n.tautau.label)
    
    add_category(config, name="ee", id=4, selection="cat_ee", label=config.channels.n.ee.label)
    add_category(config, name="mumu", id=5, selection="cat_mumu", label=config.channels.n.mumu.label)
    add_category(config, name="emu", id=6, selection="cat_emu", label=config.channels.n.emu.label)

    # qcd regions
    add_category(config, name="os", id=10, selection="cat_os", label="OS", tags={"os"})
    add_category(config, name="ss", id=11, selection="cat_ss", label="SS", tags={"ss"})
    add_category(config, name="iso", id=12, selection="cat_iso", label=r"$\tau_{h,2}$ iso", tags={"iso"})
    add_category(config, name="noniso", id=13, selection="cat_noniso", label=r"$\tau_{h,2}$ non-iso", tags={"noniso"})  # noqa: E501

    # kinematic categories
    add_category(config, name="incl", id=100, selection="cat_incl", label="inclusive")

    # add_category(config, name="3b0j", id=115, selection="cat_3b0j", label=r"$=3$ b-tags, 0 jets")
    # add_category(config, name="3b1j", id=116, selection="cat_3b1j", label=r"$=3$ b-tags, $\geq 1$ jets")
    
    # add_category(config, name="1b", id=118, selection="cat_1b", label=r"$\geq 1$ b-tags")
    # add_category(config, name="2b", id=119, selection="cat_2b", label=r"$\geq 2$ b-tags")
    # add_category(config, name="3b", id=120, selection="cat_3b", label=r"$\geq 3$ b-tags")
    # add_category(config, name="4b", id=117, selection="cat_4b", label=r"$\geq 4$ b-tags")
    # add_category(config, name="5b", id=180, selection="cat_5b", label=r"$\geq 5$ b-tags")

    add_category(config, name="mh3", id=120, selection="cat_mh3", label=r"$m_{H3}<125 \ GeV$", tags={"h3_mass"})
    # add_category(config, name="mh3_orth", id=121, selection="cat_mh3_orth", label=r"$m_{H3}>125 \ GeV$", tags={"h3_mass_orth"})

    add_category(config, name="leps_dr", id=126, selection="cat_leps_dr", label=r"$\Delta R_{ll} < 2.0$")
    # add_category(config, name="leps_dr_orth", id=127, selection="cat_leps_dr_orth", label="")

    add_category(config, name="bb1_dr", id=142, selection="cat_bb1_dr", label=r"$\Delta R_{jj1} < 2.4$", tags={"bb1_dr"})
    # add_category(config, name="bb1_dr_orth", id=143, selection="cat_bb1_dr_orth", label=r"$\Delta R_{jj1} > 2.4$", tags={"bb1_dr_orth"})

    add_category(config, name="j1_pt", id=130, selection="cat_j1_pt", label=r"$j_1 p_T >100 \ GeV$")
    # add_category(config, name="j1_pt_orth", id=131, selection="cat_j1_pt_orth", label=r"$j_1 p_T<100 \ GeV$")

    # add_category(config, name="var_cuts", id=140, selection="cat_var_cuts", label=r"cuts applied")
    # add_category(config, name="var_cuts_orth", id=141, selection="cat_var_cuts_orth", label="")

    # add_category(config, name="chi2", id=132, selection="cat_chi2", label="$\chi^2 \leq 0.31$")

    #
    # build groups
    #



    # main analysis categories
    main_categories = {
        # channels first
        "channel": [
            config.get_category("etau"), 
            config.get_category("tautau"), 
            config.get_category("mutau"),
        ],
        # "kin": [
        #     config.get_category("incl"),
        #     # config.get_category("j1_pt"),

        #     # config.get_category("3b0j"), 
        #     # config.get_category("3b1j"), 
            
        #     # config.get_category("1b"), 
        #     # config.get_category("2b"), 
        #     # config.get_category("3b"), 
        #     # config.get_category("4b"), 
        #     # config.get_category("5b"), 
        # ],
        # qcd regions last
        "sign": [
            config.get_category("os"),
            config.get_category("ss"),
                ],
        "tau2": [
            config.get_category("iso"), 
            config.get_category("noniso"),
                ],
        # "var_cuts": [
        #     config.get_category("var_cuts"), 
        #     config.get_category("var_cuts_orth")
        # ],
        "h3_mass": [
            config.get_category("mh3"), 
            # config.get_category("mh3_orth"),
        ],
        "leps_dr": [
            config.get_category("leps_dr"), 
            # config.get_category("leps_dr_orth"),
            # config.get_category("leps_dr_harsh"), 
        ],
        "j1_pt":[
            config.get_category("j1_pt"), 
            # config.get_category("j1_pt_orth"),
        ],
        "bb1_dr": [
            config.get_category("bb1_dr"), 
            # config.get_category("bb1_dr_orth"),
        ],

    }
    create_category_combinations(
        config=config,
        categories=main_categories,
        name_fn=name_fn,
        kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=True),
    )

    # # control categories
    # control_categories = {
    # #     # channels first
    #     "channel": [
    #         config.get_category("ee"), config.get_category("mumu"), config.get_category("emu"),
    #     ],
    # #     # kinematic regions in the middle (to be extended)
    #     "kin": [
    #         config.get_category("incl"), 
    #         # config.get_category("2j")
    #         ],
    #     # relative sign last
    #     "sign": [config.get_category("os")],
    # }
    # # commentated out for now for speed
    # create_category_combinations(
    #     config=config,
    #     categories=control_categories,
    #     name_fn=name_fn,
    #     kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=False),
    # )



def add_categories_incl_only(config):
    add_category(config, name="incl", id=100, selection="cat_incl", label="inclusive")