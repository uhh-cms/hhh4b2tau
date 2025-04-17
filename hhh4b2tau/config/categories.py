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
            # if cat.name != "os"  # os is the default
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
    add_category(config, name="4j", id=111, selection="cat_4j", label="4 jets")

    add_category(config, name="1btag", id=112, selection="cat_1btag", label="1 b-tag(s)")
    add_category(config, name="2btag", id=113, selection="cat_2btag", label="2 b-tag(s)")
    add_category(config, name="3btag", id=114, selection="cat_3btag", label="3 b-tag(s)")

    add_category(config, name="h3_mass", id=120, selection="cat_h3_mass", label=r"$m_{H3}<125 GeV$", tags={"h3_mass"})
    add_category(config, name="h3_mass_orth", id=121, selection="cat_h3_mass_orth", label=r"$m_{H3}>125 GeV$", tags={"h3_mass_orth"})

    # add_category(config, name="leps_cos", id=124, selection="cat_leps_cos", label=r"$\cos(\delta) > -0.25$", tags={"leps_cos"})
    # add_category(config, name="leps_cos_orth", id=125, selection="cat_leps_cos_orth", label=r"$\cos(\delta) < -0.25$", tags={"leps_cos_orth"})

    add_category(config, name="leps_dr", id=126, selection="cat_leps_dr", label=r"$\Delta R_{ll} < 2.4$", tags={"leps_dr"})
    add_category(config, name="leps_dr_orth", id=127, selection="cat_leps_dr_orth", label=r"$\Delta R_{ll} > 2.4$", tags={"leps_dr_orth"})

    #
    # build groups
    #



    # main analysis categories
    main_categories = {
        # channels first
        "channel": [
            # config.get_category("etau"), 
            # config.get_category("tautau"), 
            config.get_category("mutau"),
        ],
        # kinematic regions in the middle (to be extended)
        "kin": [
            config.get_category("incl"), 
            # config.get_category("4j"), 
            # config.get_category("1btag"), 
            config.get_category("2btag"), 
            config.get_category("3btag"),
        ],
        "h3_mass": [
            config.get_category("h3_mass"), config.get_category("h3_mass_orth"),
        ],
        # "leps_cos": [
        #     config.get_category("leps_cos"), config.get_category("leps_cos_orth"),
        # ],
        "leps_dr": [
            config.get_category("leps_dr"), config.get_category("leps_dr_orth"),
        ],
        # qcd regions last
        "sign": [
            config.get_category("os"),
            #  config.get_category("ss"),
                ],
        "tau2": [
            config.get_category("iso"), 
            #  config.get_category("noniso"),
                ],
    }
    create_category_combinations(
        config=config,
        categories=main_categories,
        name_fn=name_fn,
        kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=True),
    )

    # control categories
    control_categories = {
        # channels first
        "channel": [
            config.get_category("ee"), config.get_category("mumu"), config.get_category("emu"),
        ],
        # kinematic regions in the middle (to be extended)
        "kin": [
            config.get_category("incl"), 
            # config.get_category("2j")
            ],
        # relative sign last
        "sign": [config.get_category("os")],
    }
    # commentated out for now for speed
    # create_category_combinations(
    #     config=config,
    #     categories=control_categories,
    #     name_fn=name_fn,
    #     kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=False),
    # )



def add_categories_incl_only(config):
    add_category(config, name="incl", id=100, selection="cat_incl", label="inclusive")

def add_analysis_categories(config):

    add_category(config, name="h3_mass", id=120, selection="cat_h3_mass", label=r"$m_{H3}<125 GeV$", tags={"h3_mass"})
    add_category(config, name="h3_mass_orth", id=121, selection="cat_h3_mass_orth", label=r"$m_{H3}>125 GeV$", tags={"h3_mass_orth"})

    # add_category(config, name="leps_cos", id=124, selection="cat_leps_cos", label=r"$\cos(\delta) > -0.25$", tags={"leps_cos"})
    # add_category(config, name="leps_cos_orth", id=125, selection="cat_leps_cos_orth", label=r"$\cos(\delta) < -0.25$", tags={"leps_cos_orth"})

    add_category(config, name="leps_dr", id=126, selection="cat_leps_dr", label=r"$\Delta R_{ll} < 2.4$", tags={"leps_dr"})
    add_category(config, name="leps_dr_orth", id=127, selection="cat_leps_dr_orth", label=r"$\Delta R_{ll} > 2.4$", tags={"leps_dr_orth"})

    add_category(config, name="standard", id=130, selection="cat_standard", label=r"$\mu_\tau$,\n os,\n iso,\n > 2 b-jets")

    analysis_categories = {
        "standard": [config.get_category("standard")],
        "h3_mass": [
            config.get_category("h3_mass"), config.get_category("h3_mass_orth"),
        ],
        # "leps_cos": [
        #     config.get_category("leps_cos"), config.get_category("leps_cos_orth"),
        # ],
        "leps_dr": [
            config.get_category("leps_dr"), config.get_category("leps_dr_orth"),
        ],
    }
    create_category_combinations(
        config=config,
        categories=analysis_categories,
        name_fn=name_fn,
        kwargs_fn=functools.partial(kwargs_fn, add_qcd_group=True),
    )