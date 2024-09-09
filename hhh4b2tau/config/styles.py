# coding: utf-8

"""
Style definitions.
"""

import order as od

from columnflow.util import DotDict

# recommended colors by cms for increased accessibility
cms_six_color_scheme = [
    "#5790fc", 
    "#f89c20", 
    "#e42536", 
    "#964a8b", 
    "#9c9ca1", 
    "#7a21dd",
    ]
cms_eight_color_scheme = [
    "#1845fb", 
    "#ff5e02", 
    "#c91f16", 
    "#c849a9", 
    "#adad7d", 
    "#86c8dd", 
    "#578dff", 
    "#656364",
    ]
cms_ten_color_scheme = [
    "#3f90da", # blue
    "#ffa90e", # yellow
    "#bd1f01", # red
    "#94a4a2", # grey
    "#832db6", # purple
    "#a96b59", # brown
    "#e76300", # orange
    "#b9ac70", # green-brown
    "#717581", # dark grey
    "#92dadd", # cyan
    ]

coupling_with_colors = [
# (c3, d4, color)
    (0, 0, "#000000"),
    (0, 99, "#3f90da"),
    (0, 'minus1', "#ffa90e"),
    (19, 19, "#bd1f01"),
    (1, 0, "#94a4a2"),
    (1, 2, "#832db6"),
    (2, 'minus1', "#a96b59"),
    (4, 9, "#e76300"),
    ('minus1', 0, "#b9ac70"),
    ('minus1', 'minus1', "#717581"),
    ('minus1p5', 'minus0p5', "#92dadd"),
]

# recommended cms colors
colors = DotDict(
    bright_blue="#3f90da",
    dark_blue="#011c87",
    purple="#832db6",
    aubergine="#964a8b",
    yellow="#f7c331",
    bright_orange="#ffa90e",
    dark_orange="#e76300",
    red="#bd1f01",
    teal="#92dadd",
    grey="#94a4a2",
    brown="#a96b59",
)

def stylize_processes(config: od.Config) -> None:
    """
    Adds process colors and adjust labels.
    """
    cfg = config
    cfg.x.colors = colors
    
    for kl in ["0", "1", "2p45", "5"]:
        if (p := config.get_process(f"hh_ggf_hbb_htt_kl{kl}_kt1", default=None)):
            p.color1 = cfg.x.colors.bright_blue
            p.label = (
                r"$HH_{ggf} \rightarrow bb\tau\tau$ __SCALE__"
                "\n"
                rf"($\kappa_{{\lambda}}$={kl.replace('p', '.')},$\kappa_{{t}}$=1)"
            )
            p.scale="stack"
            p.unstack=True


    for c3, d4, color in coupling_with_colors:
        if (p := config.get_process(f"hhh_4b2tau_c3{c3}_d4{d4}", default=None)):
            p.color1 = color
            p.unstack=True
            p.scale="stack"


    if (p := config.get_process("hh_vbf_hbb_htt_kv1_k2v1_kl1", default=None)):
        p.color1 = cfg.x.colors.dark_blue
        p.label = (
            r"$HH_{vbf} \rightarrow bb\tau\tau$ __SCALE__"
            "\n"
            r"($\kappa_{\lambda}$=1,$\kappa_{V}$=1,$\kappa_{2V}$=1)"
        )
        p.scale="stack"
        p.unstack=True

    if (p := config.get_process("tth", default=None)):
        p.color1 = cfg.x.colors.purple
        p.unstack=True
        p.label = r"$t\bar{t}H$"
        p.scale="stack"

    if (p := config.get_process("tt", default=None)):
        p.color1 = cfg.x.colors.bright_orange
        p.label = r"$t\bar{t}$"
    if (p := config.get_process("tt_fh", default=None)):
        p.color1 = cfg.x.colors.brown
        p.label = r"$t\bar{t} (FH)$"
    if (p := config.get_process("tt_sl", default=None)):
        p.color1 = cfg.x.colors.bright_orange
        p.label = r"$t\bar{t} (SL)$"
    if (p := config.get_process("tt_dl", default=None)):
        p.color1 = cfg.x.colors.teal
        p.label = r"$t\bar{t} (DL)$"
    

    if (p := config.get_process("st", default=None)):
        p.color1 = cfg.x.colors.aubergine

    if (p := config.get_process("dy", default=None)):
        p.color1 = cfg.x.colors.dark_orange

    if (p := config.get_process("vv", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("vvv", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("multiboson", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("w", default=None)):
        p.color1 = cfg.x.colors.teal
        p.label = "W"

    if (p := config.get_process("z", default=None)):
        p.color1 = cfg.x.colors.brown
        p.label = "Z"

    if (p := config.get_process("v", default=None)):
        p.color1 = cfg.x.colors.teal

    if (p := config.get_process("ewk", default=None)):
        p.color1 = cfg.x.colors.brown

    if (p := config.get_process("ttv", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$t\bar{t} + V$"

    if (p := config.get_process("ttvv", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$t\bar{t} + VV$"

    if (p := config.get_process("tt_multiboson", default=None)):
        p.color1 = cfg.x.colors.grey

    if (p := config.get_process("qcd", default=None)):
        p.color1 = cfg.x.colors.red
