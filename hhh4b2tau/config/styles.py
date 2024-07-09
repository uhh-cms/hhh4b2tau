# coding: utf-8

"""
Style definitions.
"""

import order as od
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
def stylize_processes(config: od.Config) -> None:
    """
    Adds process colors and adjust labels.
    """
# TODO

