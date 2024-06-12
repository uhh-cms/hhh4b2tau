# coding: utf-8

"""
Style definitions.
"""

import order as od


def stylize_processes(config: od.Config) -> None:
    """
    Adds process colors and adjust labels.
    """
    # if config.has_process("c3_0_d4_0"):
    #     config.processes.n.c3_0_d4_0.color1 = (67, 118, 201)

    # if config.has_process("h"):
    #     config.processes.n.h.color1 = (65, 180, 219)

    # if config.has_process("tt"):
    #     config.processes.n.tt.color1 = (244, 182, 66)

    # if config.has_process("st"):
    #     config.processes.n.st.color1 = (244, 93, 66)

    # if config.has_process("dy"):
    #     config.processes.n.dy.color1 = (68, 186, 104)

    # if config.has_process("vv"):
    #     config.processes.n.vv.color1 = (2, 24, 140)

    # if config.has_process("qcd"):
    #     config.processes.n.qcd.color1 = (242, 149, 99)
    comparison = [(0, 0), 
                  (0, 99), 
                  (1, 0), 
                  (19, 19)
                  ]
    colors =[(67, 118, 201), 
             (65, 180, 219), 
             (244, 182, 66), 
             (244, 93, 66), 
             (68, 186, 104), 
             (2, 24, 140), 
             ]
    # accesibility colors provided by cms
    six_colors = [
        "#5790fc", 
        "#f89c20", 
        "#e42536", 
        "#964a8b", 
        "#9c9ca1", 
        "#7a21dd"
        ]
    
    ten_colors = [
        "#3f90da", 
        "#ffa90e", 
        "#bd1f01", 
        "#94a4a2", 
        "#832db6", 
        "#a96b59", 
        "#e76300", 
        "#b9ac70", 
        "#717581", 
        "#92dadd"
        ]
