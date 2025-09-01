# coding: utf-8

"""
Examples for custom plot functions.
"""

from __future__ import annotations

from collections import defaultdict, OrderedDict

import law

from columnflow.util import maybe_import, try_float
from columnflow.plotting.plot_all import (
    plot_all, draw_error_bands, draw_stack, 
    draw_hist, draw_profile, draw_errorbars
)
# from columnflow.plotting.plot_util import (
#     # prepare_stack_plot_config,
#     prepare_style_config,
#     remove_residual_axis,
#     # apply_variable_settings,
#     new_apply_variable_settings,
#     apply_process_settings,
#     apply_density_to_hists,
#     # apply_density, # use this for new cf version
# )
from columnflow.plotting.new_plot_util import (
    # prepare_stack_plot_config,
    prepare_style_config,
    remove_residual_axis,
    apply_variable_settings,
    apply_process_settings,
    apply_density,
)

from columnflow.plotting.plot_util import (
    get_position,
    get_cms_label,
    remove_label_placeholders,
)

from hhh4b2tau.util import round_sig

hist = maybe_import("hist")
np = maybe_import("numpy")
mpl = maybe_import("matplotlib")
plt = maybe_import("matplotlib.pyplot")
mplhep = maybe_import("mplhep")
od = maybe_import("order")

logger = law.logger.get_logger(__name__)


def separate_sig_bkg_hists(hists: OrderedDict):
    # separate histograms into signal and background based on process tag 'is_signal'
    h_sig, h_bkg = None, None
    for proc_inst, h in hists.items():
        # if proc_inst.has_tag("is_signal"):
        if proc_inst.has_tag("signal"):
            h_sig = h + h_sig if h_sig else h
        else:
            h_bkg = h + h_bkg if h_bkg else h

    if not h_sig:
        raise Exception(
            "No signal processes given. Remember to add the 'signal' tag to your "
            "signal processes",
        )
    if not h_bkg:
        raise Exception("No background processes given")

    return h_sig, h_bkg

# taken from hh -> bbww analysis
def plot_roc(
    hists: OrderedDict,
    config_inst: od.Config,
    category_inst: od.Category,
    variable_insts: list[od.Variable],
    style_config: dict | None = None,
    density: bool | None = False,
    shape_norm: bool = False,
    yscale: str | None = "",
    hide_errors: bool | None = None,
    cumsum: bool = False,
    reversed_cumsum: bool = False,
    bkg_rej: bool = False,
    process_settings: dict | None = None,
    variable_settings: dict | None = None,
    cms_label: str = "simpw",
    whitespace_fraction: float = 0.3,
    magnitudes: float = 4,
    **kwargs,
) -> plt.Figure:
    """
    Plotting function to create a single line presenting signal vs background ratio.
    Some plot parameters might be supported but have not been tested.
    Exemplary task call:

    .. code-block:: bash
        law run cf.PlotVariables1D --version prod1 \
            --processes hh_ggf_hbb_hvvqqlnu_kl1_kt1,tt_sl --variables jet1_pt \
            --plot-function hhh4b2tau.plotting.s_over_b.plot_ROC \
            --general-settings reversesed_cumsum,shape_norm
    """
    if  (not cumsum and not reversed_cumsum) or (cumsum and reversed_cumsum):
        raise Exception("Exactly one cumulative cumsum direction has to be chosen!")

    hists = remove_residual_axis(hists, "shift")
    variable_inst = variable_insts[0]
    hists, variable_style_config = apply_variable_settings(hists, variable_insts, variable_settings)
    hists, process_style_config = apply_process_settings(hists, process_settings)
    hists = apply_density(hists, density)

    # separate histograms into signal and background based on process tag 'is_signal'
    h_sig, h_bkg = separate_sig_bkg_hists(hists)


    # calculate the cumulative sum of the histograms
    for h in (h_sig, h_bkg):
        h_view = h.view()
        h_view.value = np.cumsum(h_view.value) if cumsum else np.cumsum(h_view.value[::-1])[::-1]

        # TODO: implement correct variances; in the meantime, just set them to 0
        h_view.variance = [0] * len(h.view())

        # overwrite view of original histograms
        h[...] = h_view


    if shape_norm:
        h_sig = h_sig / h_sig.values()[0] if reversed_cumsum else h_sig / h_sig.values()[-1]
        h_bkg = h_bkg / h_bkg.values()[0] if reversed_cumsum else h_bkg / h_bkg.values()[-1]


    xaxis = 1 - h_bkg.values() if bkg_rej else h_bkg.values()
    yaxis = h_sig.values() 

    xlabel = "Background Rejection" if bkg_rej else "Background Efficiency"
    ylabel = "Signal Efficiency"


    default_style_config = prepare_style_config(
        config_inst, category_inst, variable_inst, density, shape_norm, yscale,
    )

    # disable autmatic setting of xlim
    default_style_config["ax_cfg"]["xlim"] = None
    default_style_config["ax_cfg"]["xlabel"] = xlabel

    # disable autmatic setting of ylim
    default_style_config["ax_cfg"]["ylim"] = None
    default_style_config["ax_cfg"]["ylabel"] = ylabel

    style_config = law.util.merge_dicts(
        default_style_config,
        process_style_config,
        variable_style_config[variable_inst],
        style_config,
        deep=True,
    )



    # emulate plot_all

    # general mplhep style
    plt.style.use(mplhep.style.CMS)

    # setup figure and axes
    rax = None
    grid_spec = {"left": 0.15, "right": 0.95, "top": 0.95, "bottom": 0.1}
    grid_spec |= style_config.get("gridspec_cfg", {})

    fig, ax = plt.subplots(gridspec_kw=grid_spec)
    axs = (ax,)

    ax.plot(xaxis, yaxis, linestyle='-', color='black')

    # if reversed_cumsum:
    #     ax.xaxis.set_inverted(True)
    #     ax.yaxis.set_inverted(True)
        
    # axis styling
    ax_kwargs = {
        "ylabel": "Counts",
        "xlabel": "Counts",
        "yscale": "linear",
        "xscale": "linear",
    }

    # some default ylim settings based on yscale
    log_y = style_config.get("ax_cfg", {}).get("yscale", "linear") == "log"

    ax_ymin = ax.get_ylim()[1] / 10**magnitudes if log_y else 0.0000001
    ax_ymax = get_position(ax_ymin, ax.get_ylim()[1], factor=1 / (1 - whitespace_fraction), logscale=log_y)
    ax_kwargs.update({"ylim": (ax_ymin, ax_ymax)})

    # prioritize style_config ax settings
    ax_kwargs.update(style_config.get("ax_cfg", {}))

    # some settings cannot be handled by ax.set
    xminorticks = ax_kwargs.pop("xminorticks", ax_kwargs.pop("minorxticks", None))
    yminorticks = ax_kwargs.pop("yminorticks", ax_kwargs.pop("minoryticks", None))
    xloc = ax_kwargs.pop("xloc", None)
    yloc = ax_kwargs.pop("yloc", None)

    del ax_kwargs['xrotation']
    # set all values
    ax.set(**ax_kwargs)

    # set manual configs
    if xminorticks is not None:
        ax.set_xticks(xminorticks, minor=True)
    if yminorticks is not None:
        ax.set_xticks(yminorticks, minor=True)
    if xloc is not None:
        ax.set_xlabel(ax.get_xlabel(), loc=xloc)
    if yloc is not None:
        ax.set_ylabel(ax.get_ylabel(), loc=yloc)


    # label alignment
    fig.align_labels()


    # custom annotation
    log_x = style_config.get("ax_cfg", {}).get("xscale", "linear") == "log"
    annotate_kwargs = {
        "text": "",
        "xy": (
            get_position(*ax.get_xlim(), factor=0.05, logscale=log_x),
            get_position(*ax.get_ylim(), factor=0.95, logscale=log_y),
        ),
        "xycoords": "data",
        "color": "black",
        "fontsize": 22,
        "horizontalalignment": "right" if bkg_rej else "left",
        "verticalalignment": "top",
    }
    annotate_kwargs.update(style_config.get("annotate_cfg", {}))

    if variable_inst.name == 'dhh':
       var_label = r"$D_{HH}$"

    if variable_inst.name == 'chi2':
        var_label = r"$\chi^2$"

    if variable_inst.name == 'chi2_1':
        var_label = r"$\chi^2_{H1}$"

    if variable_inst.name == 'chi2_2':
        var_label = r"$\chi^2_{H2}$"

    # add area under curve
    area = round(abs(np.trapz(yaxis, xaxis, dx=0.001)), 4)

    annotate_kwargs["text"] = "Category: " + annotate_kwargs["text"] + "\n" + "Variable: " + var_label + "\n" + "AUC: " + str(area)

    ax.annotate(**annotate_kwargs)
    # cms label
    if cms_label != "skip":
        cms_label_kwargs = get_cms_label(ax, cms_label)

        cms_label_kwargs.update(style_config.get("cms_label_cfg", {}))
        mplhep.cms.label(**cms_label_kwargs)

    # finalization
    fig.tight_layout()

    return fig, axs


def plot_s_over_b(
    hists: OrderedDict,
    config_inst: od.Config,
    category_inst: od.Category,
    variable_insts: list[od.Variable],
    style_config: dict | None = None,
    density: bool | None = False,
    shape_norm: bool | None = False,
    yscale: str | None = "",
    hide_errors: bool | None = None,
    sqrt_b: bool | None = None,
    cumsum: bool = False,
    reversed_cumsum: bool = False,
    process_settings: dict | None = None,
    variable_settings: dict | None = None,
    **kwargs,
) -> plt.Figure:
    """
    Plotting function to create a single line presenting signal vs background ratio.
    Some plot parameters might be supported but have not been tested.
    Exemplary task call:

    .. code-block:: bash
        law run cf.PlotVariables1D --version prod1 \
            --processes hh_ggf_hbb_hvvqqlnu_kl1_kt1,tt_sl --variables jet1_pt \
            --plot-function hhh4b2tau.plotting.s_over_b.plot_s_over_b \
            --general-settings sqrt_b
    """
    if cumsum and reversed_cumsum:
        raise Exception("We can only do the cumulative cumsum in one direction at the time!")

    hists = remove_residual_axis(hists, "shift")

    variable_inst = variable_insts[0]
    hists, variable_style_config = apply_variable_settings(hists, variable_insts, variable_settings)
    hists, process_style_config = apply_process_settings(hists, process_settings)
    hists = apply_density(hists, density)

    # separate histograms into signal and background based on process tag 'is_signal'
    h_sig, h_bkg = separate_sig_bkg_hists(hists)

    # TODO: include uncertainties
    S_over_B = h_sig[::sum].value / h_bkg[::sum].value
    S_over_sqrtB = h_sig[::sum].value / np.sqrt(h_bkg[::sum].value)
    logger.info(
        f"\n    Integrated S over B:     {round_sig(S_over_B)}" +
        f"\n    Integrated S over sqrtB: {round_sig(S_over_sqrtB)}",
    )

    if cumsum or reversed_cumsum:
        # calculate the cumulative sum of the histograms
        for h in (h_sig, h_bkg):
            h_view = h.view()
            h_view.value = np.cumsum(h_view.value) if cumsum else np.cumsum(h_view.value[::-1])[::-1]

            # TODO: implement correct variances; in the meantime, just set them to 0
            h_view.variance = [0] * len(h.view())

            # overwrite view of original histograms
            h[...] = h_view

    # NOTE: this does not take into account the variances on the background stack
    if sqrt_b:
        h_out = h_sig / np.sqrt(h_bkg.values())
        ylabel = r"$Signal / \sqrt{Background}$"
    else:
        h_out = h_sig / h_bkg.values()
        ylabel = "Signal / Background"

    # draw lines
    plot_config = {}
    line_norm = sum(h_out.values()) if shape_norm else 1
    plot_config["s_over_b"] = {
        "method": "draw_hist",
        "hist": h_out,
        "kwargs": {
            "norm": line_norm,
        },
    }

    default_style_config = prepare_style_config(
        config_inst, category_inst, variable_inst, density, shape_norm, yscale,
    )
    # disable autmatic setting of ylim
    default_style_config["ax_cfg"]["ylim"] = None
    default_style_config["ax_cfg"]["ylabel"] = ylabel

    style_config = law.util.merge_dicts(
        default_style_config,
        process_style_config,
        variable_style_config[variable_inst],
        style_config,
        deep=True,
    )
    if shape_norm:
        style_config["ax_cfg"]["ylabel"] = r"$\Delta N/N$"

    # ratio plot not used here; set `skip_ratio` to True
    kwargs["skip_ratio"] = True

    return plot_all(plot_config, style_config, **kwargs)

def cutflow_s_over_b(
    hists: OrderedDict,
    config_inst: od.Config,
    category_inst: od.Category,
    style_config: dict | None = None,
    density: bool | None = False,
    shape_norm: bool = False,
    yscale: str | None = None,
    hide_errors: bool | None = None,
    sqrt_b: bool | None = None,
    process_settings: dict | None = None,
    **kwargs,
) -> plt.Figure:
    """
    Plotting function to create a single line presenting signal vs background ratio
    after each selection step.
    Also prints a table containing all process yields after each step + S over B ratios
    Some plot parameters might be supported but have not been tested.
    Exemplary task call:

    .. code-block:: bash
        law run cf.PlotCutflow --version prod1 \
            --processes hh_ggf_hbb_hvvqqlnu_kl1_kt1,tt_sl \
            --plot-function hbw.plotting.s_over_b.cutflow_s_over_b \
            --general-settings sqrt_b
    """
    from tabulate import tabulate

    hists = remove_residual_axis(hists, "shift")

    hists = apply_process_settings(hists, process_settings)
    hists = apply_density_to_hists(hists, density)

    selector_steps = list(hists[list(hists.keys())[0]].axes["step"])
    selector_step_labels = config_inst.x("selector_step_labels", {})

    #
    # gather yields per process and step
    #
    yields = defaultdict(list)
    for process_inst, h in hists.items():
        yields["Label"].append(process_inst.name)
        for step in selector_steps:
            step_label = step
            # step_label = selector_step_labels.get(step, step)
            yields[step_label].append(round_sig(h[{"step": step}].value, 4))

    # separate histograms into signal and background based on process tag 'is_signal'
    h_sig, h_bkg = separate_sig_bkg_hists(hists)

    # NOTE: this does not take into account the variances on the background stack
    h_s_over_sqrt_b = h_sig / np.sqrt(h_bkg.values())
    h_s_over_b = h_sig / h_bkg.values()
    if sqrt_b:
        h_out = h_s_over_sqrt_b
        ylabel = r"$Signal / \sqrt{Background}$"
    else:
        h_out = h_s_over_b
        ylabel = "Signal / Background"

    yields["Label"].append("S/sqrt(B)")
    yields["Label"].append("S/B")
    for step in selector_steps:
        step_label = step
        # step_label = selector_step_labels.get(step, step)
        yields[step_label].append(round_sig(h_s_over_sqrt_b[{"step": step}].value, 4))
        yields[step_label].append(round_sig(h_s_over_b[{"step": step}].value, 4))

    # create, print and save the yield table
    yield_table = tabulate(yields, headers="keys", tablefmt="fancy_grid")
    print(yield_table)

    #
    # plotting
    #
    plot_config = {}
    plot_config["s_over_b"] = {
        "method": "draw_hist",
        "hist": h_out,
    }

    # update xticklabels based on config
    xticklabels = []

    for step in selector_steps:
        xticklabels.append(selector_step_labels.get(step, step))

    # setup style config
    if not yscale:
        yscale = "linear"

    default_style_config = {
        "ax_cfg": {
            "ylim": None,
            "ylabel": ylabel,
            "xlabel": "Selection step",
            "xticklabels": xticklabels,
            "yscale": yscale,
        },
        "legend_cfg": {
            "loc": "upper right",
        },
        "annotate_cfg": {"text": category_inst.label},
        "cms_label_cfg": {
            "lumi": config_inst.x.luminosity.get("nominal") / 1000,  # pb -> fb
        },
    }
    style_config = law.util.merge_dicts(default_style_config, style_config, deep=True)

    # ratio plot not used here; set `skip_ratio` to True
    kwargs["skip_ratio"] = True

    fig, (ax,) = plot_all(plot_config, style_config, **kwargs)

    ax.set_xticklabels(xticklabels, rotation=45, ha="right")

    return fig, (ax,)