"""
Plots a tree object as a dendrogram for visualisation with plotly.
"""

from itertools import repeat
from copy import deepcopy
import re

from plotly.graph_objs import graph_objs

import numpy as np
from scipy.cluster import hierarchy as sph


def dendrogram(
    tree,
    orientation='bottom',
    colorscale=None,
    width=np.inf,
    height=np.inf,
    xaxis='xaxis',
    yaxis='yaxis',
    title=None,
    xlabel=None,
    ylabel=None,
    hovertext=None,
    margin_scalar=12,
):
    layout = {xaxis: {}, yaxis: {}}
    if title is not None:
        layout["title"] = title

    if xlabel is not None:
        layout[xaxis]["title"] = xlabel

    if ylabel is not None:
        layout[yaxis]["title"] = ylabel

    sign = _get_sign(orientation, xaxis, yaxis)

    (dd_traces, xvals, yvals, ordered_labels, leaves) = _get_traces(
        hierarchy=tree.tree,
        labels=list(tree.onehot_df.index),
        threshold=tree.cutoff,
        orientation=orientation,
        sign=sign,
        xaxis=xaxis,
        yaxis=yaxis,
        hovertext=hovertext,
        colorscale=colorscale,
    )

    yvals_flat = yvals.flatten()
    xvals_flat = xvals.flatten()

    zero_vals = {x for x, y in zip(xvals_flat, yvals_flat) if y == 0.0}

    if len(zero_vals) > len(yvals) + 1:
        # If the length of zero_vals is larger than the length of yvals,
        # it means that there are wrong vals because of the identicial samples.
        # Three and more identicial samples will make the yvals of spliting
        # center into 0 and it will accidentally take it as leaves.

        l_border = int(min(zero_vals))
        r_border = int(max(zero_vals))

        correct_leaves_pos = range(
            l_border,
            r_border + 1,
            int((r_border - l_border) / len(yvals))
        )

        # Regenerating the leaves pos from the zero_vals with equal intervals.
        zero_vals = list(correct_leaves_pos)
    else:
        zero_vals = list(zero_vals)

    zero_vals.sort()
    layout = _figure_layout(
        layout,
        width,
        height,
        ordered_labels,
        xaxis,
        yaxis,
        sign,
        orientation,
        zero_vals,
    )

    # Adjust the margins to account for long sample labels.
    longest_label = max([len(l) for l in ordered_labels])
    if "margin" not in layout:
        ori = orientation[0]  # This spec used just l, r, t, or b
        layout["margin"] = {ori: int(margin_scalar) * (longest_label + 1)}

    if hovertext is None:
        cluster_map = dict(zip(tree.onehot_df.index, tree.clusters))
        data = _format_cluster_hovertexts(dd_traces, layout, cluster_map)
    else:
        data = dd_traces

    return graph_objs.Figure(data, layout)


def _format_cluster_hovertexts(data, layout, cluster_map):
    """ Updates data elements to add hover text for cluster names. """

    labels = dict(zip(layout["xaxis"]["tickvals"],
                      layout["xaxis"]["ticktext"]))
    hover_labs = {
        k: "cluster: {}, component: {}".format(cluster_map[v], v)
        for k, v
        in labels.items()
    }

    out = []
    for scatter in data:
        # Avoid side-effects just in case.
        scatter = deepcopy(scatter)
        xaxis = scatter["xaxis"]
        yaxis = scatter["yaxis"]
        textl = None
        textr = None

        # If left is at 0.
        if scatter[yaxis][0] == 0.0:
            textl = hover_labs.get(scatter[xaxis][0], None)

        if scatter[yaxis][3] == 0.0:
            textr = hover_labs.get(scatter[xaxis][3], None)

        scatter["text"] = [textl, textl, textr, textr]
        out.append(scatter)

    return out


def _get_sign(orientation, xaxis, yaxis):
    """ Helper to find multiplicative factors to get rotations right. """

    sign = {xaxis: 1, yaxis: 1}

    if orientation in ['left', 'bottom']:
        sign[xaxis] = 1
    else:
        sign[xaxis] = -1

    if orientation in ['right', 'bottom']:
        sign[yaxis] = 1
    else:
        sign[yaxis] = -1
    return sign


def _get_color_dict(colorscale=None):
    """ Returns colorscale used for dendrogram tree clusters.

    Keyword arguments:
    colorscale -- Colors to use for the plot in rgb format.
        Should have 8 colours.
    """

    # These are the color codes returned for dendrograms
    # We're replacing them with nicer colors
    default_colors = [
        ('b', "blue"),
        ('c', "cyan"),
        ('g', "green"),
        ('k', "black"),
        ('m', "magenta"),
        ('r', "red"),
        ('w', "white"),
        ('y', "yellow"),
    ]

    if colorscale is None:
        colorscale = [
            'rgb(0,116,217)',  # instead of blue
            'rgb(35,205,205)',  # cyan
            'rgb(61,153,112)',  # green
            'rgb(40,35,35)',  # black
            'rgb(133,20,75)',  # magenta
            'rgb(255,65,54)',  # red
            'rgb(255,255,255)',  # white
            'rgb(255,220,0)',  # yellow
        ]

    mapped_colors = {
        k: o
        for (k, v), o
        in zip(default_colors, colorscale)
    }

    return mapped_colors


def _get_traces(
    hierarchy,
    labels,
    threshold,
    orientation,
    sign,
    xaxis,
    yaxis,
    hovertext=None,
    colorscale=None,
):
    """ Format the dendrogram nodes/clades as edges in graph. """

    # Scipy does most of the heavy lifting.
    dendro = sph.dendrogram(
        hierarchy,
        orientation=orientation,
        labels=labels,
        no_plot=True,
        color_threshold=threshold
    )

    icoords = np.array(dendro["icoord"])
    dcoords = np.array(dendro["dcoord"])

    # xs and ys are arrays of 4 points that make up the '∩' shapes
    # of the dendrogram tree
    if orientation in ["top", "bottom"]:
        xs = icoords
        ys = dcoords
    else:
        ys = icoords
        xs = dcoords

    # Multiply by 1 or -1 to get correct rotation
    xs = sign[xaxis] * xs
    ys = sign[yaxis] * ys

    ordered_labels = dendro["ivl"]

    if hovertext is None:
        # Infinite generator of None values.
        hovertext = repeat(None)

    color_map = _get_color_dict(colorscale)
    colors = [color_map[k] for k in dendro["color_list"]]
    traces = _trace_as_scatter(xs, ys, colors, hovertext, xaxis, yaxis)

    return traces, icoords, dcoords, ordered_labels, dendro["leaves"]


def _trace_as_scatter(xs, ys, colors, hovertext, xaxis, yaxis):
    """ Formats the values from the scipy dendro as a list of plotly
    compatible dicts. There will be converted into Scatter objects by Figure."
    """

    traces = []
    for x, y, color, htext in zip(xs, ys, colors, hovertext):
        trace = {
            "type": "scatter",
            "x": x,
            "y": y,
            "mode": 'lines',
            "marker": {"color": color},
            "text": htext,
            "hoverinfo": 'text',
            "xaxis": _axis_index(xaxis, "x"),
            "yaxis": _axis_index(yaxis, "y"),
        }
        # The axis index stuff seems to be an anticipation of people having
        # multiple plots in the same figure.

        traces.append(trace)
    return traces


def _axis_index(string, prefix):
    """ In the case of a plot with sub-plots, the axes may have numbers after
    This captures those number for us.

    Note, I can't see why this is actually necessary, (i.e. if multiplots
    will actually break it) but I'm porting it over from the original plotly
    code just in case.
    """

    regex = re.compile(r"(?P<integer>\d+)$")

    match = regex.search(string)
    if match is None:
        return prefix
    else:
        return "".join([prefix, match.groupdict().get("integer", "")])


def _axis_layout_defaults(layout={}):
    """ Updates layout dict with some defauls.

    Returns copy of data, i.e. doesn't mutate.
    """

    axis_defaults = {
        'type': 'linear',
        'ticks': 'outside',
        'mirror': 'allticks',
        'rangemode': 'tozero',
        'showticklabels': True,
        'zeroline': False,
        'showgrid': False,
        'showline': True,
    }

    layout = deepcopy(layout)
    layout.update(axis_defaults)
    return layout


def _key_labels_layout(sign, zero_vals, labels):
    """ Format the leaf labels to be shown on the axis. """

    layout = {
        "tickvals": [sign * zv for zv in zero_vals],
        "ticktext": labels,
        "tickmode": "array",
    }
    return layout


def _figure_layout_defaults(width, height, layout={}):
    """ Updates layout dict with some defaults.

    Returns copy of data, i.e. no mutation.
    """

    figure_defaults = {
        'showlegend': False,
        'autosize': False,
        'hovermode': 'closest',
        'width': width,
        'height': height
    }

    layout = deepcopy(layout)
    layout.update(figure_defaults)
    return layout


def _figure_layout(
    layout,
    width,
    height,
    labels,
    xaxis,
    yaxis,
    sign,
    orientation,
    zero_vals
):
    """ Sets up figure and axis layouts, and adds leaf labels. """

    if len(labels) > 0:
        if orientation in ["left", "right"]:
            axis_key_labels = yaxis
        else:
            axis_key_labels = xaxis

        label_layout = _key_labels_layout(
            sign=sign[axis_key_labels],
            zero_vals=zero_vals,
            labels=labels
        )

        layout[axis_key_labels].update(label_layout)

    figure_defaults = _figure_layout_defaults(width, height)
    layout.update(figure_defaults)

    layout[xaxis].update(_axis_layout_defaults())
    layout[yaxis].update(_axis_layout_defaults())
    return layout
