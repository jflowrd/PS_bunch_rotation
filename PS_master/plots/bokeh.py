'''
Custom plot improvements using bokeh
'''

import numpy as np
import bokeh.plotting as blt
from bokeh.models import Range1d
from bokeh.palettes import Category10

list_colors = Category10[10]
global bokeh_uber_plot_index_color
bokeh_uber_plot_index_color = 0


def uber_plot(x, y, label=None, fig=None, figsize=(6.4, 4.8),
              sorting=True, xerr=None, yerr=None,
              linewidth=2, linestyle='-', marker=None, markersize=8,
              color=None, alpha=1,
              xlim=None, ylim=None,
              xlabel=None, ylabel=None, legend=None, title=None,
              sizing_mode='stretch_both', grid=True, fontsize=14,
              show=True, savefig=None):

    if figsize == 'golden':
        figsize = (2 * (1 + np.sqrt(5)) * (4.8 / 4.0), 4.8)

    plot_opt = {}

    if fig is None:
        fig = blt.figure(plot_width=int(figsize[0]*100),
                         plot_height=int(figsize[1]*100),
                         title=title,
                         x_axis_label=xlabel,
                         y_axis_label=ylabel,
                         sizing_mode=sizing_mode)

    if xerr is None and yerr is None:
        plot_func = fig.line
    else:
        pass
#         plot_func = plt.errorbar
#         if xerr is not None:
#             plot_opt['xerr'] = xerr
#         if yerr is not None:
#             plot_opt['yerr'] = yerr

    if title is not None:
        fig.title.text_font_size = '%dpt' % (fontsize)
    fig.xaxis.axis_label_text_font_size = '%dpt' % (fontsize)
    fig.xaxis.major_label_text_font_size = '%dpt' % (fontsize)
    fig.yaxis.axis_label_text_font_size = '%dpt' % (fontsize)
    fig.yaxis.major_label_text_font_size = '%dpt' % (fontsize)

    if color is None:
        global bokeh_uber_plot_index_color
        color = list_colors[bokeh_uber_plot_index_color]
        bokeh_uber_plot_index_color += 1

    if sorting:
        array_order = np.argsort(x)
#        if 'xerr' in plot_opt:
#            if isinstance(plot_opt['xerr'], list):
#                plot_opt['xerr'] = [np.array(plot_opt['xerr'][0][array_order]),
#                                    np.array(plot_opt['xerr'][1][array_order])]
#            else:
#                plot_opt['xerr'] = np.array(plot_opt['xerr'][array_order])
#        if 'yerr' in plot_opt:
#            if isinstance(plot_opt['yerr'], list):
#                plot_opt['yerr'] = [np.array(plot_opt['yerr'][0][array_order]),
#                                    np.array(plot_opt['yerr'][1][array_order])]
#            else:
#                plot_opt['yerr'] = np.array(plot_opt['yerr'][array_order])
        plot_func(x[array_order], y[array_order],
                  line_width=linewidth,
                  color=color,
                  legend=legend,
                  alpha=alpha)
    else:
        plot_func(x, y,
                  line_width=linewidth,
                  color=color,
                  legend=legend,
                  alpha=alpha)

    if xlim is not None:
        fig.x_range = Range1d(xlim[0], xlim[1])
    if ylim is not None:
        fig.y_range = Range1d(ylim[0], ylim[1])

    if not grid:
        fig.xgrid.visible = False
        fig.ygrid.visible = False

    if savefig is not None:
        blt.output_file(savefig)

    if show:
        blt.show(fig)

    return fig
