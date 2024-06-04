'''
Custom plot improvements using matplotlib
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import griddata

from .colormap import cmap_white_blue_red


def uber_plot(x, y, label=None, figname=None, clf=True, figsize=(6.4, 4.8),
              sorting=True, xerr=None, yerr=None,
              linewidth=2, linestyle='-', marker=None, markersize=8,
              color=None, alpha=None,
              xlim=None, ylim=None,
              xlabel=None, ylabel=None, legend=None, title=None,
              grid=False, fontsize=14, savefig=None, dpi=100,
              scilimits=(-2, 4)):

    if figsize == 'golden':
        figsize = (2 * (1 + np.sqrt(5)) * (4.8 / 4.0), 4.8)

    plt.figure(figname, figsize)

    if clf:
        plt.clf()

    fig = plt.gcf()
    axes = plt.gca()

    plot_opt = {
        'label': label,
        'linewidth': linewidth,
        'linestyle': linestyle,
        'marker': marker,
        'markersize': markersize,
        'color': color,
        'alpha': alpha}

    if xerr is None and yerr is None:
        plot_func = plt.plot
    else:
        plot_func = plt.errorbar
        if xerr is not None:
            plot_opt['xerr'] = xerr
        if yerr is not None:
            plot_opt['yerr'] = yerr

    if sorting:
        array_order = np.argsort(x)
        if 'xerr' in plot_opt:
            if isinstance(plot_opt['xerr'], list):
                plot_opt['xerr'] = [np.array(plot_opt['xerr'][0][array_order]),
                                    np.array(plot_opt['xerr'][1][array_order])]
            else:
                plot_opt['xerr'] = np.array(plot_opt['xerr'][array_order])
        if 'yerr' in plot_opt:
            if isinstance(plot_opt['yerr'], list):
                plot_opt['yerr'] = [np.array(plot_opt['yerr'][0][array_order]),
                                    np.array(plot_opt['yerr'][1][array_order])]
            else:
                plot_opt['yerr'] = np.array(plot_opt['yerr'][array_order])
        plot_func(x[array_order], y[array_order], **plot_opt)
    else:
        plot_func(x, y, **plot_opt)

    plot_finish(fig=fig, axes=axes,
                xlim=xlim, ylim=ylim,
                xlabel=xlabel, ylabel=ylabel, legend=legend,
                grid=grid, fontsize=fontsize, savefig=savefig, dpi=dpi,
                scilimits=scilimits, title=title)

    return fig, axes


def uber_scatter(x, y, c=None, s=None, label=None, figname=None, clf=True,
                 figsize=(6.4, 4.8), sorting=True,
                 linewidth=2, linestyle='-', marker=None,
                 alpha=None, cmap='jet', cbar=False,
                 xlim=None, ylim=None, clim=None, clabel=None,
                 xlabel=None, ylabel=None, legend=None, title=None,
                 grid=False, fontsize=14, savefig=None, dpi=100,
                 scilimits=(-2, 4)):

    if figsize == 'golden':
        figsize = (2 * (1 + np.sqrt(5)) * (4.8 / 4.0), 4.8)

    plt.figure(figname, figsize)

    if clf:
        plt.clf()

    fig = plt.gcf()
    axes = plt.gca()

    plot_opt = {
        'cmap': cmap,
        'label': label,
        'linewidth': linewidth,
        'linestyle': linestyle,
        'marker': marker,
        'alpha': alpha,
        'zorder': 2}

    if c is not None:
        plot_opt['c'] = c
    if s is not None:
        plot_opt['s'] = s

    if sorting:
        array_order = np.argsort(x)
        if 'c' in plot_opt:
            plot_opt['c'] = np.array(plot_opt['c'][array_order])
        if 's' in plot_opt:
            plot_opt['s'] = np.array(plot_opt['s'][array_order])
        plt.scatter(x[array_order], y[array_order], **plot_opt)
    else:
        plt.scatter(x, y, **plot_opt)

    if cbar:
        cbar = plt.colorbar()
    else:
        cbar = None

    plot_finish(fig=fig, axes=axes,
                xlim=xlim, ylim=ylim, clim=clim, cbar=cbar, clabel=clabel,
                xlabel=xlabel, ylabel=ylabel, legend=legend,
                grid=grid, fontsize=fontsize, savefig=savefig, dpi=dpi,
                scilimits=scilimits, title=title)

    return fig, axes


def uber_bar(x, y, label=None, figname=None, clf=True, figsize=(6.4, 4.8),
             linewidth=1, width=None, edgecolor='black',
             color=None, alpha=None,
             xlim=None, ylim=None,
             xlabel=None, ylabel=None, legend=None, title=None,
             grid=False, fontsize=14, savefig=None, dpi=100,
             scilimits=(-2, 4)):

    if figsize == 'golden':
        figsize = (2 * (1 + np.sqrt(5)) * (4.8 / 4.0), 4.8)

    plt.figure(figname, figsize)

    if clf:
        plt.clf()

    fig = plt.gcf()
    axes = plt.gca()

    if width is None:
        width = (x[1] - x[0])

    plt.bar(x, y, label=label, width=width, color=color, alpha=alpha,
            linewidth=linewidth, edgecolor=edgecolor)

    plot_finish(fig=fig, axes=axes,
                xlim=xlim, ylim=ylim,
                xlabel=xlabel, ylabel=ylabel, legend=legend,
                grid=grid, fontsize=fontsize, savefig=savefig, dpi=dpi,
                scilimits=scilimits, title=title)

    return fig, axes


def uber_imshow(img, x=None, y=None, figname=None,
                clf=True, figsize=(6.4, 4.8),
                img_interpolation='bilinear', data_interpolation=None,
                cmap=cmap_white_blue_red,
                cbar=False, clim=None, clabel=None,
                origin='lower',
                xlim=None, ylim=None, contour=None,
                xlabel=None, ylabel=None, legend=None, title=None,
                grid=False, fontsize=14, savefig=None, dpi=100,
                scilimits=(-2, 4)):

    if figsize == 'golden':
        figsize = (2 * (1 + np.sqrt(5)) * (4.8 / 4.0), 4.8)

    plt.figure(figname, figsize)

    if clf:
        plt.clf()

    fig = plt.gcf()
    axes = plt.gca()

    if x is not None and y is not None:
        dx = np.min(np.diff(x))
        dy = np.min(np.diff(y))

        # Checking if need to reinterpolate
        dx_allclose = np.all(np.isclose(np.unique(np.diff(x)), dx))
        dy_allclose = np.all(np.isclose(np.unique(np.diff(y)), dy))

        if (not dx_allclose or not dy_allclose) and \
                data_interpolation is not None:

            # Reinterpolation of data
            if data_interpolation == 'nearest':
                method = data_interpolation
            else:
                method = data_interpolation[2:]

            x_grid, y_grid = np.meshgrid(x, y)

            x_min = np.min(x)
            x_max = np.max(x)
            y_min = np.min(y)
            y_max = np.max(y)

            x_interp = np.arange(np.min(x), np.max(x)+0.01*dx, dx)
            y_interp = np.arange(np.min(y), np.max(y)+0.01*dy, dy)

            x_grid_interp, y_grid_interp = np.meshgrid(
                x_interp, y_interp)

            img = griddata(
                ((x_grid.flatten()-x_min)/(x_max-x_min),
                 (y_grid.flatten()-y_min)/(y_max-y_min)),
                img.T.flatten(),
                ((x_grid_interp-x_min)/(x_max-x_min),
                 (y_grid_interp-y_min)/(y_max-y_min)),
                method=method).T

            x_grid, y_grid = (x_grid_interp, y_grid_interp)

        if origin == 'lower':
            extent = [x[0]-dx/2, x[-1]+dx/2, y[0]-dy/2, y[-1]+dy/2]
        elif origin == 'upper':
            extent = [x[-1]+dx/2, x[0]-dx/2, y[-1]+dy/2, y[0]-dy/2]
            
    else:
        x_grid, y_grid = np.meshgrid(
            np.arange(img.shape[0]),
            np.arange(img.shape[1]))
        extent = None

    plot_opt = {
        'aspect': 'auto',
        'origin': origin,
        'interpolation': img_interpolation,
        'cmap': cmap,
        'extent': extent}

    if clim is not None:
        plot_opt['vmin'] = clim[0]
        plot_opt['vmax'] = clim[1]

    plt.imshow(img.T, **plot_opt)

    if cbar:
        cbar = plt.colorbar()
    else:
        cbar = None

    if contour is not None:
        if isinstance(contour, list) or isinstance(contour, np.ndarray):
            levels = contour
        else:
            levels = None
        CS = plt.contour(x_grid, y_grid, img.T,
                         colors='k', levels=levels)
        plt.clabel(CS, fontsize=9, inline=1)

    plot_finish(fig=fig, axes=axes,
                xlim=xlim, ylim=ylim, clim=clim, cbar=cbar, clabel=clabel,
                xlabel=xlabel, ylabel=ylabel, legend=legend,
                grid=grid, fontsize=fontsize, savefig=savefig, dpi=dpi,
                scilimits=scilimits, title=title)

    return fig, axes, cbar


def plot_finish(fig=None, axes=None,
                xlim=None, ylim=None, clim=None, cbar=None, clabel=None,
                xlabel=None, ylabel=None, legend=None, title=None,
                grid=False, fontsize=14, savefig=None, dpi=100,
                scilimits=(-2, 4),
                n_minor=5):

    if fig is None:
        fig = plt.gcf()
    if axes is None:
        axes = plt.gca()

    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)

    if cbar is not None and clabel is not None:
        cbar.set_label(clabel, fontsize=fontsize)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if clim is not None:
        plt.clim(clim)

    if grid:
        plt.grid(grid)
        plt.grid(grid, 'minor', alpha=0.3)

    if title is not None:
        plt.title(title)

    axes.spines['top'].set_linewidth(1.5)
    axes.spines['right'].set_linewidth(1.5)
    axes.spines['bottom'].set_linewidth(1.5)
    axes.spines['left'].set_linewidth(1.5)

    axes.tick_params(axis='x', labelsize=fontsize)
    axes.tick_params(axis='y', labelsize=fontsize)

    axes.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    axes.yaxis.set_minor_locator(AutoMinorLocator(n=5))

    axes.tick_params(which='major', axis='x', color='k', length=6, width=1.5)
    axes.tick_params(which='major', axis='y', color='k', length=6, width=1.5)

    axes.xaxis.major.formatter._useMathText = True
    axes.yaxis.major.formatter._useMathText = True
    axes.ticklabel_format(style='sci', axis='both', scilimits=scilimits)

    if cbar is not None:

        cbar.ax.spines['top'].set_linewidth(1.5)
        cbar.ax.spines['right'].set_linewidth(1.5)
        cbar.ax.spines['bottom'].set_linewidth(1.5)
        cbar.ax.spines['left'].set_linewidth(1.5)

        cbar.ax.tick_params(axis='x', labelsize=fontsize)
        cbar.ax.tick_params(axis='y', labelsize=fontsize)

        cbar.ax.xaxis.set_minor_locator(AutoMinorLocator(n=n_minor))
        cbar.ax.yaxis.set_minor_locator(AutoMinorLocator(n=n_minor))

        cbar.ax.tick_params(which='major', axis='x', color='k', length=6,
                            width=1.5)
        cbar.ax.tick_params(which='major', axis='y', color='k', length=6,
                            width=1.5)

        try:
            cbar.ax.xaxis.major.formatter._useMathText = True
            cbar.ax.yaxis.major.formatter._useMathText = True
            cbar.ax.ticklabel_format(style='sci', axis='both', scilimits=scilimits)
        except:
            pass

    plt.tight_layout()

    if legend is not None:
        plt.legend(loc=legend, prop={'size': fontsize})

    if savefig is not None:
        fig.savefig(savefig, dpi=dpi)
