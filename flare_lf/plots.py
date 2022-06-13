
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl

import cmasher as cmr # provides wider range of cmaps, see https://cmasher.readthedocs.io

from . import evo


def simple_fig(fig_size = 3.5):

    if type(fig_size) == float or type(fig_size) == int:
        fig_size = (fig_size, fig_size)

    fig = plt.figure(figsize = fig_size)

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax




def logerr_lo(phi, sigma, ulims):
    err = np.log10(phi) - np.log10(phi - sigma)
    for i, x in enumerate(err):
        if not np.isfinite(x):
            err[i] = 100.
        if ulims[i] == True:
            err[i] = 0.5
    return err


def logerr_hi(phi, sigma):
    err = np.log10(phi + sigma) - np.log10(phi)
    for i, x in enumerate(err):
        if not np.isfinite(x):
            err[i] = 100.
    return err

def _integ(x, a):
    return x ** (a) * np.exp(-x)


def _integ2(x, a):
    return 0.4*np.log(10)*10**(-0.4*x* (a+1)) * np.exp(-10**(-0.4*x))


def _integ_dblpow(x, a, b):
    return 1 / (10 ** (x*(a+1)) + 10 ** (x*(b+1)))





def plot_lf(z, models, lum_type = 'Lnu', cmap = 'cmr.neon', save = False, show = True):

    """ plot the published luminosity function (i.e. no interpolation) """

    fig, ax = simple_fig()

    colors = cmr.take_cmap_colors(cmap, len(models))

    for model, color in zip(models, colors):

        m = evo.read(model)

        if lum_type == 'Lnu':
            ax.step(m.log10L[z], m.log10phi[z], where = 'mid', color = color, label = rf'$\rm {m.name} $')
        if lum_type == 'M':
            ax.step(m.M[z], m.log10phi_mag[z], where = 'mid', color = color, label = rf'$\rm {m.name} $')


    ax.legend(fontsize = 8)

    if lum_type == 'Lnu':
        ax.set_xlabel(r'$\rm \log_{10}(L/erg\ s^{-1}\ Hz^{-1}) $')
        ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ dex^{-1}) $')

    if lum_type == 'M':
        ax.set_xlabel(r'$\rm \log_{10}(M) $')
        ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ mag^{-1}) $')

    if show: plt.show()
    if save: fig.savefig(f'figs/{save}.pdf')

    return fig, ax






def plot_lf_evo(model, lum_type = 'Lnu', cmap = 'cmr.neon', save = False, show = True):

    """ plot the published luminosity function (i.e. no interpolation) """

    fig, ax = simple_fig()

    if model in evo.list_models_binned():

        m = evo.read(f'models/binned/{model}')

        print(m.redshifts)

        colors = cmr.take_cmap_colors(cmap, len(m.redshifts))

        for z, color in zip(m.redshifts, colors):

            if lum_type == 'Lnu':
                ax.step(m.log10L[z], m.log10phi[z], where = 'mid', color = color, label = rf'$\rm z={z:.0f} $')
            if lum_type == 'M':
                ax.step(m.M[z], m.log10phi_mag[z], where = 'mid', color = color, label = rf'$\rm z={z:.0f} $')

    # if model in evo.list_models_binned():
    #
    #     model = 'flares_binned'
    #
    #     m = evo.read(model, scheme = 'interp')
    #
    #     print(m.p(5.0))

    ax.legend(fontsize = 8)

    if lum_type == 'Lnu':
        ax.set_xlabel(r'$\rm \log_{10}(L/erg\ s^{-1}\ Hz^{-1}) $')
        ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ dex^{-1}) $')

    if lum_type == 'M':
        ax.set_xlabel(r'$\rm \log_{10}(M) $')
        ax.set_ylabel(r'$\rm \log_{10}(\phi/Mpc^{-3}\ mag^{-1}) $')

    if show: plt.show()
    if save: fig.savefig(f'figs/{save}.pdf')

    return fig, ax









def evo_plot(bin_edges, N, cosmo=False, f_limits=False, save_file=False, vmin = -8., vmax = 0.0):
    # --- make nice plot

    if not cosmo: cosmo = flare.core.default_cosmo()

    fig = plt.figure(figsize=(6, 5))

    X, Y = np.meshgrid(bin_edges['z'], bin_edges['log10L'])

    cm = plt.get_cmap('plasma')
    plt.pcolormesh(X, Y, np.log10(N), cmap=cm, vmin = -8., vmax = 0.0)

    # --- draw lines of constant flux

    if type(f_limits) is list or type(f_limits) is np.ndarray or type(f_limits) is range:

        for f_limit in f_limits:
            plt.plot(bin_edges['z'], np.log10(flux_to_L(f_limit, cosmo, bin_edges['z'])), 'k--', alpha=0.8)

    if type(f_limits) is float:
        plt.plot(bin_edges['z'], np.log10(flux_to_L(f_limits, cosmo, bin_edges['z'])), 'k--', alpha=0.8)

    bar = plt.colorbar(orientation='vertical')
    bar.set_label(r'$\rm log_{10}(N \; / \; arcmin^{-2})$', rotation=90)

    plt.ylabel(r"$\rm log_{10}(L_{\nu} \; / \; erg\, s^{-1}\, Hz^{-1})$")
    plt.xlabel(r"$\rm z$")
    plt.xlim(min(bin_edges['z']), max(bin_edges['z']))
    plt.ylim(min(bin_edges['log10L']), max(bin_edges['log10L']))

    if save_file == False:
        return fig
    else:
        plt.savefig(save_file + '.png', dpi=300)
