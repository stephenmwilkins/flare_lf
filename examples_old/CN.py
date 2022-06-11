


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import cmasher as cmr

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import flare

import flare.plt as fplt

from flare.LF import completeness, evo, literature, plots
from flare.photom import m_to_flux

cosmo = flare.default_cosmo()

fig, ax = fplt.simple()


models = ['FLARES', 'FLARES_binned']


m_limits = [27., 28., 29., 30.]
colors = cmr.take_cmap_colors('cmr.neon', len(m_limits))


for model, ls in zip(models, ['-','-.','--',':']):

    m = getattr(literature, model)()

    bin_edges, bin_centres, N_ = m.N(redshift_limits = [5., 15.], log10L_limits = [27., 30.], dz = 0.1, dlog10L = 0.01)

    for m_limit, color in zip(m_limits, colors):

        flux_limit = flare.photom.m_to_flux(m_limit)

        # --- simple completeness
        c = completeness.completeness_cut(bin_centres, flux_limit, cosmo = flare.default_cosmo())

        # --- get expected number of galaxies
        N = np.multiply(N_, c)

        n = np.sum(N, axis=0)[::-1]
        cn = np.cumsum(n)

        ax.plot(bin_centres['z'][::-1], np.log10(cn), c=color, alpha = 0.7, lw=1., ls = ls, label = rf'$\rm m<{m_limit}$')

ax.legend()
# ax.set_ylimits([-10, 2])


ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm log_{10}(N(>z, z<15)/arcmin^2)$')

fig.savefig('figs/CN.pdf')
