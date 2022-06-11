


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare.photom import  M_to_lum
import flare.plt as fplt
from flare.LF import models, evo, plots
from flare.LF.literature import UV

fig, ax = fplt.simple()

z = 5.0


bin_w = 0.1
bin_edges = np.arange(27, 30.5, bin_w)
bin_centres = models.bin_centres(bin_edges)


m = getattr(UV, 'FLARES_binned')()

ax.step(m.log10L, np.log10(m.phi[z]), where = 'mid', label = 'original binned', lw=3, c='k', alpha = 0.3) # raw data # CORRECT


# --- rebinned on to finer luminosity

phi = m.phi_binned(z, bin_edges)
ax.step(bin_centres, np.log10(phi)-np.log10(bin_w), where = 'mid', label = 're-binned')

tot_phi = np.sum(phi)

print('total number of galaxies per Mpc:', tot_phi)


from scipy.stats import rv_histogram

hist_dist = rv_histogram((phi, bin_edges))

sample_size = 100000

sample = hist_dist.rvs(size = sample_size)

N, _ = np.histogram(sample, bins = bin_edges)

ax.step(bin_centres, np.log10(N*tot_phi/sample_size)-np.log10(bin_w), where = 'mid', label = 'sampled')


ax.legend()

ax.set_xlabel(r'$\rm log_{10}(L_{FUV}/erg\ s^{-1} Hz^{-1}) $')
ax.set_ylabel(r'$\rm log_{10}(\phi/dex^{-1}\ Mpc^{-1}) $')

fig.savefig('figs/sampled.pdf')
