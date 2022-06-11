


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare.photom import  M_to_lum
import flare.plt as fplt


from flare_lf import models, evo, plots
from flare_lf.uv import model


fig, ax = fplt.simple()

z = 5.0


bin_w = 0.1
bin_edges = np.arange(27, 30.5, bin_w)
bin_centres = models.bin_centres(bin_edges)


print(model.FLARES)

m = model.FLARES.binned()



ax.step(m.log10L, np.log10(m.phi[z]), where = 'mid', label = 'original binned', lw=3, c='k', alpha = 0.3) # raw data # CORRECT


# --- rebinned on to finer luminosity

phi = m.phi_binned(z, bin_edges)
ax.step(bin_centres, np.log10(phi)-np.log10(bin_w), where = 'mid', label = 're-binned')

# p = m.p(z)
# ax.step(p['log10L'], np.log10(p['phi']), where = 'mid', label = 're-binned')

# Schechter function

# m = getattr(literature, 'FLARES')()
#
# phi = m.phi_binned(z, bin_edges)
# ax.step(bin_centres, np.log10(phi)-np.log10(bin_w), where = 'mid', label = 'Schechter')
#
# ax.axvline(np.log10(M_to_lum(m.p(z)['M*'])), c='k', alpha=0.2, ls='--')
# ax.axhline(m.p(z)['log10phi*'], c='k', alpha=0.2, ls='--')

# m = getattr(literature, 'FLARES_DPL')()
#
# phi = m.phi_binned(z, bin_edges)
#
# # ax.plot(bin_centres, np.log10(phi))
# ax.step(bin_centres, np.log10(phi), where = 'mid', label = 'DPL')



ax.legend()




ax.set_xlabel(r'$\rm log_{10}(L_{FUV}/erg\ s^{-1} Hz^{-1}) $')
ax.set_ylabel(r'$\rm log_{10}(\phi/dex^{-1}\ Mpc^{-1}) $')

fig.savefig('figs/basic_lf.pdf')
