


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare.photom import  M_to_lum, lum_to_M
import flare.plt as fplt
from flare.LF import models, evo, plots
from flare.LF.literature import UV, SFR


Kappa_UV = 1.15E-28
# Kappa_UV = 1.4E-28
L_limit = 1E28
M_limit = lum_to_M(L_limit)
print(M_limit)



fig, ax, ax2 = fplt.twinx() # UV luminosity density



# --- UV LF models

for model in ['FLARES_binned', 'FLARES_binned_intrinsic', 'Bluetides', 'Ma2019']:

    m = getattr(UV, model)()

    rhos = []

    print(model, '-'*5)

    for z in m.redshifts:
        rho = m.density(z, L_limit)
        print(z, rho, np.log10(rho))
        rhos.append(np.log10(rho))

    ax.plot(m.redshifts, rhos, label = rf'$\rm {m.name}$', ls = '-', lw=1)



# --- observations



for model in ['Bouwens2021','Mcleod16']:

    m = getattr(UV, model)()

    rhos = []

    for z in m.redshifts:

        rho = m.density(z, L_limit)
        rhos.append(rho)

        # --- sample the observed uncertainties

    ax.scatter(m.redshifts, np.log10(np.array(rhos)), label = rf'$\rm {m.name}$', s = 10)




# --- SFR DF models

for model in ['FLARES_binned']:

    m = getattr(SFR, model)()

    sfrds = []

    print(model, '-'*5)

    for z in m.redshifts:

        sfrd = m.density(z, 1.0)
        print(z, sfrd, np.log10(sfrd))
        sfrds.append(np.log10(sfrd))

    ax2.plot(m.redshifts, sfrds, label = rf'$\rm {m.name}$', ls = '--', lw=1)



ax.legend(fontsize = 8, title = r'$\rm \rho_{FUV}\ (L_{FUV}>10^{28}\ erg\ s^{-1}\ Hz^{-1}$')
ax2.legend(fontsize = 8, title = r'$\rm SFRD\ (SFR>1\ M_{\odot}\ yr^{-1})$')

ax.set_xlim([6, 15])
ax.set_ylim([23, 26.5])
ax2.set_ylim([-5, -1])

ax.set_xlabel(r'$\rm z $')

ax.set_ylabel(r'$\rm log_{10}(\rho_{FUV}/erg\ s^{-1}\ Hz^{-1}\ Mpc^{-3}) $')
ax2.set_ylabel(r'$\rm log_{10}(\dot{\rho}_{\star}/M_{\odot}\ Mpc^{-3}) $')

fig.savefig('figs/csfd.pdf')
