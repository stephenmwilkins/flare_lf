


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import flare
from flare.LF import completeness, evo, literature, plots
from flare.photom import m_to_flux

cosmo = flare.default_cosmo() # WMAP9

f_limit_deep = flare.photom.m_to_flux(26.)

area = 60.*60.*40.

print('area of Euclid deep: {0} arcmin2'.format(area))

# load in lf parameters
m = getattr(literature, 'Bluetides')()


bin_edges, bin_centres, N = m.N(redshift_limits = [8., 15.], log10L_limits = [27., 30.], dz = 0.1, dlog10L = 0.01)


c = completeness.completeness_erf(bin_centres, f_limit_deep) # calculates the completeness with a flux limit
N = np.multiply(N, c)
n = np.sum(N)
print('density per arcmin: {0:9.2f}'.format(n))
print('number in Eudlid deep: {0:9.1f}'.format(n*area))

plots.evo_plot(bin_edges, N, f_limits=[f_limit_deep])
plt.show()
