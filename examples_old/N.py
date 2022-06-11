


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from flare.LF import models, evo, plots
from flare.LF.literature import UV


# --- load literature collection
m = getattr(UV, 'FLARES_binned')()
m = getattr(UV, 'FLARES')()

# --- get schechter function parameters, or binned LF at arbitrary redshift
print(m.p(8.5)) # by default uses a linear fit

# m = getattr(literature, 'Bluetides')(scheme='interp') # use (linear) interpolation instead
# print(m.p(8.5))


# --- get binned LF
bin_edges = np.arange(28, 29.1, 0.1)
print(m.model(m.p(8.5))._phi_binned(bin_edges))
print(m.phi_binned(8.5, bin_edges)) # alias to above, should be identical

# --- get 2D density
bin_edges, bin_centres, N = m.N(redshift_limits = [5., 15.], log10L_limits = [27., 30.], dz = 0.1, dlog10L = 0.1)

print(np.sum(N))

plots.evo_plot(bin_edges, N)
plt.show()
