


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import flare.plt as fplt
from flare.LF import models, evo, plots
from flare.LF.literature import UV


# --- load literature collection
m = getattr(UV, 'FLARES_binned')()
m = getattr(UV, 'FLARES')()


area = 100000. # sq. arcmin


redshift_limits = [5.,10.]
dz = 0.1
log10L_limits = [28., 31.]
dlog10L = 0.1

bin_edges, bin_centres, N = m.N(area = area, redshift_limits = redshift_limits, log10L_limits = log10L_limits, dz = dz, dlog10L = dlog10L)

print(N.shape)
plt.imshow(N, aspect = 'auto', origin = 'lower', vmin = 0, vmax = 100, extent = [*redshift_limits, *log10L_limits])
plt.show()

z, log10L = m.sample(area = area, redshift_limits = redshift_limits, log10L_limits = log10L_limits, dz = dz, dlog10L = dlog10L)

N_sampled,_,_ = np.histogram2d(log10L, z, bins = [bin_edges['log10L'], bin_edges['z']])

print(N_sampled.shape)
plt.imshow(N_sampled, aspect = 'auto', origin = 'lower', vmin = 0, vmax = 100, extent = [*redshift_limits, *log10L_limits])
plt.show()

plt.imshow((N-N_sampled)/N, aspect = 'auto', origin = 'lower', vmin = -0.2, vmax = 0.2, extent = [*redshift_limits, *log10L_limits])
plt.show()



# N, _ = np.histogram(sample, bins = bin_edges)
#
# ax.step(bin_centres, np.log10(N*tot_phi/sample_size)-np.log10(bin_w), where = 'mid', label = 'sampled')
