


import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr # provides wider range of cmaps, see https://cmasher.readthedocs.io

import flare_lf.evo as evo


model_names = evo.list_models_schechter()

fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize = (3,6))
fig.subplots_adjust(left=0.2, bottom=0.05, top=0.85, right=0.95, wspace=0.0, hspace=0.0)


axes[-1].set_xlabel(r'$\rm z $')
axes[0].set_ylabel(r'$\rm M^{\star} $')
axes[1].set_ylabel(r'$\rm \log_{10}\phi^{\star} $')
axes[2].set_ylabel(r'$\rm \alpha $')


# axes[0].set_ylim([-21.95, -19.55])
# axes[1].set_ylim([-6.45, -2.05])
# axes[2].set_ylim([-2.95, -1.55])

# for ax in axes:
#     ax.set_xlim([4.5,14.])


# ---
cmap = 'cmr.neon'
colors = cmr.take_cmap_colors(cmap, len(model_names))
marker_styles = ['o','^','h','d','*','v','p','s']*2
s = 10
alpha = 0.8

for model_name, color, marker_style in zip(model_names, colors, marker_styles):

    m = evo.read(f'models/schechter/{model_name}')

    axes[0].scatter(m.redshifts, m.M_star, marker = marker_style, c = [color], alpha = alpha, s=s, label = m.name)
    axes[1].scatter(m.redshifts, m.log10phi_star, marker = marker_style, c = [color], alpha = alpha, s=s)
    axes[2].scatter(m.redshifts, m.alpha, marker = marker_style, c = [color], alpha = alpha, s=s)

axes[0].legend(fontsize = 6, bbox_to_anchor=(0.0, 1.05), loc = 'lower left')

# fig.savefig('figs/compare_schechter.pdf')

plt.show()
