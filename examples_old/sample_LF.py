


import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare.photom import  M_to_lum, lum_to_M
import flare.plt as fplt
from flare.LF import models, evo, plots
from flare.LF.literature import UV, SFR



M_limit = -17.7
L_limit = M_to_lum(M_limit)
N = 1000

model = 'Mcleod16'
m = getattr(UV, model)()

z = 9.

d = m.sample_density(z, L_limit)

d = np.log10(d)

med = np.percentile(d, 50)

# print(f' {np.percentile(d, 16):.2f} {med:.2f} {np.percentile(d, 84):.2f}')
print(f' {med:.2f} {np.percentile(d, 16)-med:.2f} {np.percentile(d, 84)-med:.2f}')






# iz = m.redshifts.index(z)
#
# print(iz)
#
#
# v = randn(N)
# v[v<0] *= -m.alpha_err[iz][0]
# v[v>0] *= m.alpha_err[iz][1]
# alphas = m.alpha[iz] + v
# # print(alphas)
#
# v = randn(N)
# v[v<0] *= -m.M_star_err[iz][0]
# v[v>0] *= m.M_star_err[iz][1]
# M_stars = m.M_star[iz] + v
# L_stars = M_to_lum(M_stars)
# # print(M_stars)
# # print(L_stars)
#
# v = randn(N)
# v[v<0] *= -m.phi_star_err[iz][0]
# v[v>0] *= m.phi_star_err[iz][1]
# log10phi_stars = m.phi_star[iz] + v
# # print(log10phi_stars)
#
# d = np.array([m.model._density(log10phi_star, L_star, alpha, L_limit) for log10phi_star, L_star, alpha in zip(log10phi_stars, L_stars, alphas)])
#
# d = np.log10(d)
#
# med = np.percentile(d, 50)
#
# print(f' {np.percentile(d, 16):.2f} {med:.2f} {np.percentile(d, 84):.2f}')
#
# print(f' {med:.2f} {np.percentile(d, 16)-med:.2f} {np.percentile(d, 84)-med:.2f}')
#






# 25.31−0.14+0.1
