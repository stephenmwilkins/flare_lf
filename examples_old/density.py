
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from flare.photom import  M_to_lum

from flare.LF import models


p = {'M*': -20.6, 'log10phi*': -4.89, 'alpha': -2.28}

LF = models.Schechter(p)

L = M_to_lum(-17.)
print(L/LF.Lstar)

density = LF.density(L)

print(density)


# SFR density

Kappa_UV = 1.15E-28

print(density*Kappa_UV, np.log10(density*Kappa_UV))
