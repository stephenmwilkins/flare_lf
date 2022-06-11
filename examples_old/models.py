



import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from flare.LF import models, evo, literature



bin_edges = np.arange(28, 29.1, 0.1)

# --- simple test

p = {'alpha': -4.215248131716793, 'log10phi*': -6.680373396887468, 'M*': -20.070471399056707}

LF = models.Schechter(p)

print(LF._phi_binned(bin_edges))

print(LF.name)

# ---

p = {'log10L': [28, 28.5, 29.0, 29.5], 'phi': [1, 2, 4, 5][::-1]}

LF = models.binned(p)

print(LF._phi_binned(bin_edges))
