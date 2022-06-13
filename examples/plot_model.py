



import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import flare_lf.evo as evo
import flare_lf.plots as plots



# --- plot a particular model

model = 'flares'

# fig, ax = plots.plot_lf_evo(model) # plot basic LF for model, including Schechter function is available
# fig, ax = plots.plot_lf_evo(model, lum_type = 'M') # as above but using magnitudes


# --- plot all models at a particular redshift

z = 7


models = evo.get_models_at_redshift(z = z, model_types = ['binned']) # list the models availabe at a particular redshift (or al redshifts)

print(models)

fig, ax = plots.plot_lf(z, models) # as above but using magnitudes
