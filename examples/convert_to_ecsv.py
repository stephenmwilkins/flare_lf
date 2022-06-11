

import numpy as np
from astropy.table import Table, Column
import astropy.units as u
from flare.LF.literature import UV








for model in ['FLARES_intrinsic', 'FLARES_intrinsic_binned', 'Bluetides', 'Ma2019', 'Mason2015', 'Yung2018']



m = UV.FLARES_binned()

print(m.redshifts)
print(m.log10L)

t = Table()


redshift = np.array([])
log10L = np.array([])
phi = np.array([])

for z in m.redshifts:
    redshift = np.append(redshift, z*np.ones(len(m.log10L))) # only works if only one set of luminosities
    log10L = np.append(np.round(log10L, 2), m.log10L)
    phi = np.append(phi, m.phi[z])

t.add_column(Column(data = redshift, name = 'redshift', description = 'redshift'))
t.add_column(Column(data = log10L, name = 'log10L', description = 'log10(luminosity/erg/s/Hz)'))
t.add_column(Column(data = phi, name = f'phi', description = 'log10(phi/Mpc^3/dex)'))

t.meta['name'] = m.name
t.meta['redshifts'] = m.redshifts
t.meta['model type'] = m.type
t.meta['type'] = 'binned'

t.write('../flare_lf/data/flares_binned.ecsv', format = 'ascii.ecsv', overwrite=True)


# t.add_column(Column(data = m.log10L, name = 'log10L', description = 'log10(luminosity/erg/s/Hz)'))
# print(t['log10L'])
#
# for z in m.redshifts:
#     t.add_column(m.phi[z], name = f'phi_{z}')



# q = 10 * u.erg / u.s / u.Hz
# print(np.log10(q))
#
# q = m.log10L * u.dex(u.erg / u.s / u.Hz)
#
# print(q)
# phi = 1.0 * u.Mpc**-3 / u.dex
#
# print(phi)
#
# print(phi.to('Mpc^-3 / mag'))


# t.add_column(m.log10L, name = 'log10L')
