


import numpy as np

geo = 4.*np.pi*(100.*10.*3.0867*10**16)**2 # factor relating the L to M in cm^2
log10geo = np.log10(geo)


def fnu_to_m(fnu):

    return -2.5*np.log10(fnu/1E9) + 8.9 # -- assumes flux in nJy

def m_to_fnu(m):

    return 1E9 * 10**(-0.4*(m - 8.9)) # -- flux returned nJy

def fnu_to_Lnu(fnu, cosmo, z):

    """ convert flux to luminosity including the band stretching effect """

    return fnu*(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2)/(1E9 * 1E23 * (1.+z))

def Lnu_to_fnu(Lnu, cosmo, z):

    """ convert luminosity to flux including the band stretching effect """

    return 1E9 * 1E23 * Lnu * (1.+ z)/(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2)


def Lnu_to_M(Lnu):

    return -2.5*np.log10(Lnu/geo)-48.6

def M_to_Lnu(M):

    return 10**(-0.4*(M+48.6)) * geo


def log10Lnu_to_M(log10Lnu):

    return -2.5*log10Lnu-log10geo-48.6

def M_to_log10Lnu(M):

    return -0.4*(M+48.6) + log10geo




def DM(cosmo, z):
    luminosity_distance = cosmo.luminosity_distance(z).to('pc').value
    return 5*np.log10(luminosity_distance/(np.sqrt(1.+z)*10.))

def M_to_m(M, cosmo, z):
    return M + DM(z, cosmo = cosmo)

def m_to_M(m, cosmo, z):
    return m - DM(z, cosmo = cosmo)
