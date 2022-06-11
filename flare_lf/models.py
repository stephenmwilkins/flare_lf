

import numpy as np

import scipy.integrate as cp
import scipy.interpolate as cpi
import scipy.special as cps
# from scipy.special import gammaincc
from mpmath import gammainc

from flare.photom import flux_to_L, lum_to_flux, M_to_lum, lum_to_M

# these classes




def bin_centres(bin_edges):
    return (bin_edges[:-1] + bin_edges[1:])/2.

def bin_width(bin_edges):
    return bin_edges[1] - bin_edges[0]





class binned:

    name = 'binned'

    def __init__(self, p):

        self.log10L = p['log10L']
        self.phi = p['phi']

    def _phi_binned(self, bin_edges, kind = 'linear'):

        """ return the LF on an arbitrary log10(luminosity) grid
            bin_edges = log10L
        """

        f = cpi.interp1d(self.log10L, self.phi*bin_width(bin_edges), kind = kind, fill_value='extrapolate')

        return f(bin_centres(bin_edges))

        # return np.interp(bin_centres(bin_edges), self.log10L, self.phi*bin_width(bin_edges))  # should we interpolate phi or log10phi?



    def density(self, L):

        if np.log10(L) < self.log10L[0]:
            print('WARNING: L is smaller than lowest bin. This will not work')

        # --- define new bin edges from the minimum to the maximum
        bin_edges = np.linspace(np.log10(L), self.log10L[-1], 100)

        # --- get new phi
        phi = self._phi_binned(bin_edges)

        return np.sum(phi*10**bin_centres(bin_edges))


class Schechter:

    name = 'Schechter'


    def __init__(self, p):

        self.parameters = ['alpha', 'log10phi*', 'log10L*']

        self.p = p
        self.alpha = p['alpha']

        if 'M*' in p.keys():
            self.Lstar = M_to_lum(p['M*'])
            self.Mstar = p['M*']

        self.phistar = 10 ** p['log10phi*']


    def _phif(self, x):
        return x ** (self.alpha) * np.exp(-x)

    # def _phif(self, x):
    #     return np.log(10)*(x)**(self.alpha+1)*np.exp(-x)/(x*np.log(10))

    # def _phif(self, x):
    #     dM = 2.5*np.log10(x)
    #     return (0.4*np.log(10))*(10**(0.4*dM))**(self.alpha+1)*np.exp(-10**(0.4*dM))

    def _phi(self, L):

        return self.phistar * self._phif(L/self.Lstar)

    def density(self, L):

        """ get the density down to some limit """

        return self.phistar * self.Lstar * np.float(gammainc(self.alpha + 2, L/self.Lstar))


    def _density(log10phi_star, L_star, alpha, L_limit):

        """ get the density down to some limit. This is used for re-sampling """

        return 10**log10phi_star * L_star * np.float(gammainc(alpha + 2, L_limit/L_star))


    def _phi_binned(self, bin_edges):

        """ integrate the LF between the bin edges to get the number density of galaxies in the bin """

        y = np.zeros(len(bin_edges)-1)

        for i, (a,b) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):

            y[i] = self.phistar * cp.quad(self._phif, (10 ** a)/self.Lstar, (10 ** b)/self.Lstar)[0]

        return y




class DPL:

    name = 'DPL'

    def __init__(self, p):

        print('WARNING: not yet implemented')

        self.alpha = p['alpha']
        self.beta = p['beta']
        self.Lstar = M_to_lum(p['M*'])
        self.Mstar = p['M*']
        self.phistar = 10 ** p['log10phi*']


    def _phif(self, x):

        return 2.*10./(x**(-self.alpha) + x**(-self.beta))

    def _phi_M(self, M):

        return self.phistar / (10**(0.4*(M-self.Mstar)*(self.alpha+1)) + 10**(0.4*(M-self.Mstar)*(self.beta+1)))


    def _phi(self, L):

        return self.phistar * self._phif(L/self.Lstar)


    def _phi_binned(self, bin_edges):

        """ integrate the LF between the bin edges to get the number density of galaxies in the bin """

        y = np.zeros(len(bin_edges)-1)

        for i, (a,b) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):

            y[i] = self.phistar * cp.quad(self._phif, (10 ** a)/self.Lstar, (10 ** b)/self.Lstar)[0]

        return y
