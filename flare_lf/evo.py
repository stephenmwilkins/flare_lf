# ---
from astropy.table import Table

import os
this_dir, this_filename = os.path.split(__file__)

import numpy as np
from numpy.random import randn

from scipy.stats import linregress
import scipy.integrate as cp
import scipy.interpolate as cpi
import scipy.special as cps
from scipy.stats import rv_histogram

from mpmath import gammainc

import matplotlib.pyplot as plt

import csv

from flare.photom import flux_to_L, lum_to_flux, M_to_lum, lum_to_M
import flare.core

from .utilities import log10Lnu_to_M, M_to_log10Lnu
from . import models


def dVc(z, cosmo):
    return cosmo.differential_comoving_volume(z).value


def read(model, scheme = 'linear'):

    t = Table.read(f'{this_dir}/data/{model}.ecsv')

    if t.meta['type'] == 'Schechter':
        return Schechter(t, scheme = scheme)

    if t.meta['type'] == 'binned':
        return Binned(t)


def list_data(path):
    return list(map(lambda x: x.replace('.ecsv', ''), os.listdir(f'{this_dir}/data/{path}')))

def list_models_schechter():
    return list_data('models/schechter/')

def list_models_binned():
    return list_data('models/binned/')


def get_models_at_redshift(z = None, model_types = ['schechter', 'binned'], verbose = False):

    models = {}

    if 'schechter' in model_types:
        for model in list_models_schechter():

            model_name = f'models/schechter/{model}'
            m = read(model_name)

            for redshift in m.redshifts:

                redshift = int(redshift)

                if redshift in models.keys():
                    models[redshift].append(model_name)
                else:
                    models[redshift] = [model_name]

    if 'binned' in model_types:
        for model in list_models_binned():

            model_name = f'models/binned/{model}'
            m = read(model_name)

            for redshift in m.redshifts:

                redshift = int(redshift)

                if redshift in models.keys():
                    models[redshift].append(model_name)
                else:
                    models[redshift] = [model_name]

    if verbose:
        for k, v in models.items():
            print(k, v)

    if z:
        return(models[z])
    else:
        return(models)







class evo:


    def phi_binned(self, z, bin_edges, kind = 'linear'):

        if self.model.name == 'binned':
            return(self.model(self.p(z))._phi_binned(bin_edges, kind = kind))
        else:
            return(self.model(self.p(z))._phi_binned(bin_edges))

    def density(self, z, L):

        return(self.model(self.p(z)).density(L))




    def N_singlez(self, area=1., cosmo=False, redshift = 7, log10L_limits=[27.5, 30.], dz=0.2, dlog10L=0.05):

        """ returns a 1D grid of the number of galaxies in each luminosity """

        # calculates the number of galaxies in each bin on a grid defined by redshift_limits, log10L_limits, dz, dlog10L
        # and area based on a luminosity function evolution model.

        area_sm = area  # Area in square arcmin
        area_sd = area_sm / 3600.  # Area in square degrees
        area_sr = (np.pi / 180.) ** 2 * area_sd  # Area in steradian

        if not cosmo: cosmo = flare.core.default_cosmo()

        # Setting the bin edges as well as centres for later operations
        bin_edges = np.arange(log10L_limits[0], log10L_limits[-1] + dlog10L, dlog10L)
        bin_centres = 0.5*(bin_edges[:-1]+bin_edges[1:])

        # Using astropy.cosmology to calculate the volume in each redshift bin
        vol = cosmo.comoving_volume(z+dz/2).value - cosmo.comoving_volume(z-dz/2).value

        N = self.phi_binned(z, bin_edges['log10L']) * vol * area_sr


        return bin_edges, bin_centres, N



    def N(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27.5, 30.], dz=0.05, dlog10L=0.05, per_arcmin=False, return_volumes=False):

        """ returns a 2D grid of the number of galaxies in each luminosity/redshift bin """


        # calculates the number of galaxies in each bin on a grid defined by redshift_limits, log10L_limits, dz, dlog10L
        # and area based on a luminosity function evolution model.

        area_sm = area  # Area in square arcmin
        area_sd = area_sm / 3600.  # Area in square degrees
        area_sr = (np.pi / 180.) ** 2 * area_sd  # Area in steradian

        if not cosmo: cosmo = flare.core.default_cosmo()

        # Setting the bin edges as well as centres for later operations
        bin_edges = {'log10L': np.arange(log10L_limits[0], log10L_limits[-1] + dlog10L, dlog10L),
                     'z': np.arange(redshift_limits[0], redshift_limits[-1] + dz, dz)}
        bin_centres = {'log10L': bin_edges['log10L'][:-1] + dlog10L / 2., 'z': bin_edges['z'][:-1] + dz / 2.}

        # Using astropy.cosmology to calculate the volume in each redshift bin
        volumes = np.asarray([cp.quad(dVc, bin_edges['z'][i - 1], bin_edges['z'][i], args=cosmo)[0] for i in
                              range(1, len(bin_edges['z']))])




        N = np.zeros(( len(bin_centres['log10L']), len(bin_centres['z']) ) )

        for i, (z, vol) in enumerate(zip(bin_edges['z'], volumes)):

            N[:,i] = self.phi_binned(z, bin_edges['log10L']) * vol * area_sr


        if per_arcmin:
            N /= area_sm

        if return_volumes:
            return bin_edges, bin_centres, volumes, N

        else:
            return bin_edges, bin_centres, N




    def sample(self, area=1., cosmo=False, redshift_limits=[8., 15.], log10L_limits=[27.5, 30.], dz=0.05, dlog10L=0.05, per_arcmin=False, return_volumes=False):

        """ returns a sample of redshifts and UV luminosities """

        bin_edges, bin_centres, N = self.N(area = area, cosmo = cosmo, redshift_limits = redshift_limits, log10L_limits = log10L_limits, dz = 0.1, dlog10L = 0.1)

        z = np.array([])
        log10L = np.array([])

        for i, (z_l, z_u) in enumerate(zip(bin_edges['z'][:-1], bin_edges['z'][1:])):

            # --- sample size in this bin
            sample_size = int(np.sum(N[:, i]))

            # --- sample redshift uniformly across the bin

            z = np.append(z, np.random.uniform(low = z_l, high = z_u, size = sample_size))

            hist_dist = rv_histogram((N[:, i], bin_edges['log10L']))

            log10L = np.append(log10L, hist_dist.rvs(size = sample_size))

        return z, log10L





class parameterised:

    """ used to provide LF parameters assuming a (pre-computed) linear fit to the available parameters """

    def __init__(self, scheme = 'linear'):
        # lp is a dictionary of the parameters of the linear evolution model

        self.scheme = scheme

        if scheme == 'linear':

            self.linear.calculate_linear_evolution_coeffs(self)
            self.p = lambda x: self.linear.p(self, x)

        if scheme == 'interp':

            self.p = lambda x: self.interp.p(self, x)


    class linear:


        def p(self, z):

            """ return the parameters that are fed to the model class to generate the function that gives phi = f(log10) """
            """ *** in this case the parameters are the Schechter/DPL parameters """

            return {param: self.lp[param][0] * (z - self.z_ref) + self.lp[param][1] for param in self.lp}


        def parameters_line(self, zr = [6.,13.]):

            p = {}
            for param in self.lp:
                p[param] = [self.lp[param][0] * (zr[0] - self.z_ref) + self.lp[param][1], self.lp[param][0] * (zr[1] - self.z_ref) + self.lp[param][1]]

            return zr, p


        def calculate_linear_evolution_coeffs(self, zr=[5., 15.], z_ref=5.):
            # Function that calculates the linear evolution coeffs
            # returns a dictionary of linear model coefficients and goodness of fit

            s = (np.array(self.redshifts) >= zr[0]) & (np.array(self.redshifts) <= zr[1])

            z_mod = np.array(self.redshifts)[s] - z_ref
            alpha_mod = np.array(self.alpha)[s]
            log10phi_mod = np.array(self.log10phi_star)[s]
            M_mod = np.array(self.M_star)[s]

            # The output contains full linregress output (0th and 1st element contain the slope and intercept respectively)
            fit_alpha = linregress(z_mod, alpha_mod)
            fit_log10phi = linregress(z_mod, log10phi_mod)
            fit_M = linregress(z_mod, M_mod)

            if self.model.name == 'DPL':
                beta_mod = np.array(self.beta)[s]
                fit_beta = linregress(z_mod, beta_mod)
                lp = {'alpha': fit_alpha, 'beta': fit_beta, 'log10phi*': fit_log10phi, 'M*': fit_M}

            else:
                lp = {'alpha': fit_alpha, 'log10phi*': fit_log10phi, 'M*': fit_M}

            self.z_ref = z_ref
            self.lp = lp



    class interp:

        """ used to provide LF parameters interpolating the available parameters. DEVELOPMENT: This should be adapted to use different interpolation schemes and ranges of parameters (e.g. including only neighbouring points). """


        def p(self, z=8.):
            # interpolates parameters as a function of z
            # returns a dictionary of the Schechter function parameters for given redshift(s)

            z_mod = self.redshifts
            alpha_mod = self.alpha
            log10phi_mod = self.log10phi_star
            log10M_mod = self.M_star

            if self.model.name == 'Double Power Law':
                beta_mod = self.beta
                p = {'alpha': np.interp(z, z_mod, alpha_mod), 'log10phi*': np.interp(z, z_mod, log10phi_mod),
                     'M*': np.interp(z, z_mod, log10M_mod), 'beta': np.interp(z, z_mod, beta_mod)}
            else:
                p = {'alpha': np.interp(z, z_mod, alpha_mod), 'log10phi*': np.interp(z, z_mod, log10phi_mod),
                     'M*': np.interp(z, z_mod, log10M_mod)}

            return p




class Schechter(evo, parameterised):

    def __init__(self, t, scheme = 'linear'):

        self.name = t.meta['name']
        self.model = models.Schechter
        self.t = t
        self.redshifts = self.t['redshift']
        self.alpha = self.t['alpha']
        self.log10phi_star = self.t['log10phi*']
        self.M_star = self.t['M*']

        super().__init__(scheme = scheme)





class Binned(evo):

    def __init__(self, t):

        self.name = t.meta['name']
        self.model = models.binned
        self.t = t

        self.M = {}
        self.log10L = {}

        self.phi = {}
        self.log10phi = {}
        self.log10phi_dex = {}
        self.log10phi_mag = {}


        # --- if data is redshift, log10L, phi ... this is most useful I think

        self.redshifts = list(set(self.t['redshift'].data))
        # --- make sure that redshift list is monotonically increasing
        if self.redshifts[0]>self.redshifts[1]:
            self.redshifts = self.redshifts[::-1]

        # --- extract luminosities or magnitudes and number densities, making conversions where necessary.
        for z in self.redshifts:

            if 'log10L' in self.t.colnames:
                self.log10L[z] = self.t['log10L'][self.t['redshift']==z].data
                self.M[z] = log10Lnu_to_M(self.log10L[z])

                if 'log10phi' in self.t.colnames:
                    self.log10phi_dex[z] = self.t['log10phi'][self.t['redshift']==z].data
                elif 'phi' in self.t.colnames:
                    self.log10phi_dex[z] = np.log10(self.t['phi'][self.t['redshift']==z].data)
                else:
                    print('WARNING: one column should be log10phi or phi')

                self.log10phi_mag[z] = self.log10phi_dex[z] + np.log10(0.4)

            elif 'M' in self.t.colnames:
                self.M[z] = self.t['M'][self.t['redshift']==z].data
                self.log10L[z] = M_to_log10Lnu(self.M[z])

                if 'log10phi' in self.t.colnames:
                    self.log10phi_mag[z] = self.t['log10phi'][self.t['redshift']==z].data
                elif 'phi' in self.t.colnames:
                    self.log10phi_mag[z] = np.log10(self.t['phi'][self.t['redshift']==z].data)
                else:
                    print('WARNING: one column should be log10phi or phi')

                self.log10phi_dex[z] = self.log10phi_mag[z] - np.log10(0.4)

            else:
                print('no luminosity/magnitude column found, use log10L or M')

            self.log10phi[z] = self.log10phi_dex[z]
            self.phi[z] = 10**self.log10phi[z]

        self.phi_log10L = {} # arrays of phi values at different redshifts for each luminosity bin

        for i, log10L in enumerate(self.log10L):
            self.phi_log10L[log10L] = np.zeros(len(self.redshifts))
            for j, z in enumerate(self.redshifts):
                self.phi_log10L[log10L][j] = self.phi[z][i]




    def p(self, z):

        """ return the parameters that are fed to the model class to generate the function that gives phi = f(log10) """
        """ *** in this case the parameters are log10L and phi bin values """

        p = {'log10L': self.log10L}

        p['phi'] = np.zeros(len(self.log10L))
        for i, log10L in enumerate(self.log10L):
            p['phi'][i] = np.interp(z, self.redshifts, self.phi_log10L[log10L]) # 1D linear interpolation

        return p
