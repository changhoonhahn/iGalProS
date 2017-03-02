'''


Galaxy spectrum/photometry model 


Author(s): ChangHoon Hahn


'''

import fsps
import numpy as np 
from sedpy.observate import getSED
from sedpy.observate import load_filters
from astropy.cosmology import FlatLambdaCDM

import util as UT


class FSPSgalaxy(object): 
    ''' FSPS galaxy object. Given a set of parameters that specify 
    the composite stellar population, generate spectra and photometry 
    of galaxies. 
    '''
    def __init__(self, **params):  
        '''
        '''
        # initiate FSPS Composite Stellar Population object
        self.pop = fsps.StellarPopulation(zcontinuous=1)
    
        # update the params with given params
        self.updateParams(**params)

        # some useful constants 
        self.lsun = 3.846e33 # erg/s
        self.pc = 3.085677581467e18 # cm
        self.lightspeed = 2.998e18 # AA/s
        self.to_cgs = self.lsun/(4.*np.pi * (self.pc * 10)**2)

        # default cosmology
        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    def _getSpectrum_iSEDfit(self, iseddict, units='cgs'): 
        ''' Given iSEDfit derived galaxy property dictionary, output galaxy spectra
        '''
        assert set(iseddict.keys()).issuperset(set(['tau', 'av', 'mu', 'age', 'mstar', 'zmetal']))

        self.pop.params['sfh'] = 1
        self.pop.params['tau'] = iseddict['tau']

        self.pop.params['dust_type'] = 2 
        self.pop.params['dust2'] = 0.92 * iseddict['av'] * iseddict['mu'] # ln(10) * 0.4 * A(V)

        self.pop.params['logzsol'] = np.log10(iseddict['zmetal']/0.019)

        self.pop.params['add_neb_emission'] = True
        self.pop.params['add_neb_continuum'] = True
        
        w, sp = self.pop.get_spectrum(tage = iseddict['age'], peraa=False) 

        a = 1. + iseddict['z'] # scale factor

        wave = w * a

        mass = 10**iseddict['mstar'] / self.pop.stellar_mass # stellar mass scaling 
        spec = mass * sp * a 
    
        d_lum = self.cosmo.luminosity_distance(iseddict['z']).value
        d_factor = (d_lum * 1e5)**2

        if units == 'cgs': 
            spec *= self.to_cgs / d_factor * self.lightspeed / wave**2 # in erg/s/cm^2/AA
        elif units == 'ABmag': # AB magnitude (think twice before messing with units) 
            spec = -2.5 * np.log10(spec * self.to_cgs / d_factor / 1e3 / (3631 * 1e-26))
        else: 
            raise NotImplementedError()
        
        return wave, spec 

    def _getPhotometry_iSEDfit(self, iseddict): 
        ''' Returns AB magnitude photometry
        '''
        assert set(iseddict.keys()).issuperset(set(['tau', 'av', 'mu', 'age', 'mstar', 'zmetal']))

        self.pop.params['sfh'] = 1
        self.pop.params['tau'] = iseddict['tau']

        self.pop.params['dust_type'] = 2 
        self.pop.params['dust2'] = 0.92 * iseddict['av'] * iseddict['mu'] # ln(10) * 0.4 * A(V)

        self.pop.params['logzsol'] = np.log10(iseddict['zmetal']/0.019)

        self.pop.params['add_neb_emission'] = True
        self.pop.params['add_neb_continuum'] = True
        
        w, sp = self.pop.get_spectrum(tage = iseddict['age'], peraa=False) 

        a = 1. + iseddict['z'] # scale factor

        wave = w * a

        mass = 10**iseddict['mstar'] / self.pop.stellar_mass # stellar mass scaling 
        
        d_lum = self.cosmo.luminosity_distance(iseddict['z']).value
        d_factor = (d_lum * 1e5)**2
        
        spec = mass * sp * a * self.to_cgs / d_factor * self.lightspeed / wave**2 # in erg/s/cm^2/AA

        return getSED(wave, spec, filterlist=load_filters(['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']))

    def getSpectrum(self, units='ABmag', **params): 
        ''' insert description here after things are finalized
        '''
        if params != {}: 
            print 'updating params'
            self.updateParams(**params)

        # get spectrum 
        wave, spec = self.pop.get_spectrum(tage=self.tage, peraa=False)
        mass = self.mass / self.pop.stellar_mass
        spec *= mass 
        
        # cosmology stuff
        if self.zred is None: 
            raise ValueError("No 'zred' specified") 
        d_lum = self.cosmo.luminosity_distance(self.zred).value # d_luminosity
        d_factor = (d_lum * 1e5)**2
        a = (1. + self.zred) # scale factor

        if units == 'ABmag': # AB magnitude
            spec = -2.5 * np.log10(spec * self.to_cgs * a / d_factor / 1e3 / (3631 * 1e-26)) 
        elif units == 'maggies': # maggies
            spec *= self.to_cgs * a / d_factor / 1e3 / (3631 * 1e-26)
        elif units == 'cgs': 
            spec *= self.to_cgs * a / d_factor * self.lightspeed / wave**2 # in erg/s/cm^2/AA
        else: 
            raise ValueError("In unspecified units")

        return wave, spec

    def getPhotometry(self, units='ABmag', filters='sdss', **params): 
        ''' insert description here after things are finalized
        '''
        if params != {}: 
            print 'updating params'
            self.updateParams(**params)

        # get photometry of specified bands
        if filters == 'sdss': 
            sdss_mags = self.pop.get_mags(tage=self.tage, 
                    bands=['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z'])
        else: 
            raise NotImplementedError('filter bands not implemented') 
        mass = self.mass / self.pop.stellar_mass
        
        if self.zred is None: 
            raise ValueError("No 'zred' specified") 
        d_lum = self.cosmo.luminosity_distance(self.zred).value
        d_factor = (d_lum * 1e5)**2
        a = (1. + self.zred) 

        if units == 'ABmag': # *apparent* AB magnitude
            photo = -2.5*np.log10(mass * 10**(-0.4 * sdss_mags) * a / d_factor)
        elif units == 'maggies': # maggies
            photo = 10**mass * 10**(-0.4 * sdss_mags) * a / d_factor
        else: 
            raise ValueError("In unspecified units")
        return photo 

    def updateParams(self, **params): 
        '''
        '''
        # parse through specified parameters 
        if 'tage' in params.keys(): 
            self.tage = params['tage']
            params.pop('tage')
        if 'mass' in params.keys(): 
            self.mass = params['mass']
            params.pop('mass')
        if 'zred' in params.keys(): 
            self.zred = params['zred']
        else: 
            self.zred = None 
        
        # if default params haven't been loaded yet
        if 'default_params' not in self.__dict__.keys(): 
            self.default_params = UT.Default_FSPSparams()
        
        # reset parameters
        for k in self.pop.params.all_params:
            self.pop.params[k] = self.default_params[k]
        
        for k in params.keys():  
            self.pop.params[k] = params[k]

        return None
