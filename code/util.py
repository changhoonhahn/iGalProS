'''
'''
import os 
import fsps
import wget 
import shutil
import numpy as np 
import pickle

from ChangTools.fitstables import mrdfits


def Dir(folder): 
    # Directory address to the code directory 
    if folder == 'code': 
        return os.path.dirname(os.path.realpath(__file__))
    elif folder == 'dat': 
        return os.path.dirname(os.path.realpath(__file__)).split('code')[0]+'dat/'
    elif folder == 'local': # ../local, which is symlinked to a local path that contains big data
        return os.path.dirname(os.path.realpath(__file__)).split('code')[0]+'local/iGalPros/'


def Default_FSPSparams(): 
    ''' Return default FSPS parameters. 
    Mainly used to refresh the parameter lists
    '''
    pickle_file = ''.join([Dir('dat'), 'default_FSPSparams.p' ]) 
    def_param = pickle.load(open(pickle_file, 'rb'))
    return def_param


def Make_Default_FSPSparams(): 
    ''' Save default FSPS parameters to pickle file 
    '''
    pop = fsps.StellarPopulation(zcontinuous=1)

    # load default parameters
    default_params = dict([(k, pop.params[k]) for k in pop.params.all_params])
    default_params['sfh'] = 1
    default_params['dust_type'] = 2
    default_params['add_neb_emission'] = True

    # save to pickle file 
    pickle_file = ''.join([Dir('dat'), 'default_FSPSparams.p' ]) 
    pickle.dump(default_params, open(pickle_file, 'wb'))
    return None


def GetSpectra_SDSS(mjd, plate, fiberid): 
    ''' Get SDSS spectra given mjd, plate, and fiber ID
    '''
    spec_file = '-'.join(['spec', str(plate).zfill(4), str(mjd), str(fiberid).zfill(4)])+'.fits'
    spec_dir = '../local/iGalProS/spectra/'
    if not os.path.isfile(spec_dir+spec_file): 
        print 'https://data.sdss.org/sas/dr8/sdss/spectro/redux/26/spectra/'+str(plate).zfill(4)+'/'+spec_file
        wget.download('https://data.sdss.org/sas/dr8/sdss/spectro/redux/26/spectra/'+str(plate).zfill(4)+'/'+spec_file)

        shutil.move(spec_file, spec_dir+spec_file) 

    spectra = mrdfits(spec_dir+spec_file) 
    return spectra 
