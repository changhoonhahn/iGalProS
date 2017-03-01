'''
'''
import os 
import fsps
import numpy as np 
import pickle


def Dir(folder): 
    # Directory address to the code directory 
    if folder == 'code': 
        return os.path.dirname(os.path.realpath(__file__))
    elif folder == 'dat': 
        return os.path.dirname(os.path.realpath(__file__)).split('code')[0]+'dat/'


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
