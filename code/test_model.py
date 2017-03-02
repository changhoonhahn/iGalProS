'''

test model.py: 
    test the FSPS galaxy model 


'''
import fsps 
import numpy as np 
import model as MDL 
import util as UT 

# plotting 
import matplotlib.pyplot as plt
from matplotlib import gridspec
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
from ChangTools.fitstables import mrdfits 


def test_iSEDfit():
    ''' Condensing fsps_galrecon.ipynb and spectra.ipynb into a small self contained test. 
    '''
    # lambda_SDSS bands
    sdss_bands = fsps.find_filter('sdss')
    L_sdss = [fsps.get_filter(sdss_band).lambda_eff/1.e4 for sdss_band in ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']]

    # import iSEDfit data into convenient dictionaries 
    iSEDfit = mrdfits(UT.Dir('local')+'NSA_iSEDfit/nsa_v1_2_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits.gz')

    gal_ids = [53999, 18965, 81806] # hardcoded IDs of some randomly chosen galaxies 
    
    iSED_gals = [{} for i in range(len(gal_ids))]
    for i, gal_dict in enumerate(iSED_gals): 
        for key in iSEDfit.__dict__.keys(): 
            gal_dict[key] = getattr(iSEDfit, key)[gal_ids[i]]
        
        # load in the spectra obtained manually from SDSS SpecObjAll
        sdss_spec = mrdfits(UT.Dir('local')+'gal'+str(i+1)+'.fits')
        gal_dict['lambda'] = 10**sdss_spec.loglam # Angstrom
        gal_dict['spectra'] = sdss_spec.flux
         
    # compare spectra and photometry
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(12,9))
    bkgd = fig.add_subplot(111, frameon=False) 
    gs = gridspec.GridSpec(2,2, width_ratios=[2,1]) 
    
    model = MDL.FSPSgalaxy()
    for i, gal_dict in enumerate(iSED_gals): 
        wa, sp = model._getSpectrum_iSEDfit(gal_dict, units='cgs')
        
        sub_spec = plt.subplot(gs[2*i]) 
        sub_spec.plot(gal_dict['lambda']/1e4, gal_dict['spectra'], c='k') 
        if i == 0: 
            sub_spec.plot(wa/1e4, sp*1e17/1.8, c=pretty_colors[3], ls='--', lw=3) 
            sub_spec.set_ylim([0., 30])
        elif i == 1: 
            sub_spec.plot(wa/1e4, sp*1e17/4.74, c=pretty_colors[3], ls='--', lw=3) 
            sub_spec.set_ylim([0., 200])
        sub_spec.set_xlim([0.38, 0.92]) 
        sub_spec.set_ylabel('$erg/s/cm^2/A$')
        
        sub_photo = plt.subplot(gs[2*i+1]) 
        photo_ab = -2.5*np.log10(gal_dict['maggies'][1:-1])
        sub_photo.scatter(L_sdss, photo_ab, c='k', s=60, label='Obs.') 
        model_photo = model._getPhotometry_iSEDfit(gal_dict)
        sub_photo.scatter(L_sdss, model_photo, c=pretty_colors[3], s=60, lw=0, label='FSPS Model') 
        # x-axis 
        sub_photo.set_xlim([0.3, 1.]) 
        sub_photo.set_xticks([0.3, 0.5, 0.7, 0.9]) 
        # y-axis
        sub_photo.set_ylim([1.1*np.max(photo_ab), 0.9*np.min(photo_ab)])
        sub_photo.yaxis.tick_right()
        sub_photo.set_ylabel('AB Magnitude') 
        sub_photo.yaxis.set_label_position("right")
        if i == 1: 
            sub_photo.legend(loc='lower right', scatterpoints=1, markerscale=1.5, handletextpad=-0.5)

    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_xlabel('Wavelength $\lambda\; (\mu m)$', fontsize=25)

    fig.savefig('iSEDfit.model_obs_comparison.png', bbox_inches='tight')  
    return None




def test_FSPSgalaxy(): 
    ''' Test FSPS galaxy model by comparing it to iSEDfit derived galaxy properties. 
    Using iSEDfit derived galaxy properties, generates galaxy spectrum and photometry. 
    '''
    model = MDL.FSPSgalaxy()
    
    # iSEDfit derived galaxy properties
    isedfit_dict = [
            {'sfh':1, 'tau':1., 'logzsol': np.log10(0.84), 'tage':3.97, 'mass':10**9.16, 'zred': 0.0309},       # ID: 053999
            {'sfh':1, 'tau':0.6, 'logzsol': np.log10(0.39), 'tage':12.66, 'mass':10**11.30, 'zred': 0.0273},    # ID: 133295 
            {'sfh':1, 'tau':2.6, 'logzsol': np.log10(0.74), 'tage':7.94, 'mass':10**10.82, 'zred': 0.0447},     # ID: 081806 
            {'sfh':1, 'tau':4.8, 'logzsol': np.log10(0.32), 'tage':4.44, 'mass':10**8.38, 'zred': 0.0251}       # ID: 011418
            ]

    band_wavelength = [fsps.get_filter(band).lambda_eff/1e4 for band in ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']]
    
    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(12,9)) 
    bkgd = fig.add_subplot(111, frameon=False) 

    for ii, isedfit in enumerate(isedfit_dict): 
        model.updateParams(**isedfit)
        w, spec = model.getSpectrum(units='ABmag')
        photo = model.getPhotometry(units='ABmag', filters='sdss')
    
        sub = fig.add_subplot(2,2,ii+1)
        sub.plot(w/1e4, spec, c='k', lw=2, alpha=0.5, label='FSPS spectrum') 
        sub.scatter(band_wavelength, photo, c=pretty_colors[3], lw=0, s=20, label='FSPS photometry') 

        # axes
        sub.set_xlim([0.1, 1.0])
        sub.set_xscale("log")
        sub.set_xticks([0.2, 0.5, 1.])
        sub.set_ylim([np.ceil(spec[np.where((w/1e4 > 0.1) & (w/1e4 < 1.0))].max()), np.floor(spec[np.where((w/1e4 > 0.1) & (w/1e4 < 1.0))].min())])
        if ii == 0: 
            sub.legend(loc='upper left', scatterpoints=1) 

    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_xlabel('Wavelength $\lambda\; (\mu m)$', fontsize=25)
    bkgd.set_ylabel('AB magnitude', fontsize=25)

    fig.savefig('testing.png', bbox_inches='tight')  
    return None



if __name__=='__main__': 
    test_iSEDfit()
    #test_FSPSgalaxy()
