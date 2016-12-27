'''

test model.py: 
    test the FSPS galaxy model 


'''
import fsps 
import numpy as np 
import model as MDL 

# plotting 
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors



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
    test_FSPSgalaxy()
