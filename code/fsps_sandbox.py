'''

Playing around with FSPS

Author : ChangHoon Hahn 


'''
import fsps
import time
import numpy as np 
import matplotlib.pyplot as plt 

from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors 

def plot_a_spectrum(): 
    t_start = time.time() 
    sps = fsps.StellarPopulation(zcontinuous=1) 
    spec = sps.get_spectrum() 
    print time.time()-t_start, ' seconds'    # takes 3 mins... wtf
    
    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure() 
    sub = fig.add_subplot(111)
    for i in range(spec[1].shape[0]): 
        sub.plot(spec[0], 1.e17*(spec[1][i,:]), c=pretty_colors[i % 10])

    sub.set_xlabel(r'$\lambda$', fontsize=25)
    sub.set_xlim([1000, 10000])
    fig.savefig('aspectrum.png')
    plt.close() 
    return None

def plot_spectrum_tage(): 
    t_start = time.time() 
    sps = fsps.StellarPopulation(zcontinuous=1) 
    print time.time()-t_start, ' seconds'    # takes 3 mins... wtf
    
    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure() 
    sub = fig.add_subplot(111)
    for i, tage in enumerate(np.arange(1., 15., 2.)): 
        wave, spec = sps.get_spectrum(tage=tage) 
        sub.plot(wave, 1.e17*spec, c=pretty_colors[i], label=r'$\mathtt{t_{age} = '+str(tage)+'}$')

    sub.set_xlabel(r'$\lambda$', fontsize=25)
    sub.set_xlim([1000, 10000])
    sub.set_xlabel(r'$f_\lambda$ (unknown scale)', fontsize=25)
    sub.legend(loc='upper left') 
    fig.savefig('spectrum_tage.png', bbox_inches='tight')
    plt.close() 
    return None



if __name__=='__main__': 
    #plot_a_spectrum()
    plot_spectrum_tage()
