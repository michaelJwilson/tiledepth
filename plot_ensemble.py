import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import pylab as pl

from   astropy.convolution import convolve, Box1DKernel

import speclite.filters

# root  = '/project/projectdirs/desi/users/mjwilson/tsnr-ensemble/' 
root    = '/global/cscratch1/sd/mjwilson/trash/'
# root  = '/global/homes/m/mjwilson/sandbox/desimodel/trunk/data/tsnr/'
# root  = '/global/homes/m/mjwilson/sandbox/desimodel/trunk/data/tsnr/10/'

filters = speclite.filters.load_filters('decam2014-*')

tracers = ['bgs', 'lrg', 'elg', 'qso']

colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

nfig = np.maximum(len(tracers), 2)
fig, axes = plt.subplots(nfig, 1, figsize=(7.5, 5 * nfig))

for i, (tracer, color) in enumerate(zip(tracers, colors)):
    dat = fits.open(root + '/tsnr-ensemble-{}.fits'.format(tracer))
    nmodel=dat[0].header['NMODEL']
    zlo=dat[0].header['ZLO']
    zhi=dat[0].header['ZHI']
    smoothing=dat[0].header['SMOOTH']
    
    for band in ['B', 'R', 'Z']:
        wave   = dat['WAVE_{}'.format(band)].data
        dflux  = dat['DFLUX_{}'.format(band)].data[0,:]

        dwave  = np.mean(np.diff(wave))
        smooth = convolve(dflux, Box1DKernel(np.ceil(100. / dwave)), boundary='extend')

        if band == 'B':
            label='{} ({})'.format(tracer, nmodel)
        else:
            label=''
            
        axes[i].plot(wave, smooth, label=label, alpha=0.5, color=color)

    for f in filters:
        axes[i].plot(f._wavelength, f.response * smooth.max() / f.response.max(), '-', c='k', alpha=0.5, lw=0.2)
        
    if tracer == 'bgs':
        colors = ['blue', 'red', 'green']
        les    = [3727., 6560., 5008.240]

    if tracer == 'lrg':
        colors = ['blue', 'green']
        les    = [3934.8, 5176.]
        
    if tracer == 'elg':
        colors = ['blue', 'green']
        les    = [3727., 5008.240]

    if tracer == 'qso':
        colors = ['blue', 'red', 'cyan', 'green']
        les    = [1215.24, 1240.81, 1397.61, 1908.734] 
        
    for color, le in zip(colors, les):
        axes[i].axvspan(le * (1.+zlo), le * (1.+zhi), alpha=0.1, color=color)

    axes[i].set_xlim(3500., 1.e4)
    axes[i].set_ylabel(r'$\sqrt{\langle \Delta F^2 \rangle}$')
    axes[i].legend(frameon=False)        

axes[-1].set_xlabel('Wavelength [AA]')

fig.suptitle('TSNR Ensemble: {} AA'.format(smoothing), y=0.95)

pl.show()
