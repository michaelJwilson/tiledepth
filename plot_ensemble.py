import os
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import pylab as pl

from   astropy.convolution import convolve, Box1DKernel

import speclite.filters

filters = speclite.filters.load_filters('decam2014-*')

tracers = ['bgs', 'lrg', 'elg', 'qso']

colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, axes = plt.subplots(len(tracers), 1, figsize=(7.5, 5 * len(tracers)))

for i, (tracer, color) in enumerate(zip(tracers, colors)):
    dat = fits.open('/project/projectdirs/desi/users/mjwilson/tsnr-ensemble/tsnr-ensemble-{}.fits'.format(tracer))
    nmodel=dat[0].header['NMODEL']

    for band in ['B', 'R', 'Z']:
        wave   = dat['WAVE_{}'.format(band)].data
        dflux  = dat['DFLUX_{}'.format(band)].data[0,:]
        
        smooth = convolve(dflux, Box1DKernel(125), boundary='extend')

        if band == 'B':
            label='{} ({})'.format(tracer, nmodel)
        else:
            label=''
            
        axes[i].plot(wave, smooth, label=label, alpha=0.5, color=color)

    for f in filters:
        axes[i].plot(f._wavelength, f.response * smooth.max() / f.response.max(), '-', c='k', alpha=0.5, lw=0.2)
        
    axes[i].set_ylabel(r'$\langle \Delta F^2 \rangle$')
    axes[i].legend(frameon=False)        

axes[-1].set_xlabel('Wavelength [AA]')
    
pl.show()
