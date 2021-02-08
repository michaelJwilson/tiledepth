import pylab as pl
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from astropy.table import Table

fig, axes = plt.subplots(1,3,figsize=(15,5))

for i, band in enumerate(['b', 'r', 'z']):
    dat = Table.read('cframe-{}0-00067972.fits'.format(band), 'SCORES')
    # dat.pprint()

    for tracer in ['LRG', 'ELG', 'QSO']:
        tsnr = dat['{}TSNR_{}'.format(tracer,band.upper())]
        axes[i].hist(tsnr, bins=25, label=tracer)

    axes[i].set_ylabel(band)
    axes[i].legend()
pl.show()
