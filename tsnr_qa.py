import os
import pylab as pl
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   astropy.table import Table


scr = os.environ['CSCRATCH']

alpha=1.0
fig, axes = plt.subplots(3, 2, figsize=(10, 15))
plt.subplots_adjust(hspace=0.3)

for i, camera in enumerate(['b0', 'r0', 'z0']):
    dat = Table.read(scr + '/desi/tsnr/blanc/exptables/summary_{}.fits'.format(camera))
    dat = dat[dat['EXPID'] != '00069416']
    dat.pprint()

    # bail = dat[dat['EXPID'] == 67972]
    # bail.pprint()

    # ALPHA      ELGTSNR     BGSTSNR     QSOTSNR     LRGTSNR
    # axes[i,0].hist(dat['BGSTSNR'], bins=20, label='BGS', alpha=alpha, histtype='step')
    axes[i,0].hist(dat['LRGTSNR'], bins=20, label='LRG', alpha=alpha, histtype='step')
    axes[i,0].hist(dat['ELGTSNR'], bins=20, label='ELG', alpha=alpha, histtype='step')
    axes[i,0].hist(dat['QSOTSNR'], bins=20, label='QSO', alpha=alpha, histtype='step')
    
    axes[i,1].hist(dat['ALPHA'], bins=20, alpha=alpha, histtype='step')

    axes[i,0].set_xlabel(r'TSNR$^2$')
    axes[i,0].set_ylabel(camera)

    axes[i,1].set_xlabel(r'$\alpha$')
    
    axes[i,0].legend()
    
    axes[i,0].set_yscale('log')
    axes[i,1].set_yscale('log')
    
pl.show()
