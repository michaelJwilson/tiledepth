import os
import pylab as pl
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   astropy.table import Table


scr = os.environ['CSCRATCH']
dat = Table.read(scr + '/desi/tsnr/blanc/exptable.fits')

# dat.pprint()

bail = dat[dat['EXPID'] == 67972]
# bail.pprint()

fig, axes = plt.subplots(3, 1, figsize=(5, 10))

axes[0].hist(dat['LRGTSNR'], bins=100)
axes[1].hist(dat['ELGTSNR'], bins=100)
axes[2].hist(dat['QSOTSNR'], bins=100)

axes[0].set_title('LRG')
axes[1].set_title('ELG')
axes[2].set_title('QSO')

pl.show()
