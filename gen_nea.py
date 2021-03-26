import os
import glob
import numpy as np
import pylab as pl
import itertools
import matplotlib.pyplot as plt
import astropy.io.fits as fits

size   = 'x-small'
params = {'axes.labelsize': size, 'axes.titlesize': size, 'xtick.labelsize': size, 'ytick.labelsize': size}

pl.rcParams.update(params)

'''
scr = os.environ['CSCRATCH']
dmodel = os.environ['DESIMODEL']

for band in ['b', 'r', 'z']:
    for petal in np.arange(10).astype(str):
        cam = band + petal
        cmd = 'python desispec/scripts/master_nea.py -i /global/homes/m/mjwilson/blanc/exposures/20201223/00069580/psf-{}-00069580.fits --outdir {}/desi/nea/'.format(cam, scr)

        print(cmd)
        
        # os.system(cmd)

        exit(0)
'''        
fig, axes = plt.subplots(5, 6, figsize=(20,20))

# root = '{}/desi/nea/'.format(scr)
# root = dmodel + '/data/specpsf/nea/'
root = '/global/homes/m/mjwilson/desimodel/trunk/data/specpsf/nea/'

col  = 0

petals = np.arange(10) 
cameras = [x[0] + x[1] for x in itertools.product(['b', 'r', 'z'], petals.astype(np.str))]

for i, cam in enumerate(cameras):
    nea = root + '/masternea_{}.fits'.format(cam)
    
    row = i % 5

    cam = nea.split('_')[-1].replace('.fits', '')
      
    dat = fits.open(nea)
    wave = dat['WAVELENGTH'].data
    nea = dat['NEA'].data

    axes[row, col].imshow(nea, aspect='auto', vmin=3.3, vmax=4.5, extent=(wave.min(), wave.max(), 0, 499))
    axes[row, col].set_title(cam, y=0.8, c='white')
    
    print(cam, nea.min(), nea.max())
    
    if row == 4:
        col += 1
        
plt.tight_layout()
pl.show()

