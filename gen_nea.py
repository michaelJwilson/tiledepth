import os
import glob
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt


scr = os.environ['CSCRATCH']

'''
for band in ['b', 'r', 'z']:
    for petal in np.arange(10).astype(str):
        cam = band + petal
        cmd = 'python desispec/scripts/master_nea.py -i /global/homes/m/mjwilson/blanc/exposures/20201223/00069580/psf-{}-00069580.fits --outdir {}/desi/nea/'.format(cam, scr)

        os.system(cmd)
'''

neas = glob.glob('{}/desi/nea/*'.format(scr))

fig, axes = plt.subplots(5, 6, figsize=(10,10))

pl.show()

