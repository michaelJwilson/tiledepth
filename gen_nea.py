import os
import numpy as np


scr = os.environ['CSCRATCH']

for band in ['b', 'r', 'z']:
    for petal in np.arange(10).astype(str):
        cam = band + petal
        cmd = 'python desispec/scripts/master_nea.py -i /global/homes/m/mjwilson/blanc/exposures/20201223/00069580/psf-{}-00069580.fits --outdir {}/trash/'.format(cam, scr)

        os.system(cmd)

