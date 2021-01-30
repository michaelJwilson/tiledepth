import os
import glob
import astropy.io.fits as fits
import numpy as np

from pathlib import Path
from   desispec.io.meta import findfile, specprod_root
from   desispec.calibfinder import CalibFinder
from   desispec.io import read_frame

targets = 'ELG' # 'QSO+LRG'
tile = '80608'

scr = os.environ['CSCRATCH']
blanc = '/global/cfs/cdirs/desi/spectro/redux/blanc/'
cframes = list(glob.glob('/global/cfs/cdirs/desi/spectro/redux/blanc/exposures/*/*/cframe-r0-*.fits'))
calibdir = '/global/cfs/cdirs/desi/spectro/desi_spectro_calib/trunk/'

os.environ['DESI_SPECTRO_CALIB'] = '/global/cfs/cdirs/desi/spectro/desi_spectro_calib/trunk/'
os.environ['SPECPROD'] = 'blanc'

for x in cframes:
    hdul = fits.open(x)
    hdr  = hdul[0].header

    flavor = hdr['FLAVOR']
    prog = hdr['PROGRAM']

    if flavor == 'science':
        parts = prog.split(' ')

        if parts[0] == 'SV1':
            if parts[1] == targets:        
                parts  = x.split('/')
                
                night  = parts[9]
                expid  = np.int(parts[10])

                name   = parts[-1]

                parts  = name.split('-')
                camera = parts[1] 

                calib  = findfile('fluxcalib', night=night, expid=expid, camera=camera, specprod_dir=None)
                
                # findfile('cframe', night=night, expid=expid, camera=camera, specprod_dir=blanc)
                cframe  = fits.open(x)
                hdr    = cframe[0].header['FIBERFLT']
                hdr    = hdr.replace('SPECPROD', blanc)

                tileid = cframe[0].header['TILEID']

                if tile != tile:
                    continue
                                    
                # print(night, expid, camera, calib, hdr)

                iin = x.replace('cframe', 'frame')
                sky = x.replace('cframe', 'sky')
                psf = sky.replace('sky', 'psf')
                nea = '/project/projectdirs/desi/users/mjwilson/master_nea/masternea_{}.fits'.format(camera)  
                ens = '/project/projectdirs/desi/users/mjwilson/tsnr-ensemble/'
                out = x.replace(blanc, scr + '/trash/')

                Path(os.path.dirname(out)).mkdir(parents=True, exist_ok=True)
                
                # desi_process_exposure --infile $INFILE -o $OUTFILE --calib $CALIB --sky $SKY --fiberflat $FIBERFLAT --psf $PSF --nea $NEA --ensembledir $ENSEMBLE
                cmd = 'desi_process_exposure --infile {} -o {} --calib {} --sky {} --fiberflat {} --psf {} --nea {} --ensembledir {}'.format(iin, out, calib, sky, hdr, psf, nea, ens)

                os.system(cmd)
