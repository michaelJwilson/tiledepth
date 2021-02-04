import numpy as np
import astropy.io.fits as fits
import glob
import numpy as np

from desispec.io.spectra import Spectra
from astropy.convolution import convolve, Box1DKernel
from scipy.interpolate import RectBivariateSpline
from specter.psf.gausshermite  import  GaussHermitePSF
from desispec.specscore import append_frame_scores
from desispec.tsnr import read_nea, var_model, fb_rdnoise
        
def calc_var(bands, neadir, psfpath, frame, fluxcalib, fiberflat, skymodel, components=False):
    '''
    Stripped down tsnr. 
    '''

    assert len(bands) == 1

    psf=GaussHermitePSF(psfpath)
    
    nea, angperpix=read_nea(neadir)

    nspec, nwave = fluxcalib.calib.shape
    
    fibers = np.arange(nspec)
    rdnoise = fb_rdnoise(fibers, frame, psf)

    npix = nea(fibers, frame.wave)
    angperpix = angperpix(fibers, frame.wave)
    
    for label, x in zip(['RDNOISE', 'NEA', 'ANGPERPIX'], [rdnoise, npix, angperpix]):
        print('{} \t {:.3f} +- {:.3f}'.format(label.ljust(10), np.median(x), np.std(x)))

    for band in bands:                    
        denom = var_model(rdnoise, npix, angperpix, 0.8, fiberflat, skymodel, components=components)

    return  denom
