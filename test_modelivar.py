from desispec.io import read_frame, write_frame
from desispec.io import read_fiberflat
from desispec.io import read_sky
from desispec.io import read_fibermap
from desispec.io import shorten_filename
from desispec.io.fluxcalibration import read_flux_calibration
from desispec.fiberflat import apply_fiberflat
from desispec.sky import subtract_sky
from desispec.fluxcalibration import apply_flux_calibration
from desiutil.log import get_logger
from desispec.cosmics import reject_cosmic_rays_1d
from desispec.specscore import compute_and_append_frame_scores
from desispec.fiberbitmasking import get_fiberbitmasked_frame
from model_var import calc_var
from scipy.optimize import minimize

import argparse
import sys
import copy
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

def parse(options=None):
    parser = argparse.ArgumentParser(description="Apply fiberflat, sky subtraction and calibration.")
    parser.add_argument('-i','--infile', type = str, default = None, required=True,
                        help = 'path of DESI exposure frame fits file')
    parser.add_argument('--fiberflat', type = str, default = None,
                        help = 'path of DESI fiberflat fits file')
    parser.add_argument('--sky', type = str, default = None,
                        help = 'path of DESI sky fits file')
    parser.add_argument('--calib', type = str, default = None,
                        help = 'path of DESI calibration fits file')
    parser.add_argument('--psf', type = str, default=None,
                        help = 'path of DESI calibration psf file (triggers tsnr)')
    parser.add_argument('--nea', type = str, default=None,
                        help = 'path of DESI measter nea file (required for tsnr)')
    args = None
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args):
    log = get_logger()

    frame=read_frame(args.infile, skip_resolution=True)
    fibermap=read_fibermap(args.infile)
    fiberflat=read_fiberflat(args.fiberflat)
    skymodel=read_sky(args.sky)
    fluxcalib=read_flux_calibration(args.calib)

    cam=args.infile.split('/')[-1].split('-')[1]
    band=cam[0]
    bands=[band]

    # Indices of sky fibers. 
    sky_indx = np.where(fibermap['OBJTYPE'] == 'SKY')[0]
    
    rd_var, sky_var = calc_var(bands, args.nea, args.psf, frame, fluxcalib, fiberflat, skymodel, components=True)   
    var = calc_var(bands, args.nea, args.psf, frame, fluxcalib, fiberflat, skymodel, components=False)
    
    nsky = 4 
    fig, axes = plt.subplots(1, nsky, figsize=(5 * nsky, 5))

    for i in range(nsky):
        def calc_alphavar(alpha):
            return alpha * rd_var[sky_indx,:] + sky_var[sky_indx,:]
        
        def alpha_fit(alpha):
            _var = calc_alphavar(alpha)
            ivar =  1. / _var
            X2 = (frame.ivar[sky_indx,:] - ivar)**2.
            return np.sum(X2)
        
        res = minimize(alpha_fit, x0=[1.])
        alpha = res.x[0]

        indx = sky_indx[i]
        
        axes[i].plot(skymodel.wave, frame.ivar[indx,:], lw=0.4, label='Sky frame IVAR', alpha=0.4)
        axes[i].plot(skymodel.wave, 1./rd_var[indx,:], lw=0.4, label='Model rd. IVAR', alpha=0.4)
        axes[i].plot(skymodel.wave, 1./sky_var[indx,:], lw=0.4, label='Model Sky IVAR', alpha=0.4)
        axes[i].plot(skymodel.wave, 1./var[indx,:], lw=0.4, label=r'Model IVAR', alpha=0.4)
        axes[i].plot(skymodel.wave, 1./calc_alphavar(alpha)[i,:], lw=0.4, label=r'$\alpha$ Model IVAR', alpha=0.4)
        axes[i].set_title(r'Fiber {:d} ($\alpha$ = {:.6f})'.format(indx, alpha))
        axes[i].set_xlabel(r'Wavelength [$AA$]')
        axes[i].set_yscale('log')
        axes[i].set_ylim(bottom=2.e-3, top=3.e-2)
        axes[i].legend(frameon=False, loc=2)
        
    axes[0].set_ylabel('e/A')
        
    pl.show()
    
    
if __name__ == '__main__':
    args = parse()
    
    main(args)
