import argparse
import sys
import copy

from desispec.scripts.procexp import main


# export DESI_SPECTRO_CALIB=/global/cfs/cdirs/desi/spectro/desi_spectro_calib/trunk/  
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
    parser.add_argument('-o','--outfile', type = str, default = None, required=True,
                        help = 'path of output fits file')
    parser.add_argument('--cosmics-nsig', type = float, default = 0, required=False,
                        help = 'n sigma rejection for cosmics in 1D (default, no rejection)')
    parser.add_argument('--no-sky-throughput-correction', action='store_true',
                        help = 'Do NOT apply a throughput correction when subtraction the sky')
    parser.add_argument('--no-zero-ivar', action='store_true',
                        help = 'Do NOT set ivar=0 for masked pixels')
    parser.add_argument('--no-tsnr', action='store_true',
                        help = 'Do not compute template SNR')

    args = None
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

args = parse()

main(args)

# python bad_alpha.py --infile /global/cfs/cdirs/desi/spectro/redux/cascades/exposures/20210224/00077950/frame-z7-00077950.fits --calib /global/cfs/cdirs/desi/spectro/redux/cascades/exposures/20210224/00077950/fluxcalib-z7-00077950.fits --sky /global/cfs/cdirs/desi/spectro/redux/cascades/exposures/20210224/00077950/sky-z7-00077950.fits --fiberflat /global/cfs/cdirs/desi/spectro/redux/cascades/calibnight/20210224/fiberflatnight-z7-20210224.fits --cosmics-nsig 6 --outfile $CSCRATCH/trash/blah.fits
