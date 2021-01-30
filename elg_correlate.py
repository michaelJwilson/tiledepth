import glob
import pylab as pl
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   astropy.table import Table
from   scipy         import stats

from   desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
from   desitarget.sv1.sv1_targetmask import bgs_mask as sv1_bgs_mask
from   scipy.optimize import minimize


zbests = glob.glob('/global/cscratch1/sd/raichoor/spectro/redux/blanc/tiles/80608/*/zbest-?-80608-*.fits')
expids = [x.split('/')[-1].split('-')[3].replace('.fits', '') for x in zbests]

expids = [np.int(x) for x in expids if x[0] == '0']

paths  = glob.glob('/global/cscratch1/sd/mjwilson/trash/exposures/*/*/cframe-r?-*.fits')

conds  = Table.read('sv1-exposures.fits')

fig,  axes  = plt.subplots(2, 4, figsize=(20, 20))

result = []

for i, path in enumerate(paths):
    try:
        dat   = fits.open(path)
        rtsnr = dat['SCORES'].data['LRGTSNR_R']
        
        cam   = path.split('/')[-1].split('-')[1]
        petal = cam[1]
    
        dat   = fits.open(path.replace('r0', 'b0'))
        btsnr = dat['SCORES'].data['LRGTSNR_B']

        dat   = fits.open(path.replace('r0', 'z0'))
        ztsnr = dat['SCORES'].data['LRGTSNR_Z']

        # Add up. tsnr^2.
        tsnr  = btsnr + rtsnr + ztsnr

    except:
        continue
        
    expid = np.int(dat[0].header['EXPID'])

    if expid not in expids:
        continue

    night = dat[0].header['NIGHT']
    flavor = dat[0].header['FLAVOR']

    if flavor != 'science':
        continue

    zbest  = Table.read('/global/cscratch1/sd/raichoor/spectro/redux/blanc/tiles/80608/{:08d}/zbest-{}-80608-{:08d}.fits'.format(expid, petal, expid), 'ZBEST')
    fmap   = Table.read('/global/cscratch1/sd/raichoor/spectro/redux/blanc/tiles/80608/{:08d}/zbest-{}-80608-{:08d}.fits'.format(expid, petal, expid), 'FIBERMAP') 

    # print('Loaded {}'.format(expid))

    elg    = (fmap['SV1_DESI_TARGET'] & sv1_desi_mask['ELG_FDR_GFIB']) != 0
    elg    =  elg | (fmap['SV1_DESI_TARGET'] & sv1_desi_mask['ELG_FDR_GTOT']) != 0

    # elg  = (fmap['SV1_DESI_TARGET'] & sv1_desi_mask['ELG']) != 0
    
    elg    =  elg & (fmap['FIBERSTATUS'] == 0)
    
    tids   = fmap['TARGETID'].data[elg]
    
    goodz  = (zbest['ZWARN'] == 0) & (zbest['DELTACHI2'] > 25.) & np.isin(zbest['TARGETID'], tids) 
    
    zids   = zbest['TARGETID'].data[goodz]

    zeff   = 100. * len(zids) / len(tids)
    
    # print(len(tids), len(zids), 100. * len(zids) / len(tids))
    
    axes[1,0].plot(np.median(btsnr), zeff, marker='.', c='k', markersize=2)
    axes[1,1].plot(np.median(rtsnr), zeff, marker='.', c='k', markersize=2)
    axes[1,2].plot(np.median(ztsnr), zeff, marker='.', c='k', markersize=2)
    axes[1,3].plot(np.median(tsnr),  zeff, marker='.', c='k', markersize=2)

    sample = conds[conds['EXPID'] == expid]

    bdepth = sample['B_DEPTH']
    rdepth = sample['R_DEPTH']
    zdepth = sample['Z_DEPTH']
    
    axes[0,0].plot(bdepth, zeff, marker='.', c='k', markersize=2)
    axes[0,1].plot(rdepth, zeff, marker='.', c='k', markersize=2)
    axes[0,2].plot(zdepth, zeff, marker='.', c='k', markersize=2)
    
    for i in range(4):
        axes[0,i].set_ylim(0., 75.)
        axes[1,i].set_ylim(0., 75.)

        axes[0,i].set_ylabel('ELG z success [%]', fontsize=9)
        axes[1,i].set_ylabel('ELG z success [%]', fontsize=9)

    axes[0,0].set_xlabel('B DEPTH')
    axes[0,1].set_xlabel('R DEPTH')
    axes[0,2].set_xlabel('Z DEPTH')
        
    axes[1,0].set_xlabel('ELG TSNR $b$')
    axes[1,1].set_xlabel('ELG TSNR $r$')
    axes[1,2].set_xlabel('ELG TSNR $z$')
    axes[1,3].set_xlabel('ELG TSNR tot')

    result.append([np.median(btsnr), np.median(rtsnr), np.median(ztsnr), np.median(tsnr), bdepth, rdepth, zdepth, zeff])

result = np.array(result).astype(np.float)

btsnr = result[:,0]
rtsnr =	result[:,1]
ztsnr =	result[:,2]
tsnr  =	result[:,3]

bdep  = result[:,4] 
rdep  = result[:,5]
zdep  = result[:,6]

for i, xvals in enumerate([bdep, rdep, zdep]):
    indx = np.argsort(xvals)

    xvals = xvals[indx]
    data = result[indx,-1]

    def model(a, xvals):
        _scaled = -xvals / a
        return  65. * (1. - np.exp(_scaled))

    def X2(a, xvals):
        _model  = model(a, xvals)
        return  np.sum((data - _model)**2.)

    res   = minimize(X2, x0=[500.], args=(xvals))
    const = res.x[0]

    axes[0,i].plot(xvals, model(const, xvals), c='k', lw=0.25)
    axes[0,i].set_title('1. - exp(x / {:.3f}) r.m.s. {:.3f}'.format(const, np.std(data - model(const, xvals))), fontsize=9.)

for i, xvals in enumerate([btsnr, rtsnr, ztsnr, tsnr]):
    indx = np.argsort(xvals)

    xvals = xvals[indx]
    data = result[indx,-1]
    
    def model(a, xvals):
        _scaled = -xvals / a

        return  65. * (1. - np.exp(_scaled))

    def X2(a, xvals):
        _model  = model(a, xvals)
        return  np.sum((data - _model)**2.)

    res   = minimize(X2, x0=[.5], args=(xvals))
    const = res.x[0]

    axes[1,i].plot(xvals, model(const, xvals), c='k', lw=0.25)
    axes[1,i].set_title('1. - exp(x / {:.3f}) r.m.s. {:.3f}'.format(const, np.std(data - model(const, xvals))), fontsize=9.)

fig.suptitle('Lynx singles: Blanc')
    
pl.show()
# pl.savefig('elgs.pdf')
