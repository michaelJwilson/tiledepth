import os
import glob
import pylab as pl
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   astropy.table import Table, join
from   scipy import stats

from   desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
from   desitarget.sv1.sv1_targetmask import bgs_mask as sv1_bgs_mask
from   scipy.optimize import minimize


tracer = 'ELG' # ['ELG', 'LRG+QSO', 'BGS+MWS']
subset = 'ELG'
svtype = 'ELG_SV_GTOT'

zbests = glob.glob('/global/cfs/cdirs/desi/spectro/redux/blanc/tiles/*/exposures/zbest-*.fits')
expids = [os.path.basename(x).split('-')[-1].replace('.fits', '') for x in zbests]
expids = [np.int(x) for x in expids if x[0:2] == '00']

# root = '/global/cscratch1/sd/mjwilson/trash/exposures/'
root   = '/global/cscratch1/sd/mjwilson/desi/tsnr/blanc/cframes/exposures/'
paths  = glob.glob(root + '/*/*/cframe-r?-*.fits')

conds  = Table.read('sv1-exposures.fits')
# print(conds.dtype.names)

fig,  axes  = plt.subplots(2, 4, figsize=(20, 20))

result = []

for i, path in enumerate(paths):
    try:
        # Missing arm.
        dat   = fits.open(path)
        # tids  = dat['FIBERMAP'].data['TARGETID'] 
        
        rtsnr = dat['SCORES'].data['{}TSNR_R'.format(subset)]
        
        cam   = path.split('/')[-1].split('-')[1]
        petal = cam[1]
    
        dat   = fits.open(path.replace('r0', 'b0'))
        btsnr = dat['SCORES'].data['{}TSNR_B'.format(subset)]
    
        dat   = fits.open(path.replace('r0', 'z0'))
        ztsnr = dat['SCORES'].data['{}TSNR_Z'.format(subset)]

        # Add up. tsnr^2.
        tsnr  = btsnr + rtsnr + ztsnr

    except:
        continue
    
    expid  = np.int(dat[0].header['EXPID'])

    if expid not in expids:
        continue

    night  = dat[0].header['NIGHT']
    flavor = dat[0].header['FLAVOR']
    tileid = dat[0].header['TILEID']
    prog   = dat[0].header['PROGRAM']
    
    if flavor != 'science':
        continue

    parts  = prog.split(' ')

    if parts[1] != tracer:
        print('Skipping {}'.format(parts[1]))
        continue

    else:
        print('Reducing {}'.format(parts[1]))
        
    zbest  = Table.read('/global/cfs/cdirs/desi/spectro/redux/blanc/tiles/{}/exposures/zbest-0-{}-{:08d}.fits'.format(tileid, tileid, expid), 'ZBEST')
    fmap   = Table.read('/global/cfs/cdirs/desi/spectro/redux/blanc/tiles/{}/exposures/zbest-0-{}-{:08d}.fits'.format(tileid, tileid, expid), 'FIBERMAP') 
    
    isin   = (fmap['SV1_DESI_TARGET'] & sv1_desi_mask[svtype]) != 0
    isin   =  isin & (fmap['FIBERSTATUS'] == 0)
    
    tids   = fmap['TARGETID'].data[isin]
    
    goodz  = (zbest['ZWARN'] == 0) & (zbest['DELTACHI2'] > 25.) & np.isin(zbest['TARGETID'], tids) 
    
    zids   = zbest['TARGETID'].data[goodz]

    print(expid, len(tids), len(zids))

    if len(tids) == 0:
        continue
    
    zeff   = 100. * len(zids) / len(tids)

    # print(len(tids), len(zids), 100. * len(zids) / len(tids))
    
    axes[1,0].plot(np.median(btsnr), zeff, marker='.', c='k', markersize=2)
    axes[1,1].plot(np.median(rtsnr), zeff, marker='.', c='k', markersize=2)
    axes[1,2].plot(np.median(ztsnr), zeff, marker='.', c='k', markersize=2)
    axes[1,3].plot(np.median(tsnr),  zeff, marker='.', c='k', markersize=2)

    sample = conds[conds['EXPID'] == expid]

    bdepth = sample['B_DEPTH_EBVAIR']
    rdepth = sample['R_DEPTH_EBVAIR']
    zdepth = sample['Z_DEPTH_EBVAIR']
    
    axes[0,0].plot(bdepth, zeff, marker='.', c='k', markersize=2)
    axes[0,1].plot(rdepth, zeff, marker='.', c='k', markersize=2)
    axes[0,2].plot(zdepth, zeff, marker='.', c='k', markersize=2)

    for i in range(4):
        axes[0,i].set_ylim(0., 75.)
        axes[1,i].set_ylim(0., 75.)

        axes[0,i].set_ylabel('{} z success [%]'.format(subset), fontsize=9)
        axes[1,i].set_ylabel('{} z success [%]'.format(subset), fontsize=9)
        
    axes[0,0].set_xlabel('B DEPTH_EBVAIR')
    axes[0,1].set_xlabel('R DEPTH_EBVAIR')
    axes[0,2].set_xlabel('Z DEPTH_EBVAIR')
        
    axes[1,0].set_xlabel('{} TSNR $b$'.format(subset))
    axes[1,1].set_xlabel('{} TSNR $r$'.format(subset))
    axes[1,2].set_xlabel('{} TSNR $z$'.format(subset))
    axes[1,3].set_xlabel('{} TSNR tot'.format(subset))

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

    def model(x, xvals):
        a = x[0]
        b = x[1]
        
        _scaled = -xvals / a
        return  b * (1. - np.exp(_scaled))

    def X2(x, xvals):        
        _model  = model(x, xvals)
        return  np.sum((data - _model)**2.)

    res    = minimize(X2, x0=[500., 90.], args=(xvals))
    consta = res.x[0]
    constb = res.x[1]
    
    axes[0,i].plot(np.sort(result[:,4+i]), model(res.x, np.sort(result[:,4+i])), c='k', lw=0.25)
    axes[0,i].set_title('{:.3f} (1. - exp(x / {:.3f})) r.m.s. {:.3f}'.format(constb, consta, np.std(data - model(res.x, xvals))), fontsize=9.)
    axes[0,i].set_ylim(0., 105.)

for i, (xvals, x0) in enumerate(zip([btsnr, rtsnr, ztsnr, tsnr], [0.5, 40., 400., 400.])):
    indx = np.argsort(xvals)

    xvals = xvals[indx]
    data = result[indx,-1]
    
    def model(x, xvals):
        a = x[0]
        b = x[1]
        
        _scaled = -xvals / a

        return  b * (1. - np.exp(_scaled))

    def X2(x, xvals):        
        _model  = model(x, xvals)
        return  np.sum((data - _model)**2.)

    res   = minimize(X2, x0=[x0, 90.], args=(xvals))
    consta = res.x[0]
    constab = res.x[1]
    
    axes[1,i].plot(np.sort(result[:,i]), model(res.x, np.sort(result[:,i])), c='k', lw=0.25)
    axes[1,i].set_title('{:.3f} (1. - exp(x / {:.3f})) r.m.s. {:.3f}'.format(constb, consta, np.std(data - model(res.x, xvals))), fontsize=9.)
    axes[1,i].set_ylim(0., 105.)

fig.suptitle('Blanc singles {}: {}'.format(tracer, svtype))
pl.show()
