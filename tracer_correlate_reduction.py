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
from   desitarget.targets import main_cmx_or_sv

prod     = 'cascades'
program  = 'ELG'              # ['ELG', 'QSO+LRG', 'BGS+MWS']
tsnrtype = 'ELG'
svtype   = 'ELG_FDR_GFIB'     # ['ELG_FDR_GTOT', 'ELG_FDR_GFIB']
                              # ['LRG_OPT', 'LRG_IR']
                              # ['QSO_COLOR_4PASS', 'QSO_COLOR_8PASS']

summary = Table.read('/project/projectdirs/desi/survey/observations/SV1/sv1-exposures.fits')
summary = summary[summary['TARGETS'] == program]
# summary.pprint()

expids = summary['EXPID'].data
nights = np.unique(summary['NIGHT'].data)

zbests = glob.glob('/global/cfs/cdirs/desi/spectro/redux/{}/tiles/*/exposures/zbest-*.fits'.format(prod))

tsnrs  = Table.read('/global/cscratch1/sd/mjwilson/trash/denali/tsnr-exposures-{}.fits'.format(prod), 'TSNR2_EXPID')
tsnrs  = tsnrs[tsnrs['SEEING_FWHM']  > 0.0] # Not missing.
tsnrs  = tsnrs[tsnrs['SEEING_FWHM'] != 1.1] # Not default nominal.
tsnrs.sort('NIGHT')
tsnrs.pprint()

# night/expid/
root   = '/global/cfs/cdirs/desi/spectro/redux/{}/exposures/'.format(prod)
paths  = []

for night in nights:
    nexps   = glob.glob(root + '/{}/*'.format(night))
    nexps   = np.array([x.split('/')[-1].replace('.fits','') for x in nexps])
    nexps   = np.array([x for x in nexps if x[0] == '0']).astype(np.int)
    
    keep    = np.isin(nexps, expids)
    nexps   = nexps[keep]
    
    for x in nexps:
        npaths  = glob.glob(root + '/{}/{:08d}/cframe-r0-*'.format(night, x))
        paths  += npaths
        
fig,  axes  = plt.subplots(3, 4, figsize=(20, 30))
plt.subplots_adjust(hspace=0.3)

result = []

for i, path in enumerate(paths):
    try:
        dat   = fits.open(path)
        '''
        rtsnr = dat['SCORES'].data['TSNR2_{}_R'.format(tsnrtype)]
        
        cam   = path.split('/')[-1].split('-')[1]
        petal = cam[1]
    
        dat   = fits.open(path.replace('r0', 'b0'))
        btsnr = dat['SCORES'].data['TSNR2_{}_B'.format(tsnrtype)]
    
        dat   = fits.open(path.replace('r0', 'z0'))
        ztsnr = dat['SCORES'].data['TSNR2_{}_Z'.format(tsnrtype)]

        # Add up. tsnr^2.
        tsnr  = btsnr + rtsnr + ztsnr
        '''
    except Exception as e:
        print('Failed on {}'.format(path))
        print(e)

        continue
    
    expid  = np.int(dat[0].header['EXPID'])

    night  = dat[0].header['NIGHT']
    flavor = dat[0].header['FLAVOR']
    tileid = dat[0].header['TILEID']
    prog   = dat[0].header['PROGRAM']

    if expid not in tsnrs['EXPID'].data:
        print('Warning: no new tsnr reduction for {}'.format(expid))
        continue

    tsnr   = tsnrs[tsnrs['EXPID'] == expid]['TSNR2_{}'.format(tsnrtype)]

    btsnr = rtsnr = ztsnr = tsnr / 3.
    
    try:
        zbest  = Table.read('/global/cfs/cdirs/desi/spectro/redux/{}/tiles/{}/exposures/zbest-0-{}-{:08d}.fits'.format(prod, tileid, tileid, expid), 'ZBEST')
        fmap   = Table.read('/global/cfs/cdirs/desi/spectro/redux/{}/tiles/{}/exposures/zbest-0-{}-{:08d}.fits'.format(prod, tileid, tileid, expid), 'FIBERMAP') 

    except:
        print('Failed to find: {}'.format('/global/cfs/cdirs/desi/spectro/redux/{}/tiles/{}/exposures/zbest-0-{}-{:08d}.fits'.format(prod, tileid, tileid, expid)))
        continue
        
    if svtype.split('_')[0] == 'BGS':
         isin = (fmap['SV1_DESI_TARGET'] & sv1_desi_mask['BGS_ANY']) != 0
         isin = isin & ((fmap['SV1_BGS_TARGET'] & sv1_bgs_mask[svtype]) != 0)
        
    else:
        isin = (fmap['SV1_DESI_TARGET'] & sv1_desi_mask[svtype]) != 0

    isin   =  isin & (fmap['FIBERSTATUS'] == 0)
    
    tids   = fmap['TARGETID'].data[isin]
    
    goodz  = (zbest['ZWARN'] == 0) & (zbest['DELTACHI2'] > 25.) & np.isin(zbest['TARGETID'], tids) 
    
    zids   = zbest['TARGETID'].data[goodz]

    if len(tids) == 0:
        continue
    
    zeff   = 100. * len(zids) / len(tids)

    print(expid, len(tids), len(zids), 100. * len(zids) / len(tids))

    sample = summary[summary['EXPID'] == np.int(expid)]
    
    fwhm   = sample['GFA_FWHM_ASEC']
    ebv    = np.log10(np.median(fmap['EBV']))
    
    # color  = fwhm
    # vmin   =  1.0
    # vmax   =  2.0

    color  =  ebv
    vmin   =   -3.
    vmax   =   -1.

    # axes[2,0].scatter(np.median(btsnr), zeff, marker='.', c=color, s=2, vmin=vmin, vmax=vmax)
    # axes[2,1].scatter(np.median(rtsnr), zeff, marker='.', c=color, s=2, vmin=vmin, vmax=vmax)
    # axes[2,2].scatter(np.median(ztsnr), zeff, marker='.', c=color, s=2, vmin=vmin, vmax=vmax)
    axes[2,3].scatter(np.median(tsnr),  zeff, marker='.', c=color, s=2, vmin=vmin, vmax=vmax)

    dark_depth   = sample['EFFTIME_DARK']
    bright_depth = sample['EFFTIME_BRIGHT']

    axes[1,0].plot(dark_depth, zeff, marker='.', c='k', markersize=2)
    axes[1,1].plot(bright_depth, zeff, marker='.', c='k', markersize=2)

    bdepth = sample['B_DEPTH_EBVAIR']
    rdepth = sample['R_DEPTH_EBVAIR']
    zdepth = sample['Z_DEPTH_EBVAIR']
    
    axes[0,0].plot(bdepth, zeff, marker='.', c='k', markersize=2)
    axes[0,1].plot(rdepth, zeff, marker='.', c='k', markersize=2)
    axes[0,2].plot(zdepth, zeff, marker='.', c='k', markersize=2)

    for i in range(4):
        axes[0,i].set_ylim(0., 75.)
        axes[1,i].set_ylim(0., 75.)
        axes[2,i].set_ylim(0., 75.)
        
        axes[0,i].set_ylabel('{} z success [%]'.format(tsnrtype), fontsize=9)
        axes[1,i].set_ylabel('{} z success [%]'.format(tsnrtype), fontsize=9)
        axes[2,i].set_ylabel('{} z success [%]'.format(tsnrtype), fontsize=9)
        
    axes[0,0].set_xlabel('B DEPTH_EBVAIR')
    axes[0,1].set_xlabel('R DEPTH_EBVAIR')
    axes[0,2].set_xlabel('Z DEPTH_EBVAIR')
    
    axes[2,0].set_xlabel('{} TSNR$^2$ $b$'.format(tsnrtype))
    axes[2,1].set_xlabel('{} TSNR$^2$ $r$'.format(tsnrtype))
    axes[2,2].set_xlabel('{} TSNR$^2$ $z$'.format(tsnrtype))
    axes[2,3].set_xlabel('{} TSNR$^2$ tot'.format(tsnrtype))

    axes[1,0].set_xlabel('EFFTIME DARK')
    axes[1,1].set_xlabel('EFFTIME BRIGHT')
    
    result.append([np.median(btsnr), np.median(rtsnr), np.median(ztsnr), np.median(tsnr), bdepth, rdepth, zdepth, dark_depth, bright_depth, zeff])

    
result = np.array(result).astype(np.float)

btsnr  = result[:,0]
rtsnr  = result[:,1]
ztsnr  = result[:,2]
tsnr   = result[:,3]

bdep   = result[:,4] 
rdep   = result[:,5]
zdep   = result[:,6]

dark   = result[:,7]
bright = result[:,8]

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
   
    axes[0,i].plot(np.arange(0., 10. * consta, 2.), model(res.x, np.arange(0., 10. * consta, 2.)), c='k', lw=0.25)
    axes[0,i].set_title('{:.3f} (1. - exp(x / {:.3f})) r.m.s. {:.3f}'.format(constb, consta, np.std(data - model(res.x, xvals))), fontsize=9.)
    axes[0,i].set_xlim(-1., 10. * consta)
    axes[0,i].set_ylim(0., 105.)

for i, xvals in enumerate([dark, bright]):
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

    axes[1,i].plot(np.arange(0., 10. * consta, 2.), model(res.x, np.arange(0., 10. * consta, 2.)), c='k', lw=0.25)
    axes[1,i].set_title('{:.3f} (1. - exp(x / {:.3f})) r.m.s. {:.3f}'.format(constb, consta, np.std(data - model(res.x, xvals))), fontsize=9.)
    axes[1,i].set_xlim(-1., 10. * consta)
    axes[1,i].set_ylim(0., 105.)

for i, (xvals, x0) in enumerate(zip([btsnr, rtsnr, ztsnr, tsnr], [0.5, 40., 100., 100.])):
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
    
    axes[2,i].plot(np.arange(0., 10. * consta, 2.), model(res.x, np.arange(0., 10. * consta, 2.)), c='k', lw=0.25)
    axes[2,i].set_title('{:.3f} (1. - exp(x / {:.3f})) r.m.s. {:.3f}'.format(constb, consta, np.std(data - model(res.x, xvals))), fontsize=9.)

    axes[1,i].set_xlim(left=-5.)
    axes[2,i].set_xlim(-1., 10. * consta)

    axes[1,i].set_ylim(0., 105.)
    axes[2,i].set_ylim(0., 105.)

fig.suptitle('{} singles {}: {}'.format(prod, program, svtype))
pl.show()
