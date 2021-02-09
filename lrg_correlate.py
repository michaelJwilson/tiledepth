import glob
import pylab as pl
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

from   astropy.table import Table
from   scipy         import stats


zhou  = Table.read('/global/cscratch1/sd/rongpu/desi/sv1/single_exp_coadd_blanc/lrg_redshift_efficiency.fits') 

# root = '/global/cscratch1/sd/mjwilson/trash/exposures/'
root = '/global/cscratch1/sd/mjwilson/desi/tsnr/blanc/cframes/exposures/'   

paths = glob.glob(root + '/*/*/cframe-r?-*.fits')

# ['fail_frac_lrg_sv', 'fail_frac_lrg_opt', 'fail_frac_lrg_ir', 'fail_n_lrg_sv', 'fail_n_lrg_opt', 'fail_n_lrg_ir']

fig,  axes  = plt.subplots(2, 3, figsize=(15, 10))
fig2, axes2 = plt.subplots(1, 1, figsize=(5, 5))

result = []

for i, path in enumerate(paths):
    try:
        dat   = fits.open(path)
        rtsnr = dat['SCORES'].data['LRGTSNR_R']
        
        cam   = path.split('/')[-1].split('-')[1]

        assert cam == 'r0'
    
        dat   = fits.open(path.replace('r0', 'b0'))
        btsnr = dat['SCORES'].data['LRGTSNR_B']

        dat   = fits.open(path.replace('r0', 'z0'))
        ztsnr = dat['SCORES'].data['LRGTSNR_Z']

        # Add up. tsnr^2.
        tsnr  = btsnr + rtsnr + ztsnr

    except:
        print('Missing [brz] for {} of {}'.format(i, len(paths)))
        continue
        
    expid = np.int(dat[0].header['EXPID'])
    night = dat[0].header['NIGHT']
    
    flavor = dat[0].header['FLAVOR']

    if flavor != 'science':
        continue

    zexp  = zhou[zhou['EXPID'] == expid]
    
    print('Solving for {} of {}'.format(i, len(paths)))

    if len(zexp) > 0:
        # In Rongpus LRG reductions. 
        axes[0,0].plot(np.median(btsnr), zexp['B_DEPTH'].data, marker='.', c='k', markersize=2)
        axes[0,1].plot(np.median(rtsnr), zexp['R_DEPTH'].data, marker='.', c='k', markersize=2)
        axes[0,2].plot(np.median(ztsnr), zexp['Z_DEPTH'].data, marker='.', c='k', markersize=2)

        axes[1,0].plot(np.median(tsnr), zexp['fail_frac_lrg_sv'].data,  marker='.', c='k', markersize=2)
        axes[1,1].plot(np.median(tsnr), zexp['fail_frac_lrg_opt'].data, marker='.', c='k', markersize=2)
        axes[1,2].plot(np.median(tsnr), zexp['fail_frac_lrg_ir'].data,  marker='.', c='k', markersize=2)

        axes2.plot(np.median(tsnr), zexp['deltachi2_ratio'].data, marker='.', c='k', markersize=2)

        result.append([np.median(tsnr), zexp['B_DEPTH'].data[0], zexp['R_DEPTH'].data[0], zexp['Z_DEPTH'].data[0], zexp['fail_frac_lrg_sv'].data[0],\
                       zexp['fail_frac_lrg_opt'].data[0], zexp['fail_frac_lrg_ir'].data[0], zexp['deltachi2_ratio'].data[0]])
        
    else:
        print('Unmatched: {:d}.'.format(expid))

result = np.array(result)
        
slope, intercept, r_value, p_value, std_err = stats.linregress(result[:,0], result[:,-1])
axes2.plot(result[:,0], slope * result[:,0] + intercept, c='k',lw=0.5)
axes2.text(0.1, 0.85, r'$y=$ {:.3f} $x$ + {:.3f} with $r=${:.3f}'.format(slope, intercept, r_value), transform=axes2.transAxes)

axes[0,0].set_xlabel(r'LRG <TSNR> $b$')
axes[0,1].set_xlabel(r'LRG <TSNR> $r$')
axes[0,2].set_xlabel(r'LRG <TSNR> $z$')

axes[0,0].set_ylabel(r'B DEPTH')
axes[0,1].set_ylabel(r'R DEPTH')
axes[0,2].set_ylabel(r'Z DEPTH')

axes[1,0].set_xlabel(r'LRG <TSNR>')
axes[1,1].set_xlabel(r'LRG <TSNR>')
axes[1,2].set_xlabel(r'LRG <TSNR>')

axes[1,0].set_ylabel(r'LRG SV: failure fraction')
axes[1,1].set_ylabel(r'LRG OPT: failure fraction')
axes[1,2].set_ylabel(r'LRG IR: failure fraction')

axes2.set_xlabel(r'LRG <TSNR>')
axes2.set_ylabel(r'$\Delta \chi^2$')

fig.savefig('plots/fig1.pdf')
fig2.savefig('plots/fig2.pdf')
