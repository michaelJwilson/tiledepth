import os
import glob
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from   astropy.table import Table, vstack, join


scr  = os.environ['CSCRATCH']
root = scr+'/desi/tsnr/blanc/'

cpath = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'
conds = Table.read(cpath)
conds = conds['EXPID', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']

petals = np.arange(10).astype(str)
bands = ['b', 'r', 'z']

result = None

for petal in petals:
    print(petal)

    # b{petal}
    bdat = Table.read(root + 'summary_{}.fits'.format('b' + petal))
    rdat = Table.read(root + 'summary_{}.fits'.format('r' + petal))
    zdat = Table.read(root + 'summary_{}.fits'.format('z' + petal))

    assert  np.all(bdat['EXPID'].data == rdat['EXPID'].data)
    assert  np.all(rdat['EXPID'].data == zdat['EXPID'].data)

    ssum = bdat['EXPID', 'NIGHT']
    ssum['PETAL'] = petal
    
    for t in ['ELGTSNR', 'BGSTSNR', 'QSOTSNR', 'LRGTSNR']:
        ssum[t]     = bdat[t] + rdat[t] + zdat[t]
        ssum[t+'B'] = bdat[t]
        ssum[t+'R'] = rdat[t]
        ssum[t+'Z'] = zdat[t]
        
    # ssum.pprint()

    if result == None:
        result = ssum

    else:
        result = vstack((result, ssum))

result_grouped = result.group_by('EXPID')
result_binned  = result_grouped.groups.aggregate(np.median)
result         = result_binned

result['EXPID'] = result['EXPID'].data.astype(np.int)

result = join(result, conds, keys='EXPID', join_type='left')
result.pprint()

result.write(root + '/exptable.fits', format='fits', overwrite=True)

data = Table(result, copy=True)    

data['BDEPTH'] = data['B_DEPTH']
data['RDEPTH'] = data['R_DEPTH']
data['ZDEPTH'] = data['Z_DEPTH']

del data['EXPID']
del data['B_DEPTH']
del data['R_DEPTH']
del data['Z_DEPTH']


names = data.dtype.names

fig, axes = plt.subplots(len(names), len(names), figsize=(40,40))
plt.subplots_adjust(hspace=1., wspace=1.)

for i, x in enumerate(names):
    for j, y in enumerate(names):
        if i >= j:
            axes[i,j].remove()
        
        axes[i,j].plot(data[x].data, data[y].data, marker='.', lw=0.0, markersize=1, c='k')
        axes[i,j].set_xlabel(x, fontsize=5)
        axes[i,j].set_ylabel(y, fontsize=5)

# plt.tight_layout()
pl.savefig('depth_corr.pdf', bbox_inches=0, pad_inches=0.0)
