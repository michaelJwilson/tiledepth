import os
import glob
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from   astropy.table import Table, vstack, join


scr  = os.environ['CSCRATCH']
root = scr+'/desi/tsnr/blanc/exptables/'

cpath = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'
conds = Table.read(cpath)
conds = conds['EXPID', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH', 'GFA_FWHM_ASEC', 'GFA_TRANSPARENCY', 'GFA_SKY_MAG_AB']

petals = np.arange(1).astype(str)
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

# Remove transparency ~0.2 which blows up ELGTSNRB.
result         = result[result['EXPID'] != 69416]

result = join(result, conds, keys='EXPID', join_type='left')

select=result[result['LRGTSNRB'] < 0.1]
select.sort('B_DEPTH')
select.pprint()

badlist = select['EXPID']
sconds = Table.read(cpath)
sconds=sconds[np.isin(sconds['EXPID'], badlist)]

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
        
        axes[i,j].scatter(data[x].data, data[y].data, c= data['GFA_SKY_MAG_AB'].data, marker='.', s=.5)
        axes[i,j].set_xlabel(x, fontsize=5)
        axes[i,j].set_ylabel(y, fontsize=5)

# plt.tight_layout()
pl.savefig('depth_corr.pdf', bbox_inches=0, pad_inches=0.0)
