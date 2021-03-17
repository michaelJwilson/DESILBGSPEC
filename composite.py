import os
import sys
import fitsio
import redrock.templates
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

from   astropy.convolution import convolve, Gaussian1DKernel
from   astropy.table import Table, join
from   desitarget.sv1.sv1_targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask


dchi2     = 25.
tile      = 80869

nlbg      = 0
compwave  = np.arange(1000., 5000., 1.)
comp      = []

for petal in [0,1,2,3,4,5,7,8,9]:
    for band in ['B']:
        dat    = Table.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/zbest-{:d}-{:d}-deep.fits'.format(tile, petal, tile))
        dat    = dat[dat['DELTACHI2'] > 25.]
    
        fmap   = fitsio.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/coadd-{:d}-{:d}-deep.fits'.format(tile, petal, tile), 'FIBERMAP')
        wave   = fitsio.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/coadd-{:d}-{:d}-deep.fits'.format(tile, petal, tile), '{}_WAVELENGTH'.format(band))
        flux   = fitsio.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/coadd-{:d}-{:d}-deep.fits'.format(tile, petal, tile), '{}_FLUX'.format(band))

        isin   = (fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_MAIN']) != 0
        isin  |= (fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_HP'])   != 0	

        isin  &= (fmap['FIBERSTATUS'] == 0)

        flux   = flux[isin]
        fmap   = fmap[isin]

        fmap   = join(fmap, dat, keys='TARGETID', join_type='left')

        # Cut on dchisq.
        flux   = flux[~fmap['Z'].mask]
        fmap   = fmap[~fmap['Z'].mask]

        # Cut on dchisq.
        flux   = flux[fmap['Z'] > 2.]
        fmap   = fmap[fmap['Z'] > 2.]

        rwave  = np.ones_like(fmap['Z'])[:,None] * wave[None,:] 

        for i in range(len(fmap)):
            rwave[i,:] /= (1. + fmap['Z'][i])

        mask   = (1250. < rwave) & (rwave < 1500.)
        mask   =  mask.astype(np.float)

        mflux  =  mask * flux
        
        for i in range(len(fmap)):
            row  = mflux[i,:]
            row  = row[row > 0.0]

            norm = np.median(row)

            flux[i] /= norm 

            idx   = np.digitize(rwave[i,:], bins=compwave)
            toadd = np.zeros_like(compwave)
            toadd[idx] = flux[i,:]
            
            comp.append(toadd)

        nlbg += len(flux)
                       
    print(petal, len(flux))

comp  = np.array(comp)    
stds  = np.std(comp, axis=0) 

for i in range(len(compwave)):
    if stds[i] > 0.0:
        std = stds[i]

        outlier = comp[:,i] > 10. * std

        comp[outlier,i] = 0.0

comp  = np.sum(comp, axis=0) 
        
good  = stds > 0.0 

comp /= nlbg

for ab in [1260.4221, 1302.1685, 1304.3702, 1334.5323]:
    pl.axvline(ab, c='b', lw=0.1)
    
for em in [1264.738, 1309.276, 1533.431]:
    pl.axvline(em, c='k', lw=0.1)

pl.plot(compwave[good], comp[good], lw=0.2, c='k')
pl.title('{:d} LBG composite'.format(nlbg))
pl.xlim(1000., 1650.)
pl.ylim(-0.5,  13.00)
pl.xlabel('Restframe wavelength [A]')
    
pl.show()
