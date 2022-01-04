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
from   linefit.cframe_postage    import cframe_postage
from   linefit.postage_seriesfit import series_fit, plot_postages
from   desispec.io               import read_spectra
from   linefit.doublet_obs       import doublet_obs
from   linefit.doublet_postagefit import doublet_fit
from   linefit.lines import lines


dchi2     = 25.
tile      = 80869

nlbg      = 0
compwave  = np.arange(1000., 5000., 1.)
comp      = []

cframes   = {}

for petal in [0,1,2,3,4,5,7,8,9]:
    dat    = Table.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/zbest-{:d}-{:d}-deep.fits'.format(tile, petal, tile), 'ZBEST')
    dat    = dat[dat['DELTACHI2'] > 25.]
        
    spectra = read_spectra('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/coadd-{:d}-{:d}-deep.fits'.format(tile, petal, tile))

    fmap    = spectra.fibermap

    fmap['ROW'] = np.arange(len(fmap))
        
    isin   = (fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_MAIN']) != 0
    isin  |= (fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_HP'])   != 0	

    isin  &= (fmap['FIBERSTATUS'] == 0)

    fmap   = fmap[isin]
    fmap   = join(fmap, dat, keys='TARGETID', join_type='left')
    
    # Cut on dchisq
    fmap   = fmap[~fmap['Z'].mask]

    # Cut on dchisq.
    fmap   = fmap[fmap['Z'] > 2.]

    for index in range(len(fmap)):
        rrz    = fmap['Z'][index]
        rrzerr = fmap['ZERR'][index]
        fiber  = fmap['ROW'][index]
        
        wave, flux, tflux, ilineflux, ilineflux_err = cframe_postage(spectra, fiber, rrz)
        
        # print('{:d}\t{:d}\t{:.4f}\t{:.4f}\t{:.4f}'.format(petal, fiber, rrz, ilineflux, ilineflux_err))

        x0 = np.array([rrz, 140., np.log(10.)])

        doublet_fit(x0, fiber, rrz, rrzerr, spectra, lineb=1215.67, plot=True)
