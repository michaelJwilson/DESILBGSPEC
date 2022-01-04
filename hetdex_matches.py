import numpy as np
import pylab as pl
import astropy.io.fits as fits

from   astropy.table import Table, join, vstack
from   astropy.coordinates import SkyCoord
from   astropy import units as u


cat = None

for petal in range(10):
    fpath = '/global/cscratch1/sd/mjwilson/DESILBGSPEC/80869/cumulative/20210408/v1.1/zbest-{}-80869-thru20210408.fits'.format(petal)
    zbest = Table.read(fpath)
    fmap  = Table.read(fpath, 'FIBERMAP') 

    zbest = join(zbest, fmap, join_type='left', keys='TARGETID')
    zbest['PETAL'] = petal
    
    if cat is None:
        cat = zbest

    else:
        cat = vstack((cat, zbest))
        
#
# hp   = Table.read('/project/projectdirs/desi/target/secondary/sv1/indata/HETDEX_HP.fits')
# main = Table.read('/project/projectdirs/desi/target/secondary/sv1/indata/HETDEX_MAIN.fits') 

# main.pprint()

hetdex = Table.read('/global/cscratch1/sd/mjwilson/HETDEX/hetdex_matches_80869.fits')

c    = SkyCoord(ra=cat['TARGET_RA']*u.degree, dec=cat['TARGET_DEC']*u.degree)
c2   = SkyCoord(ra=hetdex['RA'].data*u.degree, dec=hetdex['DEC'].data*u.degree)

idx, d2d, d3d = c2.match_to_catalog_sky(c)

desi = cat[idx]
isin = (desi['DELTACHI2'].data >= 25.) & (desi['FIBERSTATUS'].data == 0)

zs   = np.arange(0.0, 4.5, 0.01)

ds   = (1216. / 3727.) * (1. + zs) 
ds  -= 1.

interlopers = (desi['Z'].data[isin] < 0.5) & (desi['Z'].data[isin] > .3) & (hetdex['z_guess'].data[isin] > 3.) & (hetdex['z_guess'].data[isin] < 3.2)
# interlopers = (desi['Z'].data[isin] < 0.05) & (hetdex['z_guess'].data[isin] > 2.5) & (hetdex['z_guess'].data[isin] < 3.)
# interlopers   = (desi['Z'].data[isin] > 0.77) & (desi['Z'].data[isin] < 0.79)

pl.plot(hetdex['z_guess'].data[isin], desi['Z'].data[isin], marker=',', lw=0.0, c='k')
pl.plot(hetdex['z_guess'].data[isin][interlopers], desi['Z'].data[isin][interlopers], marker=',', lw=0.0, c='r')
pl.plot(zs, ds, lw=0.1, c='k')
pl.xlabel('HETDEX Z')
pl.ylabel('DESI Z')
pl.show()

desi[isin][interlopers]['TARGETID', 'PETAL', 'Z', 'ZWARN', 'FIBERSTATUS'].pprint()
