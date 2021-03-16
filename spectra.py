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

os.environ['RR_TEMPLATE_DIR'] = '/global/cscratch1/sd/mjwilson/DESILBGSPEC/templates/'

templates = dict()

for filename in redrock.templates.find_templates():
    t = redrock.templates.Template(filename)
    templates[(t.template_type, t.sub_type)] = t

skylines = np.array([5199.4,5578.4,5656.4,5891.4,5897.4,6302.4,6308.4,6365.4,6500.4,6546.4,\
                      6555.4,6618.4,6663.4,6679.4,6690.4,6765.4,6831.4,6836.4,6865.4,6925.4,\
                      6951.4,6980.4,7242.4,7247.4,7278.4,7286.4,7305.4,7318.4,7331.4,7343.4,\
                      7360.4,7371.4,7394.4,7404.4,7440.4,7526.4,7714.4,7719.4,7752.4,7762.4,\
                      7782.4,7796.4,7810.4,7823.4,7843.4,7855.4,7862.4,7873.4,7881.4,7892.4,\
                      7915.4,7923.4,7933.4,7951.4,7966.4,7982.4,7995.4,8016.4,8028.4,8064.4,\
                      8280.4,8284.4,8290.4,8298.4,8301.4,8313.4,8346.4,8355.4,8367.4,8384.4,\
                      8401.4,8417.4,8432.4,8454.4,8467.4,8495.4,8507.4,8627.4,8630.4,8634.4,\
                      8638.4,8652.4,8657.4,8662.4,8667.4,8672.4,8677.4,8683.4,8763.4,8770.4,\
                      8780.4,8793.4,8829.4,8835.4,8838.4,8852.4,8870.4,8888.4,8905.4,8922.4,\
                      8945.4,8960.4,8990.4,9003.4,9040.4,9052.4,9105.4,9227.4,9309.4,9315.4,\
                      9320.4,9326.4,9340.4,9378.4,9389.4,9404.4,9422.4,9442.4,9461.4,9479.4,\
                      9505.4,9521.4,9555.4,9570.4,9610.4,9623.4,9671.4,9684.4,9693.4,9702.4,\
                      9714.4,9722.4,9740.4,9748.4,9793.4,9802.4,9814.4,9820.4])

dchi2     = 25.

tile      = 80869
petal     = sys.argv[1]

for band in ['B', 'R', 'Z']:
    dat   = Table.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/zbest-{:d}-{:d}-deep.fits'.format(tile, petal, tile))
    dat   = dat[dat['DELTACHI2'] > 25.]
    
    fmap  = fitsio.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/coadd-{:d}-{:d}-deep.fits'.format(tile, petal, tile), 'FIBERMAP')
    wave  = fitsio.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/coadd-{:d}-{:d}-deep.fits'.format(tile, petal, tile), '{}_WAVELENGTH'.format(band))
    flux  = fitsio.read('/global/cscratch1/sd/mjwilson/DESILBGSPEC/{:d}/deep/coadd-{:d}-{:d}-deep.fits'.format(tile, petal, tile), '{}_FLUX'.format(band))

    isin  = (fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_MAIN']) != 0
    isin |= (fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_HP']) != 0	

    isin &= (fmap['FIBERSTATUS'] == 0)

    flux  = flux[isin]
    fmap  = fmap[isin]

    fmap  = join(fmap, dat, keys='TARGETID', join_type='left')

    nlbg  = np.count_nonzero((fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_MAIN']) != 0)
    nlbg += np.count_nonzero((fmap['SV1_SCND_TARGET'] & scnd_mask['HETDEX_HP'])   != 0)

    # Cut on dchisq.
    flux  = flux[~fmap['Z'].mask]
    fmap  = fmap[~fmap['Z'].mask]

    # fmap.pprint()
    
    if band == 'B':
        ngood        = len(fmap)
        nhiz         = np.count_nonzero(fmap['Z'] > 2.0)

        print('petal: {}'.format(petal))
        print('nlbg: {}'.format(nlbg))
        print('ngood: {}'.format(ngood))
        print('nhiz: {}'.format(nhiz))

        gauss_kernel = Gaussian1DKernel(2)
        
        nplot        = np.minimum(ngood, 50)
        
        fig, axes    = plt.subplots(nplot, 1, figsize=(25, 5*nplot))

        for ax in axes:
            for sky in skylines:
                ax.axvline(sky, c='c', lw=0.5, alpha=0.5)
    
    for i in range(nplot):
        is_main   = (fmap['SV1_SCND_TARGET'][i] & scnd_mask['HETDEX_MAIN']) != 0
        is_hp     = (fmap['SV1_SCND_TARGET'][i] & scnd_mask['HETDEX_HP'])   != 0
    
        cflux     = convolve(flux[i,:], gauss_kernel)

        if band == 'B':
            axes[i].axvspan(4300., 4500., alpha=0.2, color='red')
        
        if band != 'Z':
            axes[i].plot(wave, cflux, lw=0.5, alpha=0.75, label='', c='k')

        else:
            label = '{}  z: {:.2f} (DX2: {:.2f}, MAIN: {:d}, HP: {:d})'.format(fmap['TARGETID'][i], fmap['Z'][i], fmap['DELTACHI2'][i], is_main, is_hp)

            axes[i].axvline(1216. * (1. + fmap['Z'][i]), c='m', lw=0.1, alpha=0.1, linestyle='--')

            spectype  = fmap['SPECTYPE'][i].strip()                                                                                                                                                                                         
            subtype   = fmap['SUBTYPE'][i].strip()
            fulltype  = (spectype, subtype)                                                                                                                                                                                                  
            ncoeff    = templates[fulltype].flux.shape[0]                                                                                                                                                                                   
            coeff     = fmap['COEFF'][i][0:ncoeff]                                                                                                                                                                                        
            tflux     = templates[fulltype].flux.T.dot(coeff)                                                                                                                                                                               
            twave     = templates[fulltype].wave * (1. + fmap['Z'][i])

            axes[i].plot(wave,  cflux, lw=0.5, alpha=0.75, label=label, c='k')
            axes[i].plot(twave, tflux, lw=0.5, alpha=0.5, c='m')

for ax in axes:
    ax.set_yscale('linear')
    ax.set_xlim(left=3600., right=1.e4)
    ax.set_ylim(bottom=-0.5, top=2.)
    ax.legend(frameon=False, loc=1)            
     
fig.suptitle('DESIHETDEX: {:d} spectra (${:d} @ \Delta \chi^2 > {}$, {:d} @ z>2) '.format(nlbg, ngood, dchi2, nhiz), y=1.0)
plt.tight_layout()

# /global/cfs/cdirs/desi/www/users
pl.savefig('/project/projectdirs/desi/users/mjwilson/DESILBGSPEC/{}/redshifts_{:d}.pdf'.format(tile, petal))

# os.system('chmod --reference=/project/projectdirs/desi/users/mjwilson/plots /project/projectdirs/desi/users/mjwilson/DESILBGSPEC/{}/'.format(tile))
# os.system('chmod --reference=/project/projectdirs/desi/users/mjwilson/plots /project/projectdirs/desi/users/mjwilson/DESILBGSPEC/{}/*'.format(tile))
