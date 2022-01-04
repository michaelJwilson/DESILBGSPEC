import time
import numba
import numpy as np
import scipy

from   scipy             import optimize
from   scipy.optimize    import approx_fprime, minimize, Bounds
from   scipy.stats       import multivariate_normal
from   jax               import grad, jit, vmap, hessian, jacfwd, jacrev
from   jax.experimental  import optimizers
from   scipy.optimize    import leastsq
from   numba             import jit

from   desispec.interpolation import resample_flux

from   .doublet import doublet
from   desispec.io import read_frame
from   desispec.io.meta import findfile
from   desispec.resolution import Resolution
from   .doublet_priors import mlogprior
from   .cframe_postage import cframe_postage
from   .lines import lines, ugroups
from   .twave import twave
from   .doublet_obs import doublet_obs
from   .matchedtemp_lineflux import matchedtemp_lineflux

width  = 25.
cwidth = 10


def doublet_chi2(z, twave, wave, res, flux, ivar, mask, continuum=0.0, sigmav=5., r=0.1, line_flux=None, linea=3726.032, lineb=3728.815):    
    '''
    Given a redshift, cframe (extracted wave, res, flux, ivar) return 
    chi sq. for a doublet line model of given parameters, e.g. line flux.

    If line flux is None, estimate it first. 
    '''
    
    rflux             = doublet_obs(z, twave, wave, res, continuum=0.0, sigmav=sigmav, r=r, linea=linea, lineb=lineb)
    rflux            *= line_flux
        
    X2                = ivar * (flux - rflux)**2.

    result            = np.sum(X2[mask == 0])
    
    return  wave, rflux, result, line_flux

def doublet_fit(x0, fiber, rrz, rrzerr, spectra, lineb, plot=False):
    '''
    Passed a postage stamp, find the best fit Gaussian (doublet).
    '''

    start           = time.time()
    cams            = spectra.flux.keys()

    center = (1. + rrz) * lineb
    limits = center + np.array([-width, width])
    
    def _X2(x):  
        z           = x[0]
        v           = x[1] 
        lnA         = x[2]
    
        line_flux   = np.exp(lnA)

        result      = 0.0
        
        for cam in cams:
            wave    = spectra.wave[cam]
            res     = Resolution(spectra.resolution_data[cam][fiber,:])
            flux    = spectra.flux[cam][fiber,:]
            ivar    = spectra.ivar[cam][fiber,:]
            mask    = spectra.mask[cam][fiber,:]
            
            continuum = (wave > limits[0]) & (wave < limits[1]) & ((wave < (limits[0] + cwidth)) | (wave > (limits[1] - cwidth)))

            if np.count_nonzero(continuum):
                continuum = np.median(flux[continuum])
                flux     -= continuum

            else:
                return  0.0
                
            tw      = twave(wave.min(), wave.max())

            _, _, X2, _ = doublet_chi2(z, tw, wave, res, flux, ivar, mask, continuum=0.0, sigmav=v, r=0.0, line_flux=line_flux, linea=0.0, lineb=lineb)

            result += X2
            
        return  result

    def _res(x):
        z           = x[0]
        v           = x[1]
        lnA         = x[2]

        line_flux   = np.exp(lnA)

        residuals   = []

        for cam in cams:
            wave    = spectra.wave[cam]
            res     = Resolution(spectra.resolution_data[cam][fiber,:])
            flux    = spectra.flux[cam][fiber,:]
            ivar    = spectra.ivar[cam][fiber,:]
            mask    = spectra.mask[cam][fiber,:]

            continuum = (wave > limits[0]) & (wave < limits[1]) & ((wave < (limits[0] + cwidth)) | (wave > (limits[1] - cwidth)))

            if np.count_nonzero(continuum):
                continuum = np.median(flux[continuum])
                flux     -= continuum
	    else:
                return  0.0
            
            tw      = twave(wave.min(), wave.max())
            
            rflux   = doublet_obs(z, tw, wave, res, continuum=0.0, sigmav=v, r=0.0, linea=0.0, lineb=lineb)
            rflux  *= line_flux

            res     = np.sqrt(ivar[mask]) * np.abs(flux[mask] - rflux[mask])

            residuals += res.tolist()

        residuals = np.array(residuals) 

        print(z, v, lnA, residuals.sum())
        
        return  residuals
        
    def mloglike(x):
        return  _X2(x) / 2.

    def mlogpos(x):
        return  mloglike(x) # + mlogprior(x, rrz, rrzerr)

    def scipy_gradient(x):
        eps = 1.e-8
        
        return optimize.approx_fprime(x, mlogpos, eps)

    # 'L-BFGS-B'; hess=jax_hessian; hessp=jax_hessian_vec_prod; 'maxiter': 10, 'maxfev': 50; 'xtol': 1.e-0; 'gtol': 1e-6                                                                                                    
    # methods        = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'trust-ncg']                                                                                                                                    
    result           = scipy.optimize.minimize(mlogpos, x0, method='Nelder-Mead', options={'disp': True, 'return_all': True})                                                                
    [z, v, lnA]      = result.x   
    
    # result         = scipy.optimize.least_squares(_res, x0, verbose=1, ftol=1.e-6, max_nfev=6, method='lm')
    # [z, v, lnA]    = result.x

    print(time.time() - start)
    
    ##
    rflux           = {}
    
    for cam in cams:
        wave        = spectra.wave[cam]
        res         = Resolution(spectra.resolution_data[cam][fiber,:])

        tw          = twave(wave.min(), wave.max())
        
        rflux[cam]  = doublet_obs(z, tw, wave, res, continuum=0.0, sigmav=v, r=0.0, linea=0.0, lineb=lineb)
        rflux[cam] *= np.exp(lnA)
        
    if plot & result.success:
        import pylab as pl

        for cam in cams:
            wave    = spectra.wave[cam]
            res     = Resolution(spectra.resolution_data[cam][fiber,:])
            flux    = spectra.flux[cam][fiber,:]
            ivar    = spectra.ivar[cam][fiber,:]
            mask    = spectra.mask[cam][fiber,:]

            continuum = (wave > limits[0]) & (wave < limits[1]) & ((wave < (limits[0] + cwidth)) | (wave > (limits[1] - cwidth)))

            if np.count_nonzero(continuum):
                continuum = np.median(flux[continuum])
                flux     -= continuum

	    else:
                return  0.0

            pl.plot(wave[mask == 0], rflux[cam][mask == 0], label='Model {}'.format(cam), lw=0.5)
            pl.plot(wave[mask == 0],  flux[mask == 0], label='Observed {}'.format(cam), alpha=0.5, lw=0.5)

            mask    = (1216. * (1. + (rrz - 0.01)) < wave) & (wave < 1216. * (1. + (rrz + 0.01)))

            if np.count_nonzero(mask):
                rough = 0.8 * np.sum(flux[mask])
            
            pl.xlim(1216. * (1. + (rrz - 0.01)), 1216. * (1. + (rrz + 0.01)))
            
            pl.xlabel('Wavlength [Angstroms]')
            pl.ylabel('Flux')

            pl.title('{:.4f} ({:.4f}) Lya @ z = {:.4f}'.format(np.exp(lnA), rough, rrz))
            
            pl.legend(frameon=False)
            
            # pl.xlim(xmin, xmax)
            pl.ylim(bottom=-0.5)
        
        pl.show()
    
    return  result, rflux
