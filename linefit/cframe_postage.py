import time 
import numpy as np
import pickle
import scipy

from   scipy             import optimize
from   scipy.optimize    import approx_fprime, minimize, Bounds
from   scipy.stats       import multivariate_normal

from   .doublet import doublet
from   desispec.io import read_frame
from   desispec.io.meta import findfile
from   desispec.resolution import Resolution
from   .doublet_priors import mlogprior
from   desispec.frame import Spectrum, Frame
from   .lines import lines, ugroups
from   .matchedtemp_lineflux import matchedtemp_lineflux
from   .doublet_obs import doublet_obs
from   .plot_postages import plot_postages


width  = 25.
cwidth = 10.

def cframe_postage(spectra, fiber, redshift, ipostage=True, printit=False):    
    '''
    Given a redshift, cframe (extracted wave, res, flux, ivar) return 
    chi sq. for a doublet line model of given parameters, e.g. line flux.
    '''
    
    sample     = lines[lines['MASKED'] == 0]
        
    for i, line in enumerate(sample['WAVELENGTH']):        
        center = (1. + redshift) * line
        limits = center + np.array([-width, width])

        name   = sample['NAME'][i]
        group  = sample['GROUP'][i]
        lratio = sample['LINERATIO'][i]
            
        for band in spectra.flux.keys():        
            wave   = spectra.wave[band]
            inwave = (wave > limits[0]) & (wave < limits[1])
            
            isin   = (wave.min() < center) & (center < wave.max())
            
            if isin:
                if printit:
                    print('Reduced LINEID {:2d}:  {:16s} for {} at redshift {:.2f} ({:.3f} to {:.3f}).'.format(sample['INDEX'][i], sample['NAME'][i], cam, redshift, limits[0], limits[1]))

                res       = spectra.resolution_data[band][fiber,:,:]
                flux      = spectra.flux[band][fiber,:]
                ivar      = spectra.ivar[band][fiber,:]
                mask      = spectra.mask[band][fiber,:]

                unmask    = mask == 0
                nelem     = np.count_nonzero(unmask[inwave])

                if nelem == 0:
                    continue
                
                continuum = (wave > limits[0]) & (wave < limits[1]) & ((wave < (limits[0] + cwidth)) | (wave > (limits[1] - cwidth)))
                continuum = np.median(flux[continuum])

                wave      = wave[inwave]
                flux      = flux[inwave] - continuum
                ivar      = ivar[inwave]
                mask      = mask[inwave]
                res       =  res[:,inwave]

                lineflux, lineflux_err, rflux = matchedtemp_lineflux(redshift, wave, Resolution(res), flux, ivar, mask, sigmav=180.0, r=0.0, linea=0.0, lineb=line)

                return  wave, flux, rflux, lineflux, lineflux_err
