""" @file estimate_shear.py

    Created 27 Mar 2017

    Provides functions to measure the shape of a galaxy image.

    ---------------------------------------------------------------------

    Copyright (C) 2012-2020 Euclid Science Ground Segment      
       
    This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General    
    Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)    
    any later version.    
       
    This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied    
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more    
    details.    
       
    You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to    
    the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
"""

import numpy as np
import galsim

from SHE_GST_IceBRGpy.logging import getLogger

from SHE_GST_GalaxyImageGeneration.noise import get_var_ADU_per_pixel
from SHE_GST_GalaxyImageGeneration.unweighted_moments import get_g_from_e

from SHE_CTE_ShearEstimation import magic_values as mv

class ShearEstimate(object):
    def __init__(self, g1, g2, gerr=None):
        self.g1 = g1
        self.g2 = g2
        self.gerr = gerr
        
def get_resampled_image(subsampled_psf_image, scale):
    
    resampled_nx = int(np.shape(subsampled_psf_image.array)[0] / (scale/subsampled_psf_image.scale))
    resampled_ny = int(np.shape(subsampled_psf_image.array)[1] / (scale/subsampled_psf_image.scale))
    
    resampled_image = galsim.Image(resampled_nx,resampled_ny)
    
    galsim.InterpolatedImage(subsampled_psf_image).drawImage(resampled_image)
    
    return resampled_image

def estimate_shear(method,*args,**kwargs):
    """
    @brief
        Determines which method to use for measuring shape and calls that
    
    @param method <string> Name of the method to call
    
    @param args <list> List of arguments to pass on
    
    @param kwargs <dict> Dictionary of keyword arguments to pass on
    
    @return <measurement> Shape measurement object
    """
    
    # Check the method and call the appropriate function
    if method.lower() == 'ksb':
        return estimate_shear_gs(method="KSB",*args,**kwargs)
    if method.lower() == 'regauss':
        return estimate_shear_gs(method="REGAUSS",*args,**kwargs)
    else:
        raise Exception("Invalid shear estimation method: " + str(method))
    
def estimate_shear_gs(galaxy_image, psf_image, gain, subtracted_sky_level,
                       read_noise, shape_noise_var, method):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering estimate_shear_gs")
    
    # Calculate the sky variance
    sky_var = get_var_ADU_per_pixel(pixel_value_ADU=0.,
                                    sky_level_ADU_per_sq_arcsec=subtracted_sky_level,
                                    read_noise_count=read_noise,
                                    pixel_scale=galaxy_image.scale,
                                    gain=gain)
    
    # Get a resampled PSF image
    resampled_psf_image = get_resampled_image(psf_image, galaxy_image.scale)
    try:
        galsim_shear_estimate = galsim.hsm.EstimateShear(gal_image=galaxy_image,
                                                         PSF_image=resampled_psf_image,
                                                         sky_var=sky_var,
                                                         guess_sig_gal=0.5/galaxy_image.scale,
                                                         guess_sig_PSF=0.2/resampled_psf_image.scale,
                                                         shear_est=method)
        
        if method=="KSB":
            g1 = galsim_shear_estimate.corrected_g1
            g2 = galsim_shear_estimate.corrected_g2
            mag = g1**2 + g2**2
            if mag > 1:
                raise RuntimeError("HSM Error: Magnitude of g shear is too large: " + str(mag))
        
            if np.abs(galsim_shear_estimate.corrected_shape_err) < 1e99:
                shape_err = np.sqrt(shape_noise_var+galsim_shear_estimate.corrected_shape_err**2)
            else:
                shape_err = galsim_shear_estimate.corrected_shape_err
                
            shear_estimate = ShearEstimate(galsim_shear_estimate.corrected_g1,
                                           galsim_shear_estimate.corrected_g2,
                                           shape_err,)
        elif method=="REGAUSS":
            e1 = galsim_shear_estimate.corrected_e1
            e2 = galsim_shear_estimate.corrected_e2
            mag = e1**2 + e2**2
            if mag > 1:
                raise RuntimeError("HSM Error: Magnitude of e shear is too large: " + str(mag))
            g1, g2 = get_g_from_e(e1,e2)
            gerr = galsim_shear_estimate.corrected_shape_err * np.sqrt((g1**2+g2**2)/(e1**2+e2**2))
        
            if np.abs(galsim_shear_estimate.corrected_shape_err) < 1e99:
                shape_err = np.sqrt(shape_noise_var+gerr**2)
            else:
                shape_err = gerr
                
            shear_estimate = ShearEstimate(g1, g2, shape_err,)
        else:
            raise RuntimeError("Invalid shear estimation method for GalSim: " + str(method))
            
            
    except RuntimeError as e:
        
        if("HSM Error" not in str(e)):
            raise
        
        logger.debug(str(e))
        
        shear_estimate = ShearEstimate(0, 0, 1e99)
        
    
    logger.debug("Exiting estimate_shear_gs")
    
    return shear_estimate
    