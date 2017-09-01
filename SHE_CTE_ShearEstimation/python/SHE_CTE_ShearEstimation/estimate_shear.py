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

from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT.shear_estimates_table_format import initialise_shear_estimates_table, tf as setf
from SHE_PPT.she_image import SHEImage
from SHE_PPT.magic_values import scale_label, stamp_size_label

from SHE_CTE_ShearEstimation import magic_values as mv

class ShearEstimate(object):
    def __init__(self, g1, g2, gerr=None):
        self.g1 = g1
        self.g2 = g2
        self.gerr = gerr
        
def get_resampled_image(subsampled_image, resampled_scale):
    
    ss_scale = subsampled_image.header[scale_label]
    
    resampled_nx = int(np.shape(subsampled_image.data)[0] / (resampled_scale/ss_scale))
    resampled_ny = int(np.shape(subsampled_image.data)[1] / (resampled_scale/ss_scale))
    
    resampled_gs_image = galsim.Image(resampled_nx,resampled_ny, scale = resampled_scale)
    
    galsim.InterpolatedImage(galsim.Image(subsampled_image.data,
                                          scale=subsampled_image.header[scale_label])).drawImage(resampled_gs_image)
    
    resampled_image = SHEImage(resampled_gs_image.data)
    resampled_image.header[scale_label] = resampled_scale
    
    return resampled_image
    
def KSB_estimate_shear(*args,**kwargs):
    return GS_estimate_shear(method="KSB",*args,**kwargs)
    
def REGAUSS_estimate_shear(*args,**kwargs):
    return GS_estimate_shear(method="REGAUSS",*args,**kwargs)

def get_KSB_shear_estimate(galsim_shear_estimate):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_KSB_shear_estimate")
    
    g1 = galsim_shear_estimate.corrected_g1
    g2 = galsim_shear_estimate.corrected_g2
    mag = g1 ** 2 + g2 ** 2
    
    if mag > 1:
        raise RuntimeError("HSM Error: Magnitude of g shear is too large: " + str(mag))
    
    if np.abs(galsim_shear_estimate.corrected_shape_err) < 1e99:
        shape_err = np.sqrt(shape_noise_var + galsim_shear_estimate.corrected_shape_err ** 2)
    else:
        shape_err = galsim_shear_estimate.corrected_shape_err
        
    shear_estimate = ShearEstimate(galsim_shear_estimate.corrected_g1, 
        galsim_shear_estimate.corrected_g2, 
        shape_err)
        
    logger.debug("Exiting get_KSB_shear_estimate")
    
    return shear_estimate


def get_REGAUSS_shear_estimate(galsim_shear_estimate):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_REGAUSS_shear_estimate")
    
    e1 = galsim_shear_estimate.corrected_e1
    e2 = galsim_shear_estimate.corrected_e2
    mag = e1 ** 2 + e2 ** 2
    
    if mag > 1:
        raise RuntimeError("HSM Error: Magnitude of e shear is too large: " + str(mag))
    
    g1, g2 = get_g_from_e(e1, e2)
    gerr = galsim_shear_estimate.corrected_shape_err * np.sqrt((g1 ** 2 + g2 ** 2) / (e1 ** 2 + e2 ** 2))
    
    if np.abs(galsim_shear_estimate.corrected_shape_err) < 1e99:
        shape_err = np.sqrt(shape_noise_var + gerr ** 2)
    else:
        shape_err = gerr
        
    shear_estimate = ShearEstimate(g1, g2, shape_err)
        
    logger.debug("Exiting get_REGAUSS_shear_estimate")
    
    return shear_estimate

def get_shear_estimate(gal_stamp, psf_stamp, sky_var):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_shear_estimate")
    
    # Get a resampled PSF stamp
    resampled_psf_stamp = get_resampled_image(psf_stamp, gal_stamp.scale)
    
    try:
        
        galsim_shear_estimate = galsim.hsm.EstimateShear(gal_image=galsim.Image(gal_stamp.data.transpose(), scale=gal_stamp.scale), 
                                                         PSF_image=galsim.Image(resampled_psf_stamp.data.transpose(), 
                                                                                scale=gal_stamp.scale), 
                                                         sky_var=sky_var, 
                                                         guess_sig_gal=0.5 / galaxy_image.scale, 
                                                         guess_sig_PSF=0.2 / resampled_psf_image.scale, 
                                                         shear_est=method)
        
        if method == "KSB":
            
            shear_estimate = get_KSB_shear_estimate(galsim_shear_estimate)
            
        elif method == "REGAUSS":
            
            shear_estimate = get_REGAUSS_shear_estimate(galsim_shear_estimate)
            
        else:
            raise RuntimeError("Invalid shear estimation method for GalSim: " + str(method))
        
    except RuntimeError as e:
        if ("HSM Error" not in str(e)):
            raise
        logger.debug(str(e))
        shear_estimate = ShearEstimate(np.NaN, np.NaN, 1e99)
        
    logger.debug("Exiting get_shear_estimate")
    
    return shear_estimate

def inv_var_stack( a, a_err ):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering inv_var_stack")
    
    a_inv_var = 1/a_err**2
    
    inv_a_inv_var_sum = 1./a_inv_var.sum()
    
    a_m = (a*a_inv_var).sum()*inv_a_inv_var_sum
    
    a_m_err = inv_a_inv_var_sum
    
    return a_m, a_m_err
        
    logger.debug("Exiting inv_var_stack")

def GS_estimate_shear( data_stack ):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering estimate_shear_gs")
    
    # Get lists of exposures, PSF images, and detection tables
    data_images = data_stack.get_data_images()
    psf_images = data_stack.get_psf_images()
    detection_tables = data_stack.get_detection_tables_images()
    
    # Calculate the sky variances
    sky_vars = []
    for detection_table in detection_tables:
        sky_vars.append(get_var_ADU_per_pixel(pixel_value_ADU=0.,
                                                sky_level_ADU_per_sq_arcsec=detections_table[detf.m.subtracted_sky_level],
                                                read_noise_count=detections_table[detf.m.read_noise],
                                                pixel_scale=data_image.scale,
                                                gain=detections_table[detf.m.gain]))
    
    num_tables = len(detections_tables)
    
    # Get all unique IDs
    IDs = None
    for table_index in range(num_tables):
        if IDs is None:
            IDs = set(detections_tables[table_index][detf.ID])
        else:
            IDs = set.union(IDs,detections_tables[table_index][detf.ID])
            
    # Set the ID as an index for each table
    for detection_table in detection_tables:
        detection_table.add_index(detf.ID)
    
    shear_estimates_table = initialise_shear_estimates_table(detections_tables[0])
    
    for ID in IDs:
        
        g1s = []
        g2s = []
        gerrs = []
        
        for detections_table in detections_tables:
            
            try:
                
                # Get the row for this ID
                row = detections_table.loc[ID]
            
                # Get galaxy and PSF stamps
                gal_stamp = data_image.extract_stamp(row[detf.gal_x],
                                                     row[detf.gal_y],
                                                     data_image.header[stamp_size_label])
                psf_stamp = psf_image.extract_stamp(row[detf.psf_x],
                                                    row[detf.psf_y],
                                                    psf_image.header[stamp_size_label])
        
                shear_estimate = get_shear_estimate(gal_stamp, psf_stamp, sky_vars[table_index])
                
                g1s.append(shear_estimate.g1)
                g2s.append(shear_estimate.g2)
                gerrs.append(shear_estimate.gerr)
                
            except KeyError as e:
                if "No matches found for key" in e:
                    pass # ID isn't present in this table, so just skip it
                else:
                    raise
                
        g1s = np.array(g1s)
        g2s = np.array(g2s)
        gerrs = np.array(gerrs)
        
        g1, gerr1 = inv_var_stack(g1s,gerrs)
        g2, gerr2 = inv_var_stack(g2s,gerrs)
        
        assert np.isclose(gerr1,gerr2)
            
        # Add this row to the estimates table
        shear_estimates_table.add_row({ setf.ID : detections_table[detf.ID][i],
                                        setf.gal_g1 : g1,
                                        setf.gal_g2 : g2,
                                        setf.gal_g1_err : gerr1,
                                        setf.gal_g2_err : gerr2,
                                        setf.gal_e1_err : np.NaN,
                                        setf.gal_e2_err : np.NaN,
                                       })
        
    
    logger.debug("Exiting GS_estimate_shear")
    
    return shear_estimates_table, None # No MCMC chains for this method
    
