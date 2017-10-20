""" @file galsim_estimate_shear.py

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

from math import sqrt

import galsim

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_GST_GalaxyImageGeneration.noise import get_var_ADU_per_pixel
from SHE_GST_GalaxyImageGeneration.unweighted_moments import get_g_from_e
from SHE_GST_IceBRGpy.logging import getLogger
from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT.magic_values import scale_label, stamp_size_label
from SHE_PPT.psf_table_format import tf as pstf
from SHE_PPT.she_image import SHEImage
from SHE_PPT.shear_estimates_table_format import initialise_shear_estimates_table, tf as setf
import numpy as np


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
                                          scale=subsampled_image.header[scale_label])).drawImage(resampled_gs_image,use_true_center=False)
    
    resampled_image = SHEImage(resampled_gs_image.array)
    resampled_image.header[scale_label] = resampled_scale
    
    return resampled_image
    
def KSB_estimate_shear(data_stack,method_data):
    return GS_estimate_shear(data_stack,method="KSB")
    
def REGAUSS_estimate_shear(data_stack,method_data):
    return GS_estimate_shear(data_stack,method="REGAUSS")

def get_KSB_shear_estimate(galsim_shear_estimate):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_KSB_shear_estimate")
    
    g1 = galsim_shear_estimate.corrected_g1
    g2 = galsim_shear_estimate.corrected_g2
    mag = g1 ** 2 + g2 ** 2
    
    if mag > 1:
        raise RuntimeError("HSM Error: Magnitude of g shear is too large: " + str(mag))
        
    shear_estimate = ShearEstimate(galsim_shear_estimate.corrected_g1, 
        galsim_shear_estimate.corrected_g2, 
        galsim_shear_estimate.corrected_shape_err)
        
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
        
    shear_estimate = ShearEstimate(g1, g2, gerr)
        
    logger.debug("Exiting get_REGAUSS_shear_estimate")
    
    return shear_estimate

def get_shear_estimate(gal_stamp, psf_stamp, sky_var, method):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_shear_estimate")
    
    # Get a resampled PSF stamp
    resampled_psf_stamp = get_resampled_image(psf_stamp, gal_stamp.header[scale_label])
    
    try:
        
        galsim_shear_estimate = galsim.hsm.EstimateShear(gal_image=galsim.Image(gal_stamp.data.transpose(), scale=gal_stamp.header[scale_label]), 
                                                         PSF_image=galsim.Image(resampled_psf_stamp.data.transpose(), 
                                                                                scale=gal_stamp.header[scale_label]), 
                                                         sky_var=sky_var, 
                                                         guess_sig_gal=0.5 / gal_stamp.header[scale_label], 
                                                         guess_sig_PSF=0.2 / resampled_psf_stamp.header[scale_label], 
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
    
    a_m_err = sqrt(inv_a_inv_var_sum)
    
    return a_m, a_m_err
        
    logger.debug("Exiting inv_var_stack")

def GS_estimate_shear( data_stack, method ):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering estimate_shear_gs")
    
    # Get lists of exposures, PSF images, and detection tables
    data_images = []
    detection_tables = []
    bulge_psf_images = []
    disk_psf_images = []
    psf_tables = []
    
    num_exposures = len(data_stack.exposures)
    
    for i in range(num_exposures):
        data_images.append(data_stack.exposures[i].science_image)
        detection_tables.append(data_stack.exposures[i].detections_table)
        bulge_psf_images.append(data_stack.exposures[i].bpsf_image)
        disk_psf_images.append(data_stack.exposures[i].dpsf_image)
        psf_tables.append(data_stack.exposures[i].psf_table)
    
    # Calculate the sky variances
    sky_vars = []
    for table_index in range(num_exposures):
        detections_table = detection_tables[table_index]
        sky_vars.append(get_var_ADU_per_pixel(pixel_value_ADU=0.,
                                                sky_level_ADU_per_sq_arcsec=detections_table.meta[detf.m.subtracted_sky_level],
                                                read_noise_count=detections_table.meta[detf.m.read_noise],
                                                pixel_scale=data_images[table_index].header[scale_label],
                                                gain=detections_table.meta[detf.m.gain]))
    
    # Get all unique IDs
    IDs = None
    for table_index in range(num_exposures):
        if IDs is None:
            IDs = set(detection_tables[table_index][detf.ID])
        else:
            IDs = set.union(IDs,detection_tables[table_index][detf.ID])
            
    # Set the ID as an index for each table
    for table_index in range(num_exposures):
        detection_tables[table_index].add_index(detf.ID)
        psf_tables[table_index].add_index(pstf.ID)
    
    shear_estimates_table = initialise_shear_estimates_table(detection_tables[0],
                                                             optional_columns=[setf.e1_err,setf.e2_err])
    
    for ID in IDs:
        
        g1s = []
        g2s = []
        gerrs = []
        
        for table_index in range(num_exposures):
            
            try:
                
                # Get the rows for this ID
                g_row = detection_tables[table_index].loc[ID]
                p_row = psf_tables[table_index].loc[ID]
            
                # Get galaxy and PSF stamps
                gal_stamp = data_images[table_index].extract_stamp(g_row[detf.gal_x],
                                                                   g_row[detf.gal_y],
                                                                   data_images[table_index].header[stamp_size_label],
                                                                   keep_header=True)
                bulge_psf_stamp = bulge_psf_images[table_index].extract_stamp(p_row[pstf.psf_x],
                                                                              p_row[pstf.psf_y],
                                                                              bulge_psf_images[table_index].header[stamp_size_label],
                                                                              keep_header=True)
                disk_psf_stamp = disk_psf_images[table_index].extract_stamp(p_row[pstf.psf_x],
                                                                              p_row[pstf.psf_y],
                                                                              disk_psf_images[table_index].header[stamp_size_label],
                                                                              keep_header=True)
        
                shear_estimate = get_shear_estimate(gal_stamp, bulge_psf_stamp, sky_vars[table_index], method)
                
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
        shear_estimates_table.add_row({ setf.ID : detection_tables[table_index][detf.ID][i],
                                        setf.g1 : g1,
                                        setf.g2 : g2,
                                        setf.e1_err : gerr1,
                                        setf.e2_err : gerr2,
                                       })
        
    
    logger.debug("Exiting GS_estimate_shear")
    
    return shear_estimates_table, None # No MCMC chains for this method
    
