""" @file galsim_estimate_shear.py

    Created 27 Mar 2017

    Provides functions to measure the shape of a galaxy image.
"""

# Copyright (C) 2012-2020 Euclid Science Ground Segment      
#        
# This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General    
# Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)    
# any later version.    
#        
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied    
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more    
# details.    
#        
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to    
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

from math import sqrt

import galsim

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_PPT.noise import get_var_ADU_per_pixel
from SHE_PPT.shear_utility import get_g_from_e
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.magic_values import scale_label, stamp_size_label
from SHE_PPT.table_formats.psf import tf as pstf
from SHE_PPT.she_image import SHEImage
from SHE_PPT.table_formats.shear_estimates import initialise_shear_estimates_table, tf as setf
import numpy as np

stamp_size = 256
x_buffer = -5
y_buffer = -5
shape_noise = 0.25

class ShearEstimate(object):
    def __init__(self, g1, g2, gerr=None, re=None, snr=None, x=None, y=None):
        self.g1 = g1
        self.g2 = g2
        self.gerr = gerr
        self.re = re
        self.snr = snr
        self.x = x
        self.y = y
        
def get_resampled_image(subsampled_image, resampled_scale):
    
    ss_scale = subsampled_image.header[scale_label]
    
    resampled_nx = int(np.shape(subsampled_image.data)[0] / (resampled_scale/ss_scale))
    resampled_ny = int(np.shape(subsampled_image.data)[1] / (resampled_scale/ss_scale))
    
    resampled_gs_image = galsim.Image(resampled_nx,resampled_ny, scale = resampled_scale)
    
    galsim.InterpolatedImage(galsim.Image(subsampled_image.data,
                                          scale=subsampled_image.header[scale_label])).drawImage(resampled_gs_image,use_true_center=True)
    
    resampled_image = SHEImage(resampled_gs_image.array)
    resampled_image.header[scale_label] = resampled_scale
    
    return resampled_image
    
def KSB_estimate_shear(data_stack,method_data,workdir):
    return GS_estimate_shear(data_stack,method="KSB",workdir)
    
def REGAUSS_estimate_shear(data_stack,method_data,workdir):
    return GS_estimate_shear(data_stack,method="REGAUSS",workdir)

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
        galsim_shear_estimate.corrected_shape_err,
        galsim_shear_estimate.moments_sigma,
        galsim_shear_estimate.moments_amp,
        galsim_shear_estimate.moments_centroid.x,
        galsim_shear_estimate.moments_centroid.y)
        
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
        
    shear_estimate = ShearEstimate(g1, g2, gerr,
        galsim_shear_estimate.moments_sigma,
        galsim_shear_estimate.moments_amp,
        galsim_shear_estimate.moments_centroid.x,
        galsim_shear_estimate.moments_centroid.y)
        
    logger.debug("Exiting get_REGAUSS_shear_estimate")
    
    return shear_estimate

def get_shear_estimate(gal_stamp, psf_stamp, gal_scale, psf_scale, ID, method):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_shear_estimate")
    
    # Get a resampled PSF stamp
    resampled_psf_stamp = get_resampled_image(psf_stamp, gal_scale)
    
    gal_mask = gal_stamp.get_object_mask(ID).astype(np.uint16) # Galsim requires int array
    sky_var = np.square(gal_stamp.noisemap.transpose())
    
    try:
        
        galsim_shear_estimate = galsim.hsm.EstimateShear(gal_image=galsim.Image(gal_stamp.data.transpose(), scale=gal_scale), 
                                                         PSF_image=galsim.Image(resampled_psf_stamp.data.transpose(), 
                                                                                scale=gal_scale), 
                                                         badpix=galsim.Image(gal_mask.transpose(), scale=gal_scale),
                                                         sky_var=sky_var, 
                                                         guess_sig_gal=0.5 / gal_scale, 
                                                         guess_sig_PSF=0.2 / gal_scale, 
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
    
    a_m = (a*a_inv_var).nansum()*inv_a_inv_var_sum
    
    a_m_err = sqrt(inv_a_inv_var_sum)
    
    return a_m, a_m_err
        
    logger.debug("Exiting inv_var_stack")

def GS_estimate_shear( data_stack, method, workdir ):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering GS_estimate_shear")
    
    shear_estimates_table = initialise_shear_estimates_table(data_stack.detections_catalogue,
                                                             optional_columns=[setf.e1_err,setf.e2_err])

    gal_scale = data_stack.exposures[0].detectors[0,0].header[scale_label]
    psf_scale = data_stack.exposures[0].bulge_psf_image.header[scale_label]
    
    # Loop over galaxies and get an estimate for each one
    for row in data_stack.detections_catalogue:
        
        gal_id = row[detf.ID]
        gal_x_world = row[detf.gal_x_world]
        gal_y_world = row[detf.gal_y_world]
        
        gal_stamp_stack = data_stack.extract_stamp_stack(x_world=gal_x_world,
                                                         y_world=gal_y_world,
                                                         width=stamp_size,
                                                         x_buffer=x_buffer,
                                                         y_buffer=y_buffer,)
        
        g1s = []
        g2s = []
        gerrs = []
        res = []
        snrs = []
        x_worlds = []
        y_worlds = []
        
        for ex_i in range(len(data_stack.exposures)):
            
            gal_stamp = gal_stamp_stack.exposures[ex_i]
            
            if gal_stamp is None:
                g1s.append(np.NaN)
                g2s.append(np.NaN)
                gerrs.append(1e99)
                continue
            
            bulge_psf_stamp, disk_psf_stamp = data_stack.exposures[ex_i].extract_psf(gal_id)
    
            shear_estimate = get_shear_estimate(gal_stamp,
                                                bulge_psf_stamp,
                                                gal_scale=gal_scale,
                                                psf_scale=psf_scale,
                                                ID=gal_id,
                                                method=method)
            
            g1s.append(shear_estimate.g1)
            g2s.append(shear_estimate.g2)
            gerrs.append(shear_estimate.gerr)
            res.append(shear_estimate.re)
            snrs.append(shear_estimate.snr)
            
            x_world, y_world = data_stack.exposures[ex_i].pix2world(shear_estimate.x,shear_estimate.y)
            
            x_worlds.append(x_world)
            y_worlds.append(y_world)
                
        g1s = np.array(g1s)
        g2s = np.array(g2s)
        gerrs = np.array(gerrs)
        res = np.array(res)
        snrs = np.array(snrs)
        x_worlds = np.array(x_worlds)
        y_worlds = np.array(y_worlds)
        
        g1, gerr = inv_var_stack(g1s,gerrs)
        g2, _ = inv_var_stack(g2s,gerrs)
        re, _ = inv_var_stack(res,gerrs)
        snr, _ = inv_var_stack(snrs,gerrs)
        x_world, _ = inv_var_stack(x_worlds,gerrs)
        y_world, _ = inv_var_stack(y_worlds,gerrs)
            
        # Add this row to the estimates table
        shear_estimates_table.add_row({ setf.ID : gal_id,
                                        setf.g1 : g1,
                                        setf.g2 : g2,
                                        setf.g1_err : np.sqrt(gerr**2+shape_noise**2),
                                        setf.g2_err : np.sqrt(gerr**2+shape_noise**2),
                                        setf.re : re,
                                        setf.snr : snr,
                                        setf.x_world : x_world,
                                        setf.y_world : y_world,
                                       })
        
    
    logger.debug("Exiting GS_estimate_shear")
    
    return shear_estimates_table
    
