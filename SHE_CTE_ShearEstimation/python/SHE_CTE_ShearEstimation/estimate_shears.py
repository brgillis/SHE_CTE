""" @file estimate_shears.py

    Created 27 Mar 2017

    Primary execution loop for measuring galaxy shapes from an image file.

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

from astropy.io import fits
from astropy.table import Table
import astropy.table

import numpy as np
import galsim

from SHE_GST_IceBRGpy.logging import getLogger

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_GST_GalaxyImageGeneration import magic_values as sim_mv

from SHE_PPT.she_image import SHEImage
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.detections_table import tf as detf
from SHE_PPT.shear_estimates_table import initialise_shear_estimates_table, tf as setf

from SHE_CTE_ShearEstimation.estimate_shear import KSB_estimate_shear, REGAUSS_estimate_shear
from SHE_CTE_ShearEstimation.output_shear_estimates import output_shear_estimates

estimation_methods = {"KSB":KSB_estimate_shear,
                      "REGAUSS":REGAUSS_estimate_shear,
                      "MegaLUT":None,
                      "LensMC":None,
                      "BFD":None}

def find_value(args_value, name, label, detections_table, galaxies_hdulist):
    if args_value is not None:
        value = args_value
    else:
        try:
            value = galaxies_hdulist[0].header[label]
        except KeyError as _e1:
            try:
                value = detections_table.meta[label]
            except KeyError as _e2:
                raise KeyError("No " + name + " value available.")
    return value

def estimate_shears_from_args(kwargs):
    """
    @brief
        Perform shear estimation, given arguments from the command-line.
    
    @param kwargs <dict>
    
    @return None
    """

    logger = getLogger(mv.logger_name)
    
    logger.debug("Entering estimate_shears_from_args")
    
    # Load the detections table
    detections_table = Table.read(kwargs["detections_table"])
    
    if not is_in_format(detections_table,detf):
        raise ValueError("Detections table " + kwargs["detections_table"] + " is in incorrect format.")
    
    # Load the various images
    data_image = SHEImage.read_from_fits(filepath = kwargs["data_image"],
                                         mask_filepath = kwargs["mask_image"],
                                         noisemap_filepath = kwargs["noise_image"],
                                         segmentation_filepath = kwargs["segmentation_image"])
    
    psf_image = SHEImage.read_from_fits(filepath = kwargs["psf_image"])
    
    # Load the P(e) table if available
    default_shape_noise_var = 0.06
    if "p_of_e_table_file_name" in kwargs:
        if kwargs["p_of_e_table_file_name"] is not None:
            p_of_e_table = Table.read(kwargs["p_of_e_table_file_name"])
            e_half_step = (p_of_e_table["E_LOW"][1] - p_of_e_table["E_LOW"][0])/2.
            shape_noise_var = (((p_of_e_table["E_LOW"]+e_half_step)**2 * p_of_e_table["E_COUNT"]).sum() /
                               p_of_e_table["E_COUNT"].sum())
        else:
            shape_noise_var = default_shape_noise_var 
    else:
        shape_noise_var = default_shape_noise_var
    
    method_shear_estimates = {}
    mcmc_chains = None
    
    if len(kwargs['methods'])==0:
        methods = estimation_methods.keys()
    else:
        methods = kwargs['methods']
    
    for method in methods:
        
        estimate_shear = estimation_methods[method]
        
        try:
            
            tab, method_mcmc_chains = method( data_image, psf_image, detections_table)
            
            if not is_in_format(tab,setf):
                raise ValueError("Shear estimation table returned in invalid format for method " + method + ".")
            
            if method_mcmc_chains is not None:
                mcmc_chains = method_mcmc_chains
                
            # Rename columns to include the method name
            for col in (setf.gal_g1, setf.gal_g2,
                        setf.gal_g1_err, setf.gal_g2_err,
                        setf.gal_e1_err, setf.gal_e2_err,):
                tab.rename_column(col,method+"_"+col)
                
        except Exception as e:
            
            logger.warning(str(e))
            
            # Create an empty estimates table
            tab = initialise_shear_estimates_table(detections_table)
            
            # Fill it with NaN measurements and 1e99 errors
            tab = method_shear_estimates[method]
            tab[setf.ID] = detections_table[detf.ID]
            tab[setf.gal_x] = detections_table[detf.gal_x]
            tab[setf.gal_y] = detections_table[detf.gal_y]
            
            for col in (method + "_" + setf.gal_g1, method + "_" + setf.gal_g2,
                        method + "_" + setf.gal_e1_err, method + "_" + setf.gal_e2_err,):
                tab[col] = np.NaN
            
            for col in (method + "_" + setf.gal_g1_err, method + "_" + setf.gal_g2_err):
                tab[col] = 1e99
            
        method_shear_estimates[method] = tab
        
    # Join the tables together
    full_shear_estimates_table = None
    
    for method in methods:
        
        new_tab = method_shear_estimates[method]
        
        if full_shear_estimates_table is None:
            
            full_shear_estimates_table = new_tab
            
        else:
            
            full_shear_estimates_table = astropy.table.join(full_shear_estimates_table,new_tab)
            
        
    logger.info("Finished estimating shear. Outputting results to " + kwargs["shear_measurements_table"] + ".")
        
    full_shear_estimates_table.write(kwargs["shear_measurements_table"],format="fits")
    
    logger.debug("Exiting estimate_shears_from_args")
    
    return