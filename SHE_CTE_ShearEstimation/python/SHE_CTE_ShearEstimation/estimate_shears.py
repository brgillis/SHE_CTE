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

from SHE_PPT.she_stack import SHEStack
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.detections_table import tf as detf
from SHE_PPT.shear_estimates_table import initialise_shear_estimates_table, tf as setf
from SHE_PPT.bfd_moments_table import initialize_bfd_moments_table, bfdtf as bfdsetf

from SHE_CTE_ShearEstimation.galsim_estimate_shear import KSB_estimate_shear, REGAUSS_estimate_shear
from SHE_CTE_ShearEstimation.output_shear_estimates import output_shear_estimates
from SHE_CTE_ShearEstimation.bfd_measure_moments import bfd_measure_moments
from SHE_CTE_ShearEstimation.output_bfd_moments import output_bfd_moments

loading_methods = {"KSB":None,
                   "REGAUSS":None,
                   "MegaLUT":None,
                   "LensMC":None,
                   "BFD":None}

estimation_methods = {"KSB":KSB_estimate_shear,
                      "REGAUSS":REGAUSS_estimate_shear,
                      "MegaLUT":None,
                      "LensMC":None,
                      "BFD":BFD_measure_moments}

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
    
    # Set up a dictionary of the method data filenames
    method_data_filenames = {"KSB":kwargs["KSB_data"],
                             "REGAUSS":kwargs["REGAUSS_data"],
                             "MegaLUT":kwargs["MegaLUT_data"],
                             "LensMC":kwargs["LensMC_data"],
                             "BFD":kwargs["BFD_data"]}
    
    # Load the various images
    data_stack = SHEStack.read_from_fits(filepaths = kwargs["data_images"],
                                         mask_filepaths = kwargs["mask_images"],
                                         noisemap_filepaths = kwargs["noise_images"],
                                         segmentation_filepaths = kwargs["segmentation_images"],
                                         detection_filepaths = kwargs["detection_tables"],
                                         psf_filepaths = kwargs["psf_images"])
    
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
        
        load_method_data = loading_methods[method]
        method_data_filename = method_data_files[method]
        estimate_shear = estimation_methods[method]
        
        try:
            
            if load_method_data is not None:
                method_data = load_method_data(method_data_filename)
            else:
                method_data = None
            
            tab, method_mcmc_chains = method( data_stack, method_data )
            
            if not is_in_format(tab,setf):
                raise ValueError("Shear estimation table returned in invalid format for method " + method + ".")
            
            if method_mcmc_chains is not None:
                mcmc_chains = method_mcmc_chains
                
        except Exception as e:
            
            logger.warning(str(e))
            
            # Create an empty estimates table
            if method == "BFD":
                tab = initialize_bfd_moments_table(detections_table)
            else:
                tab = initialise_shear_estimates_table(detections_table)
            
            # Fill it with NaN measurements and 1e99 errors
            tab[setf.ID] = detections_table[detf.ID]
            
            for col in (setf.gal_g1, setf.gal_g2,
                        setf.gal_e1_err, setf.gal_e2_err,):
                tab[col] = np.NaN*np.ones_like(tab[setf.ID])
            
            for col in (setf.gal_g1_err, setf.gal_g2_err):
                tab[col] = 1e99*np.ones_like(tab[setf.ID])
            
        method_shear_estimates[method] = tab
        
    logger.info("Finished estimating shear. Outputting results to " + kwargs["shear_measurements_table"] + ".")
        
    # Output the shear estimates
    output_shear_estimates( method_shear_estimates, kwargs['shear_measurements_product'])
    
    logger.debug("Exiting estimate_shears_from_args")
    
    return
