""" @file validate_shear.py

    Created 12 Oct 2017

    Function for performing shear validation.

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

import math
from os.path import join

from astropy.io import fits
from astropy.table import Table

from ElementsKernel.Logging import getLogger
from SHE_CTE_ShearValidation import magic_values as mv
from SHE_PPT import calibration_parameters_product
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT import shear_estimates_product
from SHE_PPT.calibration_parameters_product import DpdSheCalibrationParametersProduct
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, append_hdu)
from SHE_PPT.shear_estimates_product import DpdShearEstimatesProduct
from SHE_PPT.shear_estimates_table_format import tf as setf, initialise_shear_estimates_table
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from SHE_PPT.utility import find_extension, get_detector
import numpy as np


shear_estimates_product.init()

calibration_parameters_product.init()

def validate_shear_estimates(shear_estimates_table, shear_validation_statistics_table):
    """
        Stub for validating shear estimates against validation statistics. Presently
        sets everything to "pass".
    """

    shear_estimates_table[setf.m.validated] = 1
    
    return

def inv_var_stack( a, a_err ):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering inv_var_stack")
    
    a_inv_var = 1/a_err**2
    
    inv_a_inv_var_sum = 1./a_inv_var.sum()
    
    a_m = (a*a_inv_var).sum()*inv_a_inv_var_sum
    
    a_m_err = sqrt(inv_a_inv_var_sum)
    
    return a_m, a_m_err
        
    logger.debug("Exiting inv_var_stack")

def combine_shear_estimates(detector_estimates, shape_noise_var=0.06):
    """
        Combine all shear estimates for each galaxy measurement, based on inverse-variance
        weighting.
    """
    
    detector = None
    
    # Check that all estimates are for the same detector
    for method in detector_estimates:
        if detector is None:
            detector = get_detector(detector_estimates[method])
        else:
            if detector != get_detector(detector_estimates[method]):
                raise ValueError("Not all estimates tables are for the same detector.")
    
    # Use only methods which pass validation
    validated_methods = []
        
    for method in detector_estimates:
        
        # Skip if the table doesn't pass validation
        if detector_estimates[method].meta[setf.m.validated] > 0:
            validated_methods.append(method)
            
    # Get a set of all object IDs
    IDs = set([])
    for method in validated_methods:
        IDs = set.union(IDs,detector_estimates[method][detf.ID])
            
        # Set up ID column as index
        detector_estimates[method].add_index(detf.ID)
        
    # Initialise combined table
    combined_shear_estimates = initialise_shear_estimates_table(detector=detector,
                                                                optional_columns=[setf.g1_err,setf.g2_err])
    
    for ID in IDs:
        
        g1s = []
        e1_errs = []
        g2s = []
        e2_errs = []
        
        for method in detector_estimates:
            
            tab = detector_estimates[method]
            
            # Skip if the table doesn't pass validation
            if not tab.meta[setf.m.validated] > 0:
                continue
            
            try:
                row = tab.loc[ID]
            except KeyError:
                continue
            
            m_g1 = row[setf.g1]
            m_g2 = row[setf.g2]

            if setf.e1_err in tab.columns:
                m_e1_err = row[setf.e1_err]
                m_e2_err = row[setf.e2_err]
            else:
                m_e1_var = row[setf.g1_err]**2 - shape_noise_var
                m_e2_var = row[setf.g2_err]**2 - shape_noise_var
                
                if m_e1_var < 0 or m_e2_var < 0:
                    # Bad error estimates
                    continue
                
            if m_e1_err < 1e99 and not math.isnan(m_g1):
                g1s.append(m_g1)
                e1_errs.append(m_e1_err)
                
            if m_g2_err < 1e99 and not math.isnan(m_g2):
                g2s.append(m_g2)
                e2_errs.append(m_e2_err)
                
        if len(g1s)>0 and len(g2s)>0:
                
            g1s = np.array(g1s)
            e1_errs = np.array(e1_errs)
            g2s = np.array(g2s)
            e2_errs = np.array(e2_errs)
        
            g1, e1_err = inv_var_stack(g1s, e1_errs)
            g2, e2_err = inv_var_stack(g2s, e2_errs)
            
            g1_err = math.sqrt(e1_err**2 + shape_noise_var)
            g2_err = math.sqrt(e2_err**2 + shape_noise_var)
            
            combined_shear_estimates.add_row({setf.ID: ID,
                                              setf.g1: g1,
                                              setf.g1_err: g1_err,
                                              setf.g2: g2,
                                              setf.g2_err: g2_err,
                                              })
        
        else:
            
            combined_shear_estimates.add_row({setf.ID: ID,
                                              setf.g1: np.nan,
                                              setf.g1_err: 1e99,
                                              setf.g2: np.nan,
                                              setf.g2_err: 1e99,
                                              })
        
    return combined_shear_estimates

def validate_shear(args, dry_run=False):
    """
        Main function for shear validation.
    """

    logger = getLogger(mv.logger_name)
    
    if dry_run:
        dry_label = " dry"
    else:
        dry_label = ""
    
    # Load in the files in turn to make sure there aren't any issues with them.
    
    # Shear estimates product
    
    logger.info("Reading"+dry_label+" shear estimates product...")
    
    shear_estimates_product = read_pickled_product(join(args.workdir,args.shear_estimates_product),
                                                   join(args.workdir,args.shear_estimates_listfile))

    if not isinstance(shear_estimates_product, DpdShearEstimatesProduct):
        raise ValueError("Shear estimates product from " + join(args.workdir,args.shear_estimates_product)
                          + " is invalid type.")
        
    shear_estimates_tables = {"KSB":[],
                              "REGAUSS":[],
                              "MegaLUT":[],
                              "LensMC":[],
                              "BFD":[]}
    
    for method in shear_estimates_tables:
        
        filename = shear_estimates_product.get_method_filename(method)
    
        shear_estimates_hdulist = fits.open(join(args.workdir,filename),mode='readonly',memmap=True)
        
        num_detectors = len(shear_estimates_hdulist)-1
        
        for j in range(num_detectors):
            
            table_extname = str(j) + "." + ppt_mv.shear_estimates_tag
            table_index = find_extension(shear_estimates_hdulist, table_extname)
            
            shear_estimates_table = Table.read(join(args.workdir,filename),format='fits',hdu=table_index)
            
            if not is_in_format(shear_estimates_table,setf):
                raise ValueError("Shear estimates table from " + join(args.workdir,filename) + " is in invalid format.")
            
            shear_estimates_tables[method].append(shear_estimates_table)
            
    # Shear validation statistics
    
    logger.info("Reading"+dry_label+" shear validation statistics...")
    
    shear_validation_statistics_table = Table.read(join(args.workdir,args.shear_validation_statistics_table))
            
    if not is_in_format(shear_validation_statistics_table,setf):
        raise ValueError("Shear validation statistics table from " + join(args.workdir,filename) + " is in invalid format.")
    
    combined_shear_estimates = []
    
    if not dry_run:
        
        for j in range(num_detectors):
            
            detector_estimates = {}
            
            for method in shear_estimates_tables:
                
                validate_shear_estimates(shear_estimates_tables[method][j], shear_validation_statistics_table)
                
                detector_estimates[method] = shear_estimates_tables[method][j]

            combined_shear_estimates.append(combine_shear_estimates(detector_estimates))
                
    else:
        
        for j in range(num_detectors):
        
            combined_shear_estimates.append(initialise_shear_estimates_table(detector=j))
            
    
    # Set up mock output in the correct format
    
    logger.info("Generating"+dry_label+" validated shear estimates...")
        
    for j in range(num_detectors):
        
        shm_hdu = table_to_hdu(combined_shear_estimates[j])
        append_hdu( join(args.workdir,args.validated_shear_estimates_table), shm_hdu)
    
    logger.info("Finished"+dry_label+" shear validation.")
        
    return
    
    
