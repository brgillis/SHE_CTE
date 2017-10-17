""" @file validate_shear_dry.py

    Created 12 Oct 2017

    Function for performing a dry run of shear validation.

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
from astropy.io import fits
from astropy.table import Table

from ElementsKernel.Logging import getLogger

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_PPT import magic_values as ppt_mv

from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, append_hdu)
from SHE_PPT.shear_estimates_product import DpdShearEstimatesProduct
from SHE_PPT.shear_estimates_table_format import tf as setf, initialise_shear_estimates_table
from SHE_PPT.calibration_parameters_product import DpdSheCalibrationParametersProduct
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from SHE_PPT.utility import find_extension

from SHE_PPT import shear_estimates_product
shear_estimates_product.init()

from SHE_PPT import calibration_parameters_product
calibration_parameters_product.init()

def validate_shear(args):
    """
        Dry-run of shear validation, creating only dummy files and not expecting anything of input
        aside from the correct format.
    """

    logger = getLogger(mv.logger_name)
    
    # Load in the files in turn to make sure there aren't any issues with them.
    
    # Shear estimates product
    
    logger.info("Reading mock dry shear estimates product...")
    
    shear_estimates_product = read_pickled_product(join(args.workdir,args.shear_estimates_product),
                                                        join(args.workdir,args.shear_estimates_listfile))

    if not isinstance(shear_estimates_product, DpdShearEstimatesProduct):
        raise ValueError("Shear estimates product from " + join(args.workdir,args.shear_estimates_product)
                          + " is invalid type.")
    
    for filename in shear_estimates_product.get_all_filenames():
    
        shear_estimates_hdulist = fits.open(join(args.workdir,filename),mode='readonly',memmap=True)
        
        num_detectors = len(shear_estimates_hdulist)-1
        
        for j in range(num_detectors):
            
            table_extname = str(j) + "." + ppt_mv.shear_estimates_tag
            table_index = find_extension(shear_estimates_hdulist, table_extname)
            
            shear_estimates_table = Table.read(join(args.workdir,filename),format='fits',hdu=table_index)
            
            if not is_in_format(shear_estimates_table,setf):
                raise ValueError("Shear estimates table from " + join(args.workdir,filename) + " is in invalid format.")
            
    # Shear validation statistics
    
    logger.info("Reading mock dry shear validation statistics...")
    
    shear_validation_statistics_table = Table.read(join(args.workdir,args.shear_validation_statistics_table))
            
    if not is_in_format(shear_validation_statistics_table,setf):
        raise ValueError("Shear validation statistics table from " + join(args.workdir,filename) + " is in invalid format.")
    
    
    # Set up mock output in the correct format
    
    logger.info("Generating mock dry validated shear estimates...")
        
    for j in range(num_detectors):
        
        shm_hdu = table_to_hdu(initialise_shear_estimates_table(detector=j))
        append_hdu( join(args.workdir,args.validated_shear_estimates_table), shm_hdu)
    
    logger.info("Finished mock dry shear validation.")
        
    return
    
    
