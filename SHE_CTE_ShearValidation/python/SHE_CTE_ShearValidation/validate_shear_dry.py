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

from SHE_PPT import magic_values as ppt_mv
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, append_hdu)
from SHE_PPT import shear_estimates_product as sep
from SHE_PPT.shear_estimates_table_format import tf as setf
from SHE_PPT.table_utility import is_in_format

sep.init()

def validate_shear(args):
    """
        Dry-run of shear validation, creating only dummy files and not expecting anything of input
        aside from the correct format.
    """
    
    # Load in the files in turn to make sure there aren't any issues with them.
    
    # Shear estimates product
    
    shear_estimates_product = read_pickled_product(args.shear_estimates_product, args.shear_estimates_listfile)
    
    if not isinstance(shear_estimates_product, DpdSheCalibrationParametersProduct):
        raise ValueError("Shear estimates product from " + args.shear_estimates_product + " is invalid type.")
    
    def find_extension(hdulist,extname):
        for i, hdu in enumerate(hdulist):
            if hdu.header["EXTNAME"]==extname:
                return i
        return None
    
    for filename in shear_estimates_product.get_all_filenames():
    
        shear_estimates_hdulist = fits.open(filename,mode='readonly',memmap=True)
        
        num_detectors = len(shear_estimates_hdulist)-1
        
        for j in range(num_detectors):
            
            table_extname = str(j) + "." + ppt_mv.shear_estimates_tag
            table_index = find_extension(shear_estimates_hdulist, table_extname)
            
            shear_estimates_table = Table.read(filename,format='fits',hdu=table_index)
            
            if not is_in_format(shear_estimates_table,setf):
                raise ValueError("Shear estimates table from " + filename + " is in invalid format.")
            
    # Shear validation statistics
    
    shear_validation_statistics_table = Table.read(args.shear_validation_statistics_table)
            
    if not is_in_format(shear_validation_statistics_table,setf):
        raise ValueError("Shear validation statistics table from " + filename + " is in invalid format.")
    
    
    # Set up mock output in the correct format
        
    for j in range(num_detectors):
        
        shm_hdu = table_to_hdu(initialise_shear_estimates_table(detector=j))
        append_hdu( args.validated_shear_estimates_table, shm_hdu)
        
    return
    
    