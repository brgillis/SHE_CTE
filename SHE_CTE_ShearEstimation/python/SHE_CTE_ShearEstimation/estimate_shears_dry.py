""" @file estimate_shears_dry.py

    Created 12 Oct 2017

    Function for performing a dry run of mock shear estimation.

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
from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT.file_io import read_listfile, write_listfile, read_pickled_product, write_pickled_product
from SHE_PPT.psf_table_format import tf as psft
from SHE_PPT.she_image import she_image
from SHE_PPT.table_utility import is_in_format

def estimate_shears_from_args(args):
    """
        Dry-run of shear estimation, creating only dummy files and not expecting anything of input
        aside from the correct format.
    """
    
    # Load in the files in turn to make sure there aren't any issues with them.
    
    # Data images - Read in as SHEStack object
    
    data_images = read_listfile(args.data_images)
    she_stack = SHEStack.read(data_images,
                              data_ext = ppt_mv.sci_tag,
                              mask_ext = ppt_mv.mask_tag,
                              noisemap_ext = ppt_mv.noisemap_tag)
    
    # Detections tables
    
    detection_table_filenames = read_listfile(args.detections_tables)
    detection_tables = []
    
    for i, filename in enumerate(detection_table_filenames):
        
        detections_tables_hdulist = fits.open(filename, mode="readonly", memmap=True)
        num_detectors = len(detections_tables_hdulist)-1
        
        detections_tables.append([])
        
        for j in range(num_detectors):
            
            detections_tables[i].append( Table.read( args.detections_tables, format='fits', hdu=j+1 ) )
            
            if not is_in_format(detections_tables[i][j],detf):
                raise ValueError("Detections table from " + args.detections_tables + " is in invalid format.")
            
    
    
    # Set up mock output in the correct format
    
    def append_hdu(filename, hdu):
        f = fits.open(filename, mode='append')
        try:
            f.append(hdu)
        finally:
            f.close()
    
    for i in range(num_detectors):
            
        bpsf_hdu = fits.ImageHDU(data=np.zeros((1,1)),
                                 header=fits.header.Header(("EXTNAME",str(i)+"."+ppt_mv.bulge_psf_tag)))
        append_hdu( args.psf_images_and_tables, bpsf_hdu)
        
        dpsf_hdu = fits.ImageHDU(data=np.zeros((1,1)),
                                 header=fits.header.Header(("EXTNAME",str(i)+"."+ppt_mv.disk_psf_tag)))
        append_hdu( args.psf_images_and_tables, dpsf_hdu)
        
        psfc_hdu = table_to_hdu(initialise_psf_table(detector=i))
        append_hdu( args.psf_images_and_tables, psfc_hdu)
        
    return
    
    