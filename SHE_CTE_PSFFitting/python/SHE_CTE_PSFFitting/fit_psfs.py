""" @file fit_psf.py

    Created 12 Oct 2017

    Function for performing a dry run of mock psf fitting.
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

import os
from os.path import join
from copy import deepcopy

from astropy.io import fits
from astropy.table import Table

from SHE_CTE_PSFFitting import magic_values as mv
from SHE_PPT.logging import getLogger
from SHE_PPT import products
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, find_file_in_path)
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.psf import initialise_psf_table, tf as pstf
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.table_utility import is_in_format, table_to_hdu
import numpy as np

test_mode = True

def fit_psfs(args, dry_run=False):
    """
        Mock run of PSF Fitting.
    """

    logger = getLogger(mv.logger_name)
    
    # Load in the files in turn to make sure there aren't any issues with them.
    
    # Data images
    
    if dry_run:
        dry_label = "_dry"
    else:
        dry_label = ""
    
    logger.info("Reading mock"+dry_label+" data images and detections tables...")
    
    frame_stack = SHEFrameStack.read(exposure_listfile_filename=args.data_images,
                                     detections_listfile_filename=args.detections_tables,
                                     workdir=args.workdir)
        
    # AocsTimeSeries products
    
    if args.aocs_time_series_products is not None:
    
        logger.info("Reading mock"+dry_label+" aocs time series products...")
        
        aocs_time_series_product_filenames = read_listfile(join(args.workdir,args.aocs_time_series_products))
        aocs_time_series_products = []
        
        for i, filename in enumerate(aocs_time_series_product_filenames):
            aocs_time_series_products.append(read_pickled_product(join(args.workdir,filename)))
            if not isinstance(aocs_time_series_products[i], products.aocs_time_series.DpdSheAocsTimeSeriesProduct):
                raise ValueError("AocsTimeSeries product from " + filename + " is invalid type.")
            
    else:
        
        aocs_time_series_products = None
        
    # PSFCalibration products
    
    if args.psf_calibration_product is not None:
    
        logger.info("Reading mock"+dry_label+" PSF calibration products...")
        
        psf_calibration_product = read_pickled_product(join(args.workdir,args.psf_calibration_product))
        
        if not isinstance(psf_calibration_product, products.psf_calibration.DpdShePSFCalibrationProduct):
            raise ValueError("PSFCalibration product from " + filename + " is invalid type.")
        
    else:
        
        psf_calibration_product = None
    
    # Set up mock output in the correct format
    
    logger.info("Outputting mock"+dry_label+" PSF images and tables...")
    
    num_exposures = len(frame_stack.exposures)
    psf_image_and_table_filenames = []
    
    if not dry_run:
                
        logger.debug("Searching for mock psf in " + os.environ['ELEMENTS_AUX_PATH'] + ".")
        
        psf_path = find_file_in_path("SHE_CTE_PSFFitting/sample_psf.fits",os.environ['ELEMENTS_AUX_PATH'])
        
        if psf_path is None:
            raise Exception("Cannot find mock psf.")
                
        logger.debug("Found mock psf: " + psf_path)
        
        mock_psf_data = fits.open(psf_path)[0].data
    
    # Set up for each exposure
    filenames = []
    psf_tables = []
    for x in range(num_exposures):
        filenames.append(get_allowed_filename("PSF" + dry_label,str(x)))
        
        hdulist = fits.HDUList([fits.PrimaryHDU()]) # Start with an empty primary HDU
            
        psfc = initialise_psf_table()
            
        # Add a line to the table for each galaxy - Will fill in later
        for row in frame_stack.detections_catalogue:
            psfc.add_row({pstf.ID: row[detf.ID],
                          pstf.template: -1,
                          pstf.bulge_index: -1,
                          pstf.disk_index: -1})
        
        # Add the table to the HDU list
        psfc_hdu = table_to_hdu(psfc)
        hdulist.append(psfc_hdu)
        
        psf_tables.append(psfc) # Keep a copy of the table
        psf_tables.set_index(detf.ID) # Allow it to be indexed by galaxy ID
        
        # Write out the table
        hdulist.writeto(join(args.workdir,filename),clobber=True)
        

        
    # For each galaxy, add a bulge and disk PSF image if it's in the frame
    bulge_hdu = 2
    first_galaxy = True
    for row in frame_stack.detections_catalogue:
        
        gal_id = row[detf.ID]
        
        gal_stamp_stack = frame_stack.extract_galaxy_stamp(gal_id)
        
        if gal_stamp_stack.is_empty():
            continue
        
        for x in range(num_exposures):
        
            if first_galaxy or not test_mode:
                # Get mock data for bulge and disk psfs
                bpsf_array = deepcopy(mock_psf_data)
                dpsf_array = deepcopy(mock_psf_data)
                
                # Add the bulge and disk HDUs
                bulge_psf_header = fits.header.Header(((ppt_mv.extname_label,str(gal_id)+"."+ppt_mv.bulge_psf_tag),
                                                       (ppt_mv.stamp_size_label,np.min(np.shape(bpsf_array))),
                                                       (ppt_mv.scale_label,0.02)))
                
                bpsf_hdu = fits.ImageHDU(data=bpsf_array,
                                         header=bulge_psf_header)
                    
                disk_psf_header = fits.header.Header(((ppt_mv.extname_label,str(gal_id)+"."+ppt_mv.disk_psf_tag),
                                                 (ppt_mv.stamp_size_label,np.min(np.shape(dpsf_array))),
                                                 (ppt_mv.scale_label,0.02)))
                
                dpsf_hdu = fits.ImageHDU(data=dpsf_array,
                                         header=disk_psf_header)
                
                # Append these to the proper file
                
                f = fits.open(filenames[x],mode='append')
                
                f[x].append(bpsf_hdu)
                f[x].append(dpsf_hdu)
                
                f.flush()
                f.close()
            
                # Cleanup
                del bpsf_array, dpsf_array, bpsf_hdu, dpsf_hdu, bulge_psf_header, disk_psf_header
                
                first_galaxy = False
            
            # Update the table
            output_row = psf_tables[x].loc[gal_id]
            output_row[pstf.bulge_index] = bulge_hdu
            output_row[pstf.disk_index] = bulge_hdu+1
        
        if not test_mode:
            bulge_hdu += 2
        
    # Go back and update the tables with proper values
    for x in range(num_exposures):
        
        psf_table = psf_tables[x]
        
        f = fits.open(filenames[x],memmap=True,mode='update')
        out_table = f[1].data
        
        # Update each row
        for i in range(len(psf_table)):
            out_table[pstf.bulge_index] = psf_table[i][pstf.bulge_index]
            out_table[pstf.disk_index] = psf_table[i][pstf.disk_index]
            
        f.flush()
        f.close()

    logger.info("Finished outputting mock"+dry_label+" PSF images and tables for exposure " + str(i) + ".")
        
    # Write listfile of objects
    write_listfile( join(args.workdir,args.psf_images_and_tables), psf_image_and_table_filenames )
    
    logger.info("Finished mock"+dry_label+" psf fitting.")
        
    return
    
    