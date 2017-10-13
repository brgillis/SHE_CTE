""" @file fit_psfs_dry.py

    Created 12 Oct 2017

    Function for performing a dry run of mock psf fitting.

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
from SHE_PPT.aocs_time_series_product import DpdSheAocsTimeSeriesProduct
from SHE_PPT.astrometry_product import DpdSheAstrometryProduct
from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT.file_io import read_listfile, read_pickled_product, write_pickled_product, append_hdu
from SHE_PPT.mission_time_product import DpdSheMissionTimeProduct
from SHE_PPT.psf_calibration_product import DpdShePSFCalibrationProduct
from SHE_PPT.psf_table_format import tf as psft
from SHE_PPT.she_stack import SHEStack
from SHE_PPT.table_utility import is_in_format

from SHE_PPT import aocs_time_series_product
aocs_time_series_product.init()

from SHE_PPT import astrometry_product
astrometry_product.init()

from SHE_PPT import mission_time_product
mission_time_product.init()

from SHE_PPT import psf_calibration_product
psf_calibration_product.init()

def fit_psfs(args):
    """
        Dry-run of PSF Fitting, creating only dummy files and not expecting anything of input
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
        
    # Astrometry products
    
    astrometry_product_filenames = read_listfile(args.astrometry_products)
    astrometry_products = []
    
    for i, filename in enumerate(astrometry_product_filenames):
        astrometry_products.append(read_pickled_product(filename))
        if not isinstance(astrometry_products[i], DpdSheAstrometryProduct):
            raise ValueError("Astrometry product from " + filename + " is invalid type.")
        
    # AocsTimeSeries products
    
    aocs_time_series_product_filenames = read_listfile(args.aocs_time_series_products)
    aocs_time_series_products = []
    
    for i, filename in enumerate(aocs_time_series_product_filenames):
        aocs_time_series_products.append(read_pickled_product(filename))
        if not isinstance(aocs_time_series_products[i], DpdSheAocsTimeSeriesProduct):
            raise ValueError("AocsTimeSeries product from " + filename + " is invalid type.")
        
    # MissionTime products
    
    mission_time_product_filenames = read_listfile(args.mission_time_products)
    mission_time_products = []
    
    for i, filename in enumerate(mission_time_product_filenames):
        mission_time_products.append(read_pickled_product(filename))
        if not isinstance(mission_time_products[i], DpdSheMissionTimeProduct):
            raise ValueError("MissionTime product from " + filename + " is invalid type.")
        
    # PSFCalibration products
    
    all_psf_calibration_product_filenames = read_listfile(args.psf_calibration_products)
    psf_calibration_product_filenames = all_psf_calibration_product_filenames[0]
    psf_calibration_product_sub_filenames = all_psf_calibration_product_filenames[1]
    psf_calibration_products = []
    
    for i, filename in enumerate(psf_calibration_product_filenames):
        
        psf_calibration_products.append(read_pickled_product(filename, psf_calibration_product_sub_filenames[i]))
        
        if not isinstance(psf_calibration_products[i], DpdShePSFCalibrationProduct):
            raise ValueError("PSFCalibration product from " + filename + " is invalid type.")
    
    # Set up mock output in the correct format
    
    num_exposures = len(data_images)
    psf_image_and_table_filenames = []
    
    for i in range(num_exposures):
        
        filename = get_allowed_filename("PSF_DRY",str(i))
        
        for j in range(num_detectors):
            
            bpsf_hdu = fits.ImageHDU(data=np.zeros((1,1)),
                                     header=fits.header.Header(("EXTNAME",str(j)+"."+ppt_mv.bulge_psf_tag)))
            append_hdu( filename, bpsf_hdu)
            
            dpsf_hdu = fits.ImageHDU(data=np.zeros((1,1)),
                                     header=fits.header.Header(("EXTNAME",str(j)+"."+ppt_mv.disk_psf_tag)))
            append_hdu( filename, dpsf_hdu)
            
            psfc_hdu = table_to_hdu(initialise_psf_table(detector=j))
            append_hdu( filename, psfc_hdu)
        
    write_listfile( args.psf_images_and_tables, psf_image_and_table_filenames )
        
    return
    
    