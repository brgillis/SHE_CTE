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

from SHE_PPT.calibration_parameters_product import DpdSheCalibrationParametersProduct
from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, append_hdu)
from SHE_PPT.galaxy_population_table_format import tf as gptf
from SHE_PPT.psf_table_format import tf as psft
from SHE_PPT.she_image import she_image
from SHE_PPT import shear_estimates_product as sep
from SHE_PPT.shear_estimates_table_format import tf as setf 
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import find_extension

sep.init()

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
        
        hdulist = fits.open(args.detections_tables,mode="readonly",memmap=True)
        
        for j in range(num_detectors):
            
            extname = str(j)+"."+ppt_mv.detections_tag
            table_index = find_extension(hdulist,extname)
            
            detections_tables[i].append( Table.read(hdulist(table_index)) )
            
            if not is_in_format(detections_tables[i][j],detf):
                raise ValueError("Detections table from " + args.detections_tables + " is in invalid format.")
    
    # PSF images and tables
    
    psf_images_and_table_filenames = read_listfile(args.psf_images_and_tables)
    psf_tables = []
    bulge_psf_hdus = []
    disk_psf_hdus = []
    
    for i, filename in enumerate(psf_images_and_table_filenames):
        
        psf_images_and_table_hdulist = fits.open(filename, mode="readonly", memmap=True)
        num_detectors = (len(psf_images_and_table_hdulist)-1) // 3
        
        psf_tables.append([])
        bulge_psf_hdus.append([])
        disk_psf_hdus.append([])
        
        for j in range(num_detectors):
            
            table_extname = str(j) + "." + ppt_mv.psf_cat_tag
            table_index = find_extension(psf_images_and_table_hdulist, table_extname)
            
            psf_tables[i].append( Table.read( psf_images_and_table_hdulist[table_index] ) )
            
            if not is_in_format(psf_tables[i][j],pstf):
                raise ValueError("PSF table from " + filename + " is in invalid format.")
            
            bulge_extname = str(j) + "." + ppt_mv.bulge_psf_tag
            bulge_index = find_extension(psf_images_and_table_hdulist, bulge_extname)
            
            bulge_psf_hdus[i].append(psf_images_and_table_hdulist[bulge_index])
            
            disk_extname = str(j) + "." + ppt_mv.disk_psf_tag
            disk_index = find_extension(psf_images_and_table_hdulist, disk_extname)
            
            disk_psf_hdus[i].append(psf_images_and_table_hdulist[disk_index])
    
    # Segmentation images
    
    segmentation_filenames = read_listfile(args.segmentation_images)
    segmentation_hdus = []
    
    for i, filename in enumerate(segmentation_filenames):
        
        segmentation_hdulist = fits.open(filename, mode="readonly", memmap=True)
        num_detectors = len(psf_images_and_table_hdulist)-1
        
        segmentation_hdus.append([])
        
        for j in range(num_detectors):
            
            segmentation_extname = str(j) + "." + ppt_mv.segmentation_psf_tag
            segmentation_index = find_extension(segmentation_hdulist, segmentation_extname)
            
            segmentation_hdus[i].append(segmentation_hdulist[segmentation_index])
            
    # Galaxy population priors
    galaxy_population_priors_table = Table.read(args.galaxy_population_priors_table)
            
    if not is_in_format(galaxy_population_priors_table,gptf):
        raise ValueError("Galaxy population priors table from " + args.galaxy_population_priors_table +
                         " is in invalid format.")
        
    # Calibration parameters product
    
    calibration_parameters_product = read_pickled_product(args.calibration_parameters_products,
                                                          args.calibration_parameters_listfile)
    if not isinstance(calibration_parameters_product, DpdSheCalibrationParametersProduct):
        raise ValueError("CalibrationParameters product from " + args.calibration_parameters_product + " is invalid type.")
    
    # Set up mock output in the correct format
    shear_estimates_product = sep.create_shear_estimates_product(BFD_filename = get_allowed_filename("DRY_BFD_SHM","0"),
                                                                 KSB_filename = get_allowed_filename("DRY_KSB_SHM","0"),
                                                                 LensMC_filename = get_allowed_filename("DRY_LensMC_SHM","0"),
                                                                 MegaLUT_filename = get_allowed_filename("DRY_MegaLUT_SHM","0"),
                                                                 REGAUSS_filename = get_allowed_filename("DRY_REGAUSS_SHM","0"))
    
    
    for filename in shear_estimates_product.get_all_filenames():
        shear_estimates_table = initialise
        
        for j in range(num_detectors):
            
            shm_hdu = table_to_hdu(initialise_shear_estimates_table(detector=j))
            append_hdu( filename, shm_hdu)
    
    write_pickled_product(shear_estimates_product, args.shear_estimates_product, args.shear_estimates_listfile)
    
    
   
    return
    
    