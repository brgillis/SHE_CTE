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
from os.path import join
from astropy.io import fits
from astropy.table import Table

from ElementsKernel.Logging import getLogger

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_PPT import magic_values as ppt_mv

from SHE_PPT.calibration_parameters_product import DpdSheCalibrationParametersProduct
from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, append_hdu)
from SHE_PPT.galaxy_population_table_format import tf as gptf
from SHE_PPT.psf_table_format import tf as pstf
from SHE_PPT.she_image import SHEImage
from SHE_PPT.she_stack import SHEStack
from SHE_PPT import shear_estimates_product as sep
from SHE_PPT.shear_estimates_table_format import initialise_shear_estimates_table
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from SHE_PPT.utility import find_extension

sep.init()

from SHE_PPT import calibration_parameters_product
calibration_parameters_product.init()

def estimate_shears_from_args(args):
    """
        Dry-run of shear estimation, creating only dummy files and not expecting anything of input
        aside from the correct format.
    """

    logger = getLogger(mv.logger_name)
    
    # Load in the files in turn to make sure there aren't any issues with them.
    
    # Data images - Read in as SHEStack object
    
    logger.info("Reading mock dry data images...")
    
    data_images = read_listfile(join(args.workdir,args.data_images))

    sci_hdus = []
    noisemap_hdus = []
    mask_hdus = []
    
    she_images = []
    
    for i, filename in enumerate(data_images):
        
        data_image_hdulist = fits.open(join(args.workdir,filename), mode="readonly", memmap=True)
        num_detectors = len(data_image_hdulist) // 3
        
        sci_hdus.append([])
        noisemap_hdus.append([])
        mask_hdus.append([])
        she_images.append([])
        
        for j in range(num_detectors):
            
            sci_extname = str(j) + "." + ppt_mv.sci_tag
            sci_index = find_extension(data_image_hdulist, sci_extname)
            
            sci_hdus[i].append( data_image_hdulist[sci_index] )
            
            noisemap_extname = str(j) + "." + ppt_mv.noisemap_tag
            noisemap_index = find_extension(data_image_hdulist, noisemap_extname)
            
            noisemap_hdus[i].append( data_image_hdulist[noisemap_index] )
            
            mask_extname = str(j) + "." + ppt_mv.mask_tag
            mask_index = find_extension(data_image_hdulist, mask_extname)
            
            mask_hdus[i].append( data_image_hdulist[mask_index] )
            
            she_images[i].append(SHEImage(data=sci_hdus[i][j].data,
                                          noisemap=noisemap_hdus[i][j].data,
                                          mask=mask_hdus[i][j].data,
                                          header=sci_hdus[i][j].header))
            
    num_exposures = len(data_images)
    she_stacks = []
    for j in range(num_detectors):
        detector_images = []
        for i in range(num_exposures):
            detector_images.append(she_images[i][j])
        she_stacks.append(detector_images)
    
    # Detections tables
    
    logger.info("Reading mock dry detections tables...")
    
    detections_table_filenames = read_listfile(join(args.workdir,args.detections_tables))
    detections_tables = []
    
    for i, filename in enumerate(detections_table_filenames):
        
        detections_tables_hdulist = fits.open(join(args.workdir,filename), mode="readonly", memmap=True)
        num_detectors = len(detections_tables_hdulist)-1
        
        detections_tables.append([])
        
        for j in range(num_detectors):
            
            extname = str(j)+"."+ppt_mv.detections_tag
            table_index = find_extension(detections_tables_hdulist,extname)
            
            detections_tables[i].append( Table.read(detections_tables_hdulist[table_index]) )
            
            if not is_in_format(detections_tables[i][j],detf):
                raise ValueError("Detections table from " + args.detections_tables + " is in invalid format.")
    
    # PSF images and tables
    
    logger.info("Reading mock dry PSF images and tables...")
    
    psf_images_and_table_filenames = read_listfile(join(args.workdir,args.psf_images_and_tables))
    psf_tables = []
    bulge_psf_hdus = []
    disk_psf_hdus = []
    
    for i, filename in enumerate(psf_images_and_table_filenames):
        
        psf_images_and_table_hdulist = fits.open(join(args.workdir,filename), mode="readonly", memmap=True)
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
    
    logger.info("Reading mock dry segmentation images...")
    
    segmentation_filenames = read_listfile(join(args.workdir,args.segmentation_images))
    segmentation_hdus = []
    
    for i, filename in enumerate(segmentation_filenames):
        
        segmentation_hdulist = fits.open(join(args.workdir,filename), mode="readonly", memmap=True)
        num_detectors = len(segmentation_hdulist)
        
        segmentation_hdus.append([])
        
        for j in range(num_detectors):
            
            segmentation_extname = str(j) + "." + ppt_mv.segmentation_tag
            segmentation_index = find_extension(segmentation_hdulist, segmentation_extname)
            
            segmentation_hdus[i].append(segmentation_hdulist[segmentation_index])
            
    # Galaxy population priors
    
    logger.info("Reading mock dry galaxy population priors...")
    
    galaxy_population_priors_table = Table.read(join(args.workdir,args.galaxy_population_priors_table))
            
    if not is_in_format(galaxy_population_priors_table,gptf):
        raise ValueError("Galaxy population priors table from " + join(args.workdir,args.galaxy_population_priors_table) +
                         " is in invalid format.")
        
    # Calibration parameters product
    
    logger.info("Reading mock dry calibration parameters...")
    
    calibration_parameters_product = read_pickled_product(join(args.workdir,args.calibration_parameters_product),
                                                          join(args.workdir,args.calibration_parameters_listfile))
    if not isinstance(calibration_parameters_product, DpdSheCalibrationParametersProduct):
        raise ValueError("CalibrationParameters product from " + join(args.workdir,args.calibration_parameters_product)
                         + " is invalid type.")
    
    # Set up mock output in the correct format
    
    logger.info("Generating mock dry shear estimes...")
    
    shear_estimates_product = sep.create_shear_estimates_product(BFD_filename = get_allowed_filename("DRY_BFD_SHM","0"),
                                                                 KSB_filename = get_allowed_filename("DRY_KSB_SHM","0"),
                                                                 LensMC_filename = get_allowed_filename("DRY_LensMC_SHM","0"),
                                                                 MegaLUT_filename = get_allowed_filename("DRY_MegaLUT_SHM","0"),
                                                                 REGAUSS_filename = get_allowed_filename("DRY_REGAUSS_SHM","0"))
    
    
    for filename in shear_estimates_product.get_all_filenames():
        
        hdulist = fits.HDUList()
        
        for j in range(num_detectors):
            
            shm_hdu = table_to_hdu(initialise_shear_estimates_table(detector=j))
            hdulist.append(shm_hdu)
            
        hdulist.writeto(join(args.workdir,filename),clobber=True)
    
    write_pickled_product(shear_estimates_product,
                          join(args.workdir,args.shear_estimates_product),
                          join(args.workdir,args.shear_estimates_listfile))
    
    logger.info("Finished mock dry shear estimation.")
   
    return
    
    
