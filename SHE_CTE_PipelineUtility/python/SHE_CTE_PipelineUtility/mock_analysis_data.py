""" @file fit_psfs_dry.py

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

from os.path import join

from astropy.io import fits
from astropy.table import Table

from SHE_PPT.logging import getLogger
from SHE_CTE_PipelineUtility import magic_values as mv
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT.table_formats.detections import initialise_detections_table
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             append_hdu, get_allowed_filename)
from SHE_PPT.table_formats.galaxy_population import initialise_galaxy_population_table
from SHE_PPT import products
from SHE_PPT.table_formats.shear_estimates import initialise_shear_estimates_table
from SHE_PPT.table_utility import table_to_hdu
import numpy as np

products.aocs_time_series.init()
products.calibration_parameters.init()
products.galaxy_population.init()
products.mission_time.init()
products.mosaic.init()
products.psf_calibration.init()
products.shear_validation_stats.init()

def make_mock_analysis_data(args, dry_run=False):
    """
        Dry-run of generating mock analysis data, creating only dummy files.
    """

    logger = getLogger(__name__)
    
    # Set up mock output in the correct format
    
    num_exposures = args.num_exposures
    num_detectors = args.num_detectors
    
    # Data images
    
    if args.data_images is not None:
    
        logger.info("Generating mock dry data images...")
        
        data_images = []
        
        for i in range(num_exposures):
            
            filename = get_allowed_filename("OBS_DRY",str(i))
            data_images.append(filename)
            
            hdulist = fits.HDUList()
            
            for j in range(num_detectors):
                
                # Science image
                im_hdu = fits.ImageHDU(data=np.zeros((1,1),dtype=np.dtype('>f4')),
                                       header=fits.header.Header((("EXTNAME",str(j)+"."+ppt_mv.sci_tag),)))
                hdulist.append(im_hdu)
                
                # Noise map
                rms_hdu = fits.ImageHDU(data=np.zeros((1,1),dtype=np.dtype('>f4')),
                                        header=fits.header.Header((("EXTNAME",str(j)+"."+ppt_mv.noisemap_tag),)))
                hdulist.append(rms_hdu)
                
                # Mask map
                flg_hdu = fits.ImageHDU(data=np.zeros((1,1),dtype=np.dtype('>i4')),
                                        header=fits.header.Header((("EXTNAME",str(j)+"."+ppt_mv.mask_tag),)))
                hdulist.append(flg_hdu)
                
            hdulist.writeto(join(args.workdir,filename),clobber=True)
        
            logger.info("Finished generating exposure " + str(i) + ".")
            
        write_listfile(join(args.workdir,args.data_images),data_images)
    
    # PSF Calibration products
    
    if args.psf_calibration_products is not None:
    
        logger.info("Generating mock dry psf calibration products...")
        
        psf_calibration_product_filenames = []
        
        for i in range(num_exposures):
            
            # Set up ZM object
            zernike_mode_filename = get_allowed_filename("PSFCAL_ZM_DRY",str(i))
            
            null_hdu = fits.ImageHDU(data=np.zeros((1,1)))
            append_hdu( join(args.workdir,zernike_mode_filename), null_hdu)
            
            # Set up SE object
            surface_error_filename = get_allowed_filename("PSFCAL_SE_DRY",str(i))
            
            null_hdu = fits.ImageHDU(data=np.zeros((1,1)))
            append_hdu( join(args.workdir,surface_error_filename), null_hdu)
            
            # Set up PSF Calibration object
            filename = get_allowed_filename("PSFCAL_DRY",str(i),extension=".bin")
            
            psf_calibration_product = products.psf_calibration.create_psf_calibration_product(timestamp="0",
                                                                     zernike_mode_filename=zernike_mode_filename,
                                                                     surface_error_filename=surface_error_filename)
            
            write_pickled_product(psf_calibration_product,
                                  join(args.workdir,filename))
            
            psf_calibration_product_filenames.append(filename)
            
        write_listfile(join(args.workdir,args.psf_calibration_products),
                       psf_calibration_product_filenames)
    
    # Segmentation images
    
    if args.segmentation_images is not None:
    
        logger.info("Generating mock dry segmentation images...")
        
        mosaic_product_filenames = []
        mosaic_product_sub_filenames = []
        
        for i in range(num_exposures):
            
            # Create the data product
            
            product_filename = get_allowed_filename("SEG_DRY",str(i),extension=".bin")
            data_filename = get_allowed_filename("SEG_DRY_DATA",str(i))
            
            mosaic_product_filenames.append([product_filename])
            mosaic_product_sub_filenames.append([data_filename])
            
            mosaic_product =  products.mosaic.create_mosaic_product(instrument_name="VIS",
                                                    filter="VIS",
                                                    wcs_params=None,
                                                    zeropoint=0,
                                                    data_filename=data_filename,)
            
            write_pickled_product(mosaic_product, product_filename)
            
            # Create its associated fits file
            
            hdulist = fits.HDUList()
            
            for j in range(num_detectors):
            
                seg_hdu = fits.ImageHDU(data=np.zeros((1,1)),
                                        header=fits.header.Header((("EXTNAME",str(j)+"."+ppt_mv.segmentation_tag),)))
                hdulist.append(seg_hdu)
                
            hdulist.writeto(join(args.workdir,data_filename),clobber=True)
        
            logger.info("Finished generating segmentation images for exposure " + str(i) + ".")
        
        write_listfile(join(args.workdir,args.segmentation_images),mosaic_product_filenames)
    
    # Detections tables
    
    if args.detections_tables is not None:
    
        logger.info("Generating mock dry detection tables...")
        
        detections_tables_filenames = []
        
        for i in range(num_exposures):
            
            filename = get_allowed_filename("DTC_DRY",str(i))
            
            hdulist = fits.HDUList()
            
            for j in range(num_detectors):
            
                dtc_hdu = table_to_hdu(initialise_detections_table(detector_x= j%6 + 1,
                                                                   detector_y= j//6 + 1))
                hdulist.append(dtc_hdu)
                
            hdulist.writeto(join(args.workdir,filename),clobber=True)
            
            detections_tables_filenames.append(filename)
        
            logger.info("Finished generating detections for exposure " + str(i) + ".")
        
        write_listfile(join(args.workdir,args.detections_tables),detections_tables_filenames)
    
    # AOCS time series products
    
    if args.aocs_time_series_products is not None:
    
        logger.info("Generating mock dry AOCS time series products...")
        
        aocs_time_series_product_filenames = []
        
        for i in range(num_exposures):
            
            filename = get_allowed_filename("AOCS_TIME_SERIES_DRY",str(i),extension=".bin")
            
            aocs_time_series_product = products.aocs_time_series.create_aocs_time_series_product()
            
            write_pickled_product(aocs_time_series_product,join(args.workdir,filename))
            
            aocs_time_series_product_filenames.append(filename)
            
        write_listfile(join(args.workdir,args.aocs_time_series_products),aocs_time_series_product_filenames)
    
    # Mission time products
    
    if args.mission_time_products is not None:
    
        logger.info("Generating mock dry mission time products...")
        
        mission_time_product_filenames = []
        
        for i in range(num_exposures):
            
            filename = get_allowed_filename("MISSON_TIME_DRY",str(i),extension=".bin")
            
            mission_time_product = products.mission_time.create_mission_time_product()
            
            write_pickled_product(mission_time_product,join(args.workdir,filename))
            
            mission_time_product_filenames.append(filename)
            
        write_listfile(join(args.workdir,args.mission_time_products),mission_time_product_filenames)
    
    # Galaxy population priors table
    
    if args.galaxy_population_priors_table is not None:
    
        logger.info("Generating mock dry galaxy population priors table...")
        
        filename = get_allowed_filename("GALPOP", "0", extension=".fits")
        galaxy_population_prod = products.galaxy_population.create_galaxy_population_product(filename=filename)
        
        write_pickled_product(galaxy_population_prod, join(args.workdir,args.galaxy_population_priors_table))
        
        galaxy_population_priors_table = initialise_galaxy_population_table()
        galaxy_population_priors_table.write(join(args.workdir,filename),
                                             format="fits",overwrite=True)
    
    # Calibration parameters product
    
    if args.calibration_parameters_product is not None:
    
        logger.info("Generating mock dry calibration parameters products...")
        
        ksb_calibration_parameters_filename = get_allowed_filename("KSB_CAL_PARAM_DRY","0")
            
        null_hdu = fits.ImageHDU(data=np.zeros((1,1)))
        append_hdu( join(args.workdir,ksb_calibration_parameters_filename), null_hdu)
            
        lensmc_calibration_parameters_filename = get_allowed_filename("LENSMC_CAL_PARAM_DRY","0")
            
        null_hdu = fits.ImageHDU(data=np.zeros((1,1)))
        append_hdu( join(args.workdir,lensmc_calibration_parameters_filename), null_hdu)
            
        momentsml_calibration_parameters_filename = get_allowed_filename("MOMENTSML_CAL_PARAM_DRY","0")
            
        null_hdu = fits.ImageHDU(data=np.zeros((1,1)))
        append_hdu( join(args.workdir,momentsml_calibration_parameters_filename), null_hdu)
            
        regauss_calibration_parameters_filename = get_allowed_filename("REGAUSS_CAL_PARAM_DRY","0")
            
        null_hdu = fits.ImageHDU(data=np.zeros((1,1)))
        append_hdu( join(args.workdir,regauss_calibration_parameters_filename), null_hdu)
            
        calibration_parameters_product = products.calibration_parameters.create_calibration_parameters_product(KSB_filename=ksb_calibration_parameters_filename,
                                                                               LensMC_filename=lensmc_calibration_parameters_filename,
                                                                               MomentsML_filename=momentsml_calibration_parameters_filename,
                                                                               REGAUSS_filename=regauss_calibration_parameters_filename)
        
        write_pickled_product(calibration_parameters_product,
                              join(args.workdir,args.calibration_parameters_product))
    
    # Shear validation statistics tables
    
    if args.shear_validation_statistics_table is not None:
    
        logger.info("Generating mock dry shear validation statistics tables...")
        
        shear_validation_stats_filename = get_allowed_filename("VAL_STATS", "0")
        
        shear_validation_stats_prod = products.shear_validation_stats.create_shear_validation_stats_product(
                                        shear_validation_stats_filename)
        
        write_pickled_product(shear_validation_stats_prod,
                              join(args.workdir,args.shear_validation_statistics_table))
        
        shear_validation_statistics_table = initialise_shear_estimates_table()
        shear_validation_statistics_table.write(join(args.workdir,shear_validation_stats_filename)
                                                ,format="fits",overwrite=True)
    
    logger.info("Finished generating mock dry data.")
    
    return
    
    