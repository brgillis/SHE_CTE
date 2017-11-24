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

from astropy.io import fits
from astropy.table import Table

from ElementsKernel.Logging import getLogger
from SHE_CTE_PSFFitting import magic_values as mv
from SHE_PPT import aocs_time_series_product
from SHE_PPT import astrometry_product
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT import mission_time_product
from SHE_PPT import psf_calibration_product
from SHE_PPT.aocs_time_series_product import DpdSheAocsTimeSeriesProduct
from SHE_PPT.astrometry_product import DpdSheAstrometryProduct
from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT import detector as dtc
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             append_hdu, get_allowed_filename,
                             find_file_in_path)
from SHE_PPT.magic_values import extname_label, scale_label, stamp_size_label
from SHE_PPT.mission_time_product import DpdSheMissionTimeProduct
from SHE_PPT.psf_calibration_product import DpdShePSFCalibrationProduct
from SHE_PPT.psf_table_format import initialise_psf_table, tf as pstf
from SHE_PPT.she_image import SHEImage
from SHE_PPT.she_stack import SHEStack
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from SHE_PPT.utility import find_extension
import numpy as np


aocs_time_series_product.init()

astrometry_product.init()

mission_time_product.init()

psf_calibration_product.init()

def fit_psfs(args, dry_run=False):
    """
        Mock run of PSF Fitting.
    """

    logger = getLogger(mv.logger_name)
    
    # Load in the files in turn to make sure there aren't any issues with them.
    
    # Data images
    
    if dry_run:
        dry_label = " dry"
    else:
        dry_label = ""
    
    logger.info("Reading mock"+dry_label+" data images...")
    
    data_images = read_listfile(join(args.workdir,args.data_images))

    sci_hdus = []
    noisemap_hdus = []
    mask_hdus = []
    
    she_images = []
    
    for i, filename in enumerate(data_images):
        
        qualified_filename = join(args.workdir,filename)
        
        data_image_hdulist = fits.open(qualified_filename, mode="readonly", memmap=True)
        num_detectors = len(data_image_hdulist) // 3
        
        sci_hdus.append([])
        noisemap_hdus.append([])
        mask_hdus.append([])
        she_images.append([])
        
        for j in range(num_detectors):
            
            id_string = dtc.get_id_string(j%6+1,j//6+1)
            
            sci_extname = id_string + "." + ppt_mv.sci_tag
            sci_index = find_extension(data_image_hdulist, sci_extname)
            
            sci_hdus[i].append( data_image_hdulist[sci_index] )
            
            noisemap_extname = id_string + "." + ppt_mv.noisemap_tag
            noisemap_index = find_extension(data_image_hdulist, noisemap_extname)
            
            noisemap_hdus[i].append( data_image_hdulist[noisemap_index] )
            
            mask_extname = id_string + "." + ppt_mv.mask_tag
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
    
    logger.info("Reading mock"+dry_label+" detection tables...")
    
    detections_table_filenames = read_listfile(join(args.workdir,args.detections_tables))
    detections_tables = []
    
    for i, filename in enumerate(detections_table_filenames):
        
        detections_tables_hdulist = fits.open(join(args.workdir,filename), mode="readonly", memmap=True)
        num_detectors = len(detections_tables_hdulist)-1
        
        detections_tables.append([])
        
        for j in range(num_detectors):
            
            extname = dtc.get_id_string(j%6+1,j//6+1)+"."+ppt_mv.detections_tag
            table_index = find_extension(detections_tables_hdulist,extname)
            
            detections_tables[i].append( Table.read(detections_tables_hdulist[table_index]) )
            
            if not is_in_format(detections_tables[i][j],detf):
                raise ValueError("Detections table from " + args.detections_tables + " is in invalid format.")
        
    # Astrometry products
    
    logger.info("Reading mock"+dry_label+" astrometry products...")
    
    astrometry_product_filenames = read_listfile(join(args.workdir,args.astrometry_products))
    astrometry_products = []
    
    for i, filename in enumerate(astrometry_product_filenames):
        astrometry_products.append(read_pickled_product(join(args.workdir,filename)))
        if not isinstance(astrometry_products[i], DpdSheAstrometryProduct):
            raise ValueError("Astrometry product from " + join(args.workdir,filename) + " is invalid type.")
        
    # AocsTimeSeries products
    
    logger.info("Reading mock"+dry_label+" aocs time series products...")
    
    aocs_time_series_product_filenames = read_listfile(join(args.workdir,args.aocs_time_series_products))
    aocs_time_series_products = []
    
    for i, filename in enumerate(aocs_time_series_product_filenames):
        aocs_time_series_products.append(read_pickled_product(join(args.workdir,filename)))
        if not isinstance(aocs_time_series_products[i], DpdSheAocsTimeSeriesProduct):
            raise ValueError("AocsTimeSeries product from " + filename + " is invalid type.")
        
    # MissionTime products
    
    logger.info("Reading mock"+dry_label+" mission time products...")
    
    mission_time_product_filenames = read_listfile(join(args.workdir,args.mission_time_products))
    mission_time_products = []
    
    for i, filename in enumerate(mission_time_product_filenames):
        mission_time_products.append(read_pickled_product(join(args.workdir,filename)))
        if not isinstance(mission_time_products[i], DpdSheMissionTimeProduct):
            raise ValueError("MissionTime product from " + filename + " is invalid type.")
        
    # PSFCalibration products
    
    logger.info("Reading mock"+dry_label+" PSF calibration products...")
    
    all_psf_calibration_product_filenames = read_listfile(join(args.workdir,args.psf_calibration_products))
    psf_calibration_product_filenames = all_psf_calibration_product_filenames[0]
    psf_calibration_product_sub_filenames = all_psf_calibration_product_filenames[1]
    psf_calibration_products = []
    
    for i, filename in enumerate(psf_calibration_product_filenames):
        
        psf_calibration_products.append(read_pickled_product(join(args.workdir,filename))
        
        if not isinstance(psf_calibration_products[i], DpdShePSFCalibrationProduct):
            raise ValueError("PSFCalibration product from " + filename + " is invalid type.")
    
    # Set up mock output in the correct format
    
    logger.info("Outputting mock"+dry_label+" PSF images and tables...")
    
    num_exposures = len(data_images)
    psf_image_and_table_filenames = []
    
    if not dry_run:
                
        logger.debug("Searching for mock psf in " + os.environ['ELEMENTS_AUX_PATH'] + ".")
        
        psf_path = find_file_in_path("SHE_CTE_PSFFitting/sample_psf.fits",os.environ['ELEMENTS_AUX_PATH'])
        
        if psf_path is None:
            raise Exception("Cannot find mock psf.")
                
        logger.debug("Found mock psf: " + psf_path)
    
    for i in range(num_exposures):
        
        filename = get_allowed_filename("PSF_DRY",str(i))
        psf_image_and_table_filenames.append(filename)
        
        hdulist = fits.HDUList()
        
        for j in range(num_detectors):
            
            psfc = initialise_psf_table(detector_x = j%6 + 1,
                                        detector_y = j//6 + 1)
    
            if not dry_run:
                
                sample_psf = fits.open(psf_path)[0].data
                bpsf_array = sample_psf
                dpsf_array = sample_psf
                
                for ID in detections_tables[i][j][detf.ID]:
                    psfc.add_row({pstf.ID: ID,
                                  pstf.template: 0,
                                  pstf.stamp_x: np.shape(sample_psf)[1]//2,
                                  pstf.stamp_y: np.shape(sample_psf)[0]//2,
                                  pstf.psf_x: np.shape(sample_psf)[1]//2,
                                  pstf.psf_y: np.shape(sample_psf)[0]//2,
                                  pstf.cal_time: "",
                                  pstf.field_time: ""})
                
            else:
                
                bpsf_array = np.zeros((1,1))
                dpsf_array = np.zeros((1,1))
                
            bulge_psf_header = fits.header.Header(((extname_label,dtc.get_id_string(j%6+1,j//6+1)+"."+ppt_mv.bulge_psf_tag),
                                             (stamp_size_label,np.min(np.shape(bpsf_array))),
                                             (scale_label,0.02)))
            
            bpsf_hdu = fits.ImageHDU(data=bpsf_array,
                                     header=bulge_psf_header)
            hdulist.append(bpsf_hdu)
                
            disk_psf_header = fits.header.Header(((extname_label,dtc.get_id_string(j%6+1,j//6+1)+"."+ppt_mv.disk_psf_tag),
                                             (stamp_size_label,np.min(np.shape(bpsf_array))),
                                             (scale_label,0.02)))
            
            dpsf_hdu = fits.ImageHDU(data=bpsf_array,
                                     header=disk_psf_header)
            hdulist.append(dpsf_hdu)
            
            psfc_hdu = table_to_hdu(psfc)
            hdulist.append(psfc_hdu)
            
        hdulist.writeto(join(args.workdir,filename),clobber=True)
    
        logger.info("Finished outputting mock"+dry_label+" PSF images and tables for exposure " + str(i) + ".")
        
    write_listfile( join(args.workdir,args.psf_images_and_tables), psf_image_and_table_filenames )
    
    logger.info("Finished mock"+dry_label+" psf fitting.")
        
    return
    
    