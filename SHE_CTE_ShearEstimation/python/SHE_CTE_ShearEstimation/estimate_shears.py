""" @file estimate_shears.py

    Created 27 Mar 2017

    Primary execution loop for measuring galaxy shapes from an image file.
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
import astropy.table
import galsim

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_CTE_ShearEstimation.galsim_estimate_shear import KSB_estimate_shear, REGAUSS_estimate_shear
from SHE_CTE_ShearEstimation.bfd_measure_moments import bfd_measure_moments, bfd_load_method_data
from SHE_MomentsML.estimate_shear import estimate_shear as MomentsML_estimate_shear

from SHE_PPT.logging import getLogger
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT import detector as dtc
from SHE_PPT.file_io import (read_listfile, read_pickled_product,
                             write_pickled_product, get_allowed_filename)
from SHE_PPT.table_formats.galaxy_population import tf as gptf
from SHE_PPT import products
from SHE_PPT.table_formats.psf import tf as pstf
from SHE_PPT.she_image import SHEImage
from SHE_PPT.she_image_data import SHEImageData
from SHE_PPT.she_stack import SHEStack
from SHE_PPT.table_formats.shear_estimates import initialise_shear_estimates_table, tf as setf
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from SHE_PPT.utility import find_extension
import numpy as np

products.calibration_parameters.init()
products.calibrated_frame.init()
products.detections.init()
products.galaxy_population.init()
products.shear_estimates.init()
products.mosaic.init()
products.psf_image.init()


loading_methods = {"KSB":None,
                   "REGAUSS":None,
                   "MomentsML":None,
                   "LensMC":None,
                   "BFD":bfd_load_method_data}

estimation_methods = {"KSB":KSB_estimate_shear,
                      "REGAUSS":REGAUSS_estimate_shear,
                      "MomentsML":MomentsML_estimate_shear,
                      "LensMC":None,
                      "BFD":bfd_measure_moments}

def estimate_shears_from_args(args, dry_run=False):
    """
    @brief
        Perform shear estimation, given arguments from the command-line.
    
    @param kwargs <dict>
    
    @return None
    """

    logger = getLogger(mv.logger_name)
    
    logger.debug("Entering estimate_shears_from_args")
    
    # Load in the files in turn
    
    # Data images - Read in as SHEStack object
    
    if dry_run:
        dry_label = " mock dry"
    else:
        dry_label = ""
    
    logger.info("Reading "+dry_label+"data images...")
    
    data_image_products = read_listfile(join(args.workdir,args.data_images))

    sci_hdus = []
    noisemap_hdus = []
    mask_hdus = []
    
    she_images = []

    for i, prod_filename in enumerate(data_image_products):
        
        qualified_prod_filename = join(args.workdir,prod_filename)
        
        data_image_prod = read_pickled_product(qualified_prod_filename)
        if not isinstance(data_image_prod, products.calibrated_frame.DpdSheCalibratedFrameProduct):
            raise ValueError("Data image product from " + qualified_prod_filename + " is invalid type.")
        
        qualified_filename = join(args.workdir,data_image_prod.get_filename())
        
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

    # Detections tables
    
    logger.info("Reading "+dry_label+"detections tables...")
    
    detections_table_prod_filenames = read_listfile(join(args.workdir,args.detections_tables))
    detections_tables = []
    
    for i, prod_filename in enumerate(detections_table_prod_filenames):
        
        detections_table_prod = read_pickled_product(join(args.workdir,prod_filename))
        if not isinstance(detections_table_prod, products.detections.DpdSheDetectionsProduct):
            raise ValueError("Detections product from " + prod_filename + " is invalid type.")
        
        qualified_filename = join(args.workdir,detections_table_prod.get_filename())
        
        detections_tables_hdulist = fits.open(join(args.workdir,qualified_filename), mode="readonly", memmap=True)
        
        detections_tables.append([])
        
        for j in range(num_detectors):
            
            extname = dtc.get_id_string(j%6+1,j//6+1)+"."+ppt_mv.detections_tag
            table_index = find_extension(detections_tables_hdulist,extname)
            
            detections_tables[i].append( Table.read(detections_tables_hdulist[table_index]) )
            
            if not is_in_format(detections_tables[i][j],detf):
                raise ValueError("Detections table from " + args.detections_tables + " is in invalid format.")
    
    # PSF images and tables
    
    logger.info("Reading "+dry_label+"PSF images and tables...")
    
    psf_images_and_table_product_filenames = read_listfile(join(args.workdir,args.psf_images_and_tables))

    she_image_datas = []
    
    for i, prod_filename in enumerate(psf_images_and_table_product_filenames):
        
        qualified_prod_filename = join(args.workdir,prod_filename)
        
        psf_image_prod = read_pickled_product(qualified_prod_filename)
        if not isinstance(psf_image_prod, products.psf_image.DpdShePSFImageProduct):
            raise ValueError("PSF image product from " + qualified_prod_filename + " is invalid type.")
        
        qualified_filename = join(args.workdir,psf_image_prod.get_filename())
        
        psf_images_and_table_hdulist = fits.open(qualified_filename, mode="readonly", memmap=True)
        
        she_image_datas.append([])
        
        for j in range(num_detectors):
            
            id_string = dtc.get_id_string(j%6+1,j//6+1)
            
            table_extname = id_string + "." + ppt_mv.psf_cat_tag
            table_index = find_extension(psf_images_and_table_hdulist, table_extname)
            
            psf_table = Table.read( psf_images_and_table_hdulist[table_index] )
            
            if not is_in_format(psf_table,pstf):
                raise ValueError("PSF table from " + qualified_filename + " is in invalid format.")
            
            bulge_extname = id_string + "." + ppt_mv.bulge_psf_tag
            bulge_index = find_extension(psf_images_and_table_hdulist, bulge_extname)
            
            bulge_psf_image = SHEImage(data=psf_images_and_table_hdulist[bulge_index].data,
                                       header=psf_images_and_table_hdulist[bulge_index].header)
            
            disk_extname = id_string + "." + ppt_mv.disk_psf_tag
            disk_index = find_extension(psf_images_and_table_hdulist, disk_extname)
            
            disk_psf_image = SHEImage(data=psf_images_and_table_hdulist[disk_index].data,
                                      header=psf_images_and_table_hdulist[disk_index].header)
            
            she_image_datas[i].append(SHEImageData(science_image=she_images[i][j],
                                                detections_table=detections_tables[i][j],
                                                bpsf_image=bulge_psf_image,
                                                dpsf_image=disk_psf_image,
                                                psf_table=psf_table))

    num_exposures = len(she_image_datas)
    data_stacks = []
    for j in range(num_detectors):
        detector_image_datas = []
        for i in range(num_exposures):
            detector_image_datas.append(she_image_datas[i][j])
        data_stacks.append(SHEStack(detector_image_datas))


    # Segmentation images
    
    logger.info("Reading "+dry_label+"segmentation images...")
    
    segmentation_product_filenames = read_listfile(join(args.workdir,args.segmentation_images))
    
    segmentation_hdus = []
    
    for i, filename in enumerate(segmentation_product_filenames):
        
        mosaic_prod = read_pickled_product(join(args.workdir,filename))
        
        if not isinstance(mosaic_prod, products.mosaic.DpdMerMosaicProduct):
            raise ValueError("Mosaic product from " + filename + " is invalid type.")
        
        mosaic_data_filename = mosaic_prod.get_data_filename()
        
        segmentation_hdulist = fits.open(join(args.workdir,mosaic_data_filename), mode="readonly", memmap=True)
        
        segmentation_hdus.append([])
        
        for j in range(num_detectors):
            
            segmentation_extname = dtc.get_id_string(j%6+1,j//6+1) + "." + ppt_mv.segmentation_tag
            segmentation_index = find_extension(segmentation_hdulist, segmentation_extname)
            
            segmentation_hdus[i].append(segmentation_hdulist[segmentation_index])
            
    # Galaxy population priors
    
    logger.info("Reading "+dry_label+"galaxy population priors...")
    
    filename = join(args.workdir,args.galaxy_population_priors_table)
    galaxy_population_priors_prod = read_pickled_product(filename)
        
    if not isinstance(galaxy_population_priors_prod, products.galaxy_population.DpdSheGalaxyPopulationProduct):
        raise ValueError("Galaxy population product from " + filename + " is invalid type.")
    
    qualified_galpop_filename = join(args.workdir,galaxy_population_priors_prod.get_filename())
    
    galaxy_population_priors_table = Table.read(qualified_galpop_filename)
            
    if not is_in_format(galaxy_population_priors_table,gptf):
        raise ValueError("Galaxy population priors table from " + join(args.workdir,args.galaxy_population_priors_table) +
                         " is in invalid format.")
        
    # Calibration parameters product
    
    logger.info("Reading "+dry_label+"calibration parameters...")
    
    calibration_parameters_prod = read_pickled_product(join(args.workdir,args.calibration_parameters_product))
    if not isinstance(calibration_parameters_prod, products.calibration_parameters.DpdSheCalibrationParametersProduct):
        raise ValueError("CalibrationParameters product from " + join(args.workdir,args.calibration_parameters_product)
                         + " is invalid type.")
    
    # Set up output
    
    logger.info("Generating shear estimates product...")
    
    shear_estimates_prod = products.shear_estimates.create_shear_estimates_product(
                                BFD_filename = get_allowed_filename("BFD_SHM","0"),
                                KSB_filename = get_allowed_filename("KSB_SHM","0"),
                                LensMC_filename = get_allowed_filename("LensMC_SHM","0"),
                                MomentsML_filename = get_allowed_filename("MomentsML_SHM","0"),
                                REGAUSS_filename = get_allowed_filename("REGAUSS_SHM","0"))
        
    if not dry_run:
        
        # Load the P(e) table if available
        shape_noise_var = 0.06
        
        # TODO - fill in loading of P(e)
        
        method_shear_estimates = {}
        
        if len(args.methods)==0:
            methods = estimation_methods.keys()
        else:
            methods = args.methods

        for method in methods:

            load_method_data = loading_methods[method]
            
            method_data_filename = calibration_parameters_prod.get_method_filename(method)

            shear_estimates_filename = shear_estimates_prod.get_method_filename(method)
            
            estimate_shear = estimation_methods[method]
            
            hdulist = fits.HDUList()
            
            try:

                if load_method_data is not None:
                    method_data = load_method_data(method_data_filename)
                    
                else:
                    method_data = None
                    
                for j in range(num_detectors):
                    
                    data_stack = data_stacks[j]

                    shear_estimates_table = estimate_shear( data_stack, method_data=method_data )

                    if not is_in_format(shear_estimates_table,setf):
                        raise ValueError("Shear estimation table returned in invalid format for method " + method + ".")
 
                    hdulist.append(table_to_hdu(shear_estimates_table))
                    
            except Exception as e:
                
                logger.warning(str(e))
            
                hdulist = fits.HDUList()
                    
                for j in range(num_detectors):
                
                    # Create an empty estimates table
                    shear_estimates_table = initialise_shear_estimates_table(detections_tables[i][j])
                    
                    for r in range(len(detections_tables[i][j][detf.ID])):
                        
                        # Fill it with NaN measurements and 1e99 errors
                        
                        shear_estimates_table.add_row({setf.ID:detections_tables[i][j][detf.ID][r],
                                     setf.g1:np.NaN,
                                     setf.g2:np.NaN,
                                     setf.e1_err:np.NaN,
                                     setf.e2_err:np.NaN,
                                     setf.e1_err:1e99,
                                     setf.e2_err:1e99,})
                        
                    hdulist.append(table_to_hdu(shear_estimates_table))
                
            method_shear_estimates[method] = shear_estimates_table
            
            # Output the shear estimates
            hdulist.writeto(join(args.workdir,shear_estimates_filename),clobber=True)
        
    else:
    
        for filename in shear_estimates_prod.get_all_filenames():
            
            hdulist = fits.HDUList()
            
            for j in range(num_detectors):
                
                shm_hdu = table_to_hdu(initialise_shear_estimates_table(detector=j))
                hdulist.append(shm_hdu)
                
            hdulist.writeto(join(args.workdir,filename),clobber=True)
    
    write_pickled_product(shear_estimates_prod,
                          join(args.workdir,args.shear_estimates_product))
    
    logger.info("Finished shear estimation.")
    
    return
