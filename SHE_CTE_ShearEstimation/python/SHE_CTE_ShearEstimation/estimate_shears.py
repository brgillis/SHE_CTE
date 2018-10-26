""" @file estimate_shears.py

    Created 27 Mar 2017

    Primary execution loop for measuring galaxy shapes from an image file.
"""

__updated__ = "2018-10-15"


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
import pdb
import os

from astropy.io import fits

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_CTE_ShearEstimation.bfd_measure_moments import bfd_measure_moments, bfd_perform_integration
from SHE_CTE_ShearEstimation.control_training_data import load_control_training_data
from SHE_CTE_ShearEstimation.galsim_estimate_shear import KSB_estimate_shear, REGAUSS_estimate_shear
from SHE_LensMC.SHE_measure_shear import fit_frame_stack
from SHE_MomentsML.estimate_shear import estimate_shear as ML_estimate_shear
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT import products
from SHE_PPT.file_io import (read_xml_product, write_xml_product, get_allowed_filename, get_data_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import read_config
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.bfd_moments import initialise_bfd_moments_table, tf as setf_bfd
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.table_formats.shear_estimates import initialise_shear_estimates_table, tf as setf
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from SHE_PPT.utility import hash_any

import numpy as np


loading_methods = {"KSB": load_control_training_data,
                   "REGAUSS": load_control_training_data,
                   "MomentsML": None,
                   "LensMC": None,
                   "BFD": None}

estimation_methods = {"KSB": KSB_estimate_shear,
                      "REGAUSS": REGAUSS_estimate_shear,
                      "MomentsML": ML_estimate_shear,
                      "LensMC": fit_frame_stack,
                      "BFD": bfd_measure_moments}

methods_key = "SHE_CTE_EstimateShear_methods"


def estimate_shears_from_args(args, dry_run=False):
    """
    @brief
        Perform shear estimation, given arguments from the command-line.

    @param kwargs <dict>

    @return None
    """

    logger = getLogger(__name__)

    logger.debug("Entering estimate_shears_from_args")

    # Load in the files in turn

    # Data images - Read in as SHEStack object

    if dry_run:
        dry_label = " mock dry"
    else:
        dry_label = ""

    logger.info("Reading " + dry_label + "data images...")

    data_stack = SHEFrameStack.read(exposure_listfile_filename=args.data_images,
                                    seg_listfile_filename=args.segmentation_images,
                                    stacked_image_product_filename=args.stacked_image,
                                    stacked_seg_filename=args.stacked_segmentation_image,
                                    psf_listfile_filename=args.psf_images_and_tables,
                                    detections_listfile_filename=args.detections_tables,
                                    workdir=args.workdir,
                                    clean_detections=True,
                                    apply_sc3_fix=False)

    # Calibration parameters product

    if args.calibration_parameters_product is not None:

        logger.info("Reading " + dry_label + "calibration parameters...")

        calibration_parameters_prod = read_xml_product(os.path.join(args.workdir, args.calibration_parameters_product))
        if not isinstance(calibration_parameters_prod, products.calibration_parameters.DpdSheCalibrationParametersProduct):
            raise ValueError("CalibrationParameters product from " + os.path.join(args.workdir, args.calibration_parameters_product)
                             + " is invalid type.")

    else:
        calibration_parameters_prod = None

    # Set up method data filenames

    training_data_filenames = {"KSB": args.ksb_training_data,
                               "REGAUSS": args.regauss_training_data,
                               "MomentsML": args.momentsml_training_data,
                               "LensMC": args.lensmc_training_data,
                               "BFD": args.bfd_training_data}

    # Set up output

    logger.info("Generating shear estimates product...")

    # Determine the instance ID to use for the estimates file
    qualified_stacked_image_data_filename = os.path.join(
        args.workdir, get_data_filename(args.stacked_image, workdir=args.workdir))
    with fits.open(qualified_stacked_image_data_filename, mode='denywrite', memmap=True) as f:
        header = f[0].header
        if ppt_mv.model_hash_label in header:
            estimates_instance_id = header[ppt_mv.model_hash_label][0:ppt_mv.short_instance_id_maxlen]
        elif ppt_mv.field_id_label in header:
            estimates_instance_id = header[ppt_mv.field_id_label][0:ppt_mv.short_instance_id_maxlen]
        else:
            logger.warn("Cannot determine proper instance ID for filenames. Using hash of stacked image header.")
            estimates_instance_id = hash_any(header, format="base64", max_length=ppt_mv.short_instance_id_maxlen)

    shear_estimates_prod = products.shear_estimates.create_shear_estimates_product(
        BFD_filename=get_allowed_filename("BFD-SHM", estimates_instance_id),
        KSB_filename=get_allowed_filename("KSB-SHM", estimates_instance_id),
        LensMC_filename=get_allowed_filename("LensMC-SHM", estimates_instance_id),
        MomentsML_filename=get_allowed_filename("MomentsML-SHM", estimates_instance_id),
        REGAUSS_filename=get_allowed_filename("REGAUSS-SHM", estimates_instance_id))

    if not dry_run:

        method_shear_estimates = {}

        # Check for methods in the pipeline options
        pipeline_config = read_config(args.pipeline_config, workdir=args.workdir)

        if pipeline_config is None:
            pipeline_config = {}

        if methods_key in pipeline_config:
            methods = pipeline_config[methods_key].split()
        else:
            methods = args.methods

        # If no methods are specified, use all
        if len(methods) == 0:
            methods = list(estimation_methods.keys())

        for method in methods:

            logger.info("Estimating shear with method " + method + "...")

            load_training_data = loading_methods[method]

            estimate_shear = estimation_methods[method]

            shear_estimates_filename = shear_estimates_prod.get_method_filename(method)

            hdulist = fits.HDUList()

            try:

                # Check we've supplied a method for it
                if estimate_shear is None:
                    raise NotImplementedError("No shear measurement method supplied for method " + method + ".")

                if load_training_data is not None:
                    training_data_filename = training_data_filenames[method]
                    if training_data_filename == 'None':
                        training_data_filename = None
                    if training_data_filename is None:
                        # Don't raise for KSB, REGAUSS, and LensMC which allow default behaviour here
                        if method not in ("KSB", "LensMC", "REGAUSS"):
                            raise ValueError(
                                "Invalid implementation: No training data supplied for method " + method + ".")
                    training_data = load_training_data(training_data_filename, workdir=args.workdir)

                else:
                    training_data = None

                calibration_data = None

                if calibration_parameters_prod is not None:

                    method_calibration_filename = calibration_parameters_prod.get_method_filename(method)

                    if method_calibration_filename is not None:

                        # For now just leave as handle to fits file
                        calibration_data = fits.open(os.path.join(args.workdir, method_calibration_filename))

                shear_estimates_table = estimate_shear(data_stack,
                                                       training_data=training_data,
                                                       calibration_data=calibration_data,
                                                       workdir=args.workdir,
                                                       debug=args.debug)

                if not (is_in_format(shear_estimates_table, setf) or is_in_format(shear_estimates_table, setf_bfd)):
                    raise ValueError("Invalid implementation: Shear estimation table returned in invalid format " +
                                     "for method " + method + ".")

                hdulist.append(table_to_hdu(shear_estimates_table))

                # Cleanup loaded data for this method
                del training_data, calibration_data

            except Exception as e:

                if isinstance(e, NotImplementedError) or (isinstance(e, ValueError) and
                                                          "Invalid implementation:" in str(e)):
                    logger.warn(str(e))
                else:
                    logger.warn("Failsafe exception block triggered with exception: " + str(e))

                hdulist = fits.HDUList()

                # Create an empty estimates table
                shear_estimates_table = initialise_shear_estimates_table()

                for r in range(len(data_stack.detections_catalogue[detf.ID])):

                    # Fill it with NaN measurements and 1e99 errors

                    shear_estimates_table.add_row({setf.ID: data_stack.detections_catalogue[detf.ID][r],
                                                   setf.g1: np.NaN,
                                                   setf.g2: np.NaN,
                                                   setf.g1_err: 1e99,
                                                   setf.g2_err: 1e99,
                                                   setf.g1g2_covar: np.NaN,
                                                   setf.re: np.NaN,
                                                   setf.snr: np.NaN,
                                                   setf.x_world: data_stack.detections_catalogue[detf.gal_x_world][r],
                                                   setf.y_world: data_stack.detections_catalogue[detf.gal_y_world][r], })

                hdulist.append(table_to_hdu(shear_estimates_table))

            method_shear_estimates[method] = shear_estimates_table

            # Output the shear estimates
            hdulist.writeto(os.path.join(args.workdir, shear_estimates_filename), clobber=True)

            if method == 'BFD':
                bfd_perform_integration(os.path.join(args.workdir, shear_estimates_filename))

    else:  # Dry run

        for filename in shear_estimates_prod.get_all_filenames():

            hdulist = fits.HDUList()

            shm_hdu = table_to_hdu(initialise_shear_estimates_table())
            hdulist.append(shm_hdu)

            hdulist.writeto(os.path.join(args.workdir, filename), clobber=True)

    # If we're not using all methods, don't write unused ones in the product
    for method in estimation_methods:
        if method not in methods:
            shear_estimates_prod.set_method_filename(method, None)

    write_xml_product(shear_estimates_prod,
                      os.path.join(args.workdir, args.shear_estimates_product))

    logger.info("Finished shear estimation.")

    return
