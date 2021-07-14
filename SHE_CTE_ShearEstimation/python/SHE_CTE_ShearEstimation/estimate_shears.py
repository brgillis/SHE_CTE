""" @file estimate_shears.py

    Created 27 Mar 2017

    Primary execution loop for measuring galaxy shapes from an image file.
"""

__updated__ = "2021-06-23"

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

from SHE_PPT import flags
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT import mdb
from SHE_PPT import products
from SHE_PPT.file_io import (write_xml_product, get_allowed_filename, get_data_filename,
                             read_listfile, find_file)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (AnalysisConfigKeys, CalibrationConfigKeys, get_conditional_product)
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf
from SHE_PPT.table_formats.she_ksb_measurements import initialise_ksb_measurements_table, tf as ksbm_tf
from SHE_PPT.table_formats.she_lensmc_chains import initialise_lensmc_chains_table, tf as lmcc_tf, len_chain
from SHE_PPT.table_formats.she_lensmc_measurements import initialise_lensmc_measurements_table, tf as lmcm_tf
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_formats.she_momentsml_measurements import initialise_momentsml_measurements_table, tf as mmlm_tf
from SHE_PPT.table_formats.she_regauss_measurements import initialise_regauss_measurements_table, tf as regm_tf
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import hash_any
from astropy.io import fits
from astropy.io.fits import table_to_hdu

import SHE_CTE
from SHE_CTE_ShearEstimation.control_training_data import load_control_training_data
from SHE_CTE_ShearEstimation.galsim_estimate_shear import KSB_estimate_shear, REGAUSS_estimate_shear
from SHE_LensMC import she_measure_shear
from SHE_LensMC.training_data import load_training_data
from SHE_MomentsML.estimate_shear import estimate_shear as ML_estimate_shear
import numpy as np


logger = getLogger(__name__)

loading_methods = {"KSB": load_control_training_data,
                   "REGAUSS": load_control_training_data,
                   "MomentsML": None,
                   "LensMC": load_training_data}

estimation_methods = {"KSB": KSB_estimate_shear,
                      "REGAUSS": REGAUSS_estimate_shear,
                      "MomentsML": ML_estimate_shear,
                      "LensMC": she_measure_shear}

initialisation_methods = {"KSB": initialise_ksb_measurements_table,
                          "REGAUSS": initialise_regauss_measurements_table,
                          "MomentsML": initialise_momentsml_measurements_table,
                          "LensMC": initialise_lensmc_measurements_table}

table_formats = {"KSB": ksbm_tf,
                 "REGAUSS": regm_tf,
                 "MomentsML": mmlm_tf,
                 "LensMC": lmcm_tf}

default_chains_method = "LensMC"
default_fast_mode = "normal"


def fill_measurements_table_meta(t,
                                 mer_final_catalog_products,
                                 vis_calibrated_frame_products,
                                 observation_id=None,
                                 observation_time=None,
                                 pointing_id_list=None,
                                 tile_id_list=None):
    """ A function to get the Observation ID, Observation Time, Field ID, and Tile IDs from data products
        and add them to a table's meta. If values are provided, they will be used. Otherwise values will
        be determined from the data products.
    """

    if tile_id_list is None:
        # Get a list of the tile IDs from the MER catalogs
        tile_id_list = np.empty_like(mer_final_catalog_products, dtype=int)
        for i, mer_final_catalog_product in enumerate(mer_final_catalog_products):
            tile_id_list[i] = mer_final_catalog_product.Data.TileIndex

    # Turn the Tile ID list into a spaced string
    tile_id_list_str = " ".join(map(str, tile_id_list))

    if observation_time is None:

        # Get the observation data from the exposure products

        observation_times = [None] * len(vis_calibrated_frame_products)

        for i, vis_calibrated_frame_product in enumerate(vis_calibrated_frame_products):
            observation_times[i] = vis_calibrated_frame_product.Data.ObservationDateTime

        # Get the most recent observation time
        observation_times.sort(key=lambda t: t.OBT)
        observation_time = observation_times[-1]

    observation_time_str = str(observation_time.OBT)

    if observation_id is None:
        # Get the observation ID
        observation_ids = np.empty_like(vis_calibrated_frame_products, dtype=int)
        for i, vis_calibrated_frame_product in enumerate(vis_calibrated_frame_products):
            observation_ids[i] = vis_calibrated_frame_product.Data.ObservationSequence.ObservationId

        if not (observation_ids == observation_ids[0]).all():
            logger.warning("Not all observation IDs are the same. First in list will be used.")
        observation_id = observation_ids[0]

    if pointing_id_list is None:
        pointing_id_list = [None] * len(vis_calibrated_frame_products)
        for i, vis_calibrated_frame_product in enumerate(vis_calibrated_frame_products):
            pointing_id_list[i] = vis_calibrated_frame_product.Data.ObservationSequence.PointingId
    pointing_id_list_str = " ".join(map(str, pointing_id_list))

    # Update the table's meta with the observation and tile data
    t.meta[sm_tf.m.observation_id] = observation_id
    t.meta[sm_tf.m.observation_time] = observation_time_str
    t.meta[sm_tf.m.pointing_id] = pointing_id_list_str
    t.meta[sm_tf.m.tile_id] = tile_id_list_str

    return observation_id, observation_time, pointing_id_list, tile_id_list


def estimate_shears_from_args(args, dry_run=False):
    """
    @brief
        Perform shear estimation, given arguments from the command-line.

    @param kwargs <dict>

    @return None
    """

    logger.debug("Entering estimate_shears_from_args")

    # Load in the files in turn

    # Data images - Read in as SHEStack object

    if dry_run:
        dry_label = " mock dry"
    else:
        dry_label = ""

    # Load in the MDB
    if args.mdb is None:
        logger.warning("No MDB file provided as input. Default values will be used where necessary.")
        mdb.init(find_file("WEB/SHE_PPT_8_7/sample_mdb-SC8.xml"))
    elif args.mdb[-5:] == ".json":
        mdb_files = read_listfile(os.path.join(args.workdir, args.mdb))
        mdb.init(mdb_files=mdb_files, path=args.workdir)
    elif args.mdb[-4:] == ".xml":
        mdb.init(mdb_files=args.mdb, path=args.workdir)
    else:
        logger.warning("Unrecognized format for MDB file: " + os.path.splitext(args.mdb)[-1] +
                       ". Expected '.xml' or '.json'. Will attempt to proceed with default values.")

    logger.info("Reading " + dry_label + "data images...")

    data_stack = SHEFrameStack.read(exposure_listfile_filename=args.data_images,
                                    seg_listfile_filename=args.segmentation_images,
                                    stacked_image_product_filename=args.stacked_image,
                                    stacked_seg_product_filename=args.stacked_segmentation_image,
                                    psf_listfile_filename=args.psf_images_and_tables,
                                    detections_listfile_filename=args.detections_tables,
                                    object_id_list_product_filename=args.object_ids,
                                    workdir=args.workdir,
                                    save_products=True,
                                    memmap=False,
                                    load_images=False,
                                    mode='denywrite')

    # Read in the catalog and exposure data products, which we'll need for updating metadata
    mer_final_catalog_products = data_stack.detections_catalogue_products
    vis_calibrated_frame_products = data_stack.exposure_products

    # Calibration parameters product
    calibration_parameters_prod = get_conditional_product(args.calibration_parameters_product)
    if calibration_parameters_prod is not None and not isinstance(calibration_parameters_prod,
                                                                  products.she_psf_calibration_parameters.dpdShePsfCalibrationParameters):
        raise ValueError("CalibrationParameters product from " + os.path.join(args.workdir, args.calibration_parameters_product)
                         + " is invalid type.")

    # Set up method data filenames

    training_data_filenames = {"KSB": args.ksb_training_data,
                               "REGAUSS": args.regauss_training_data,
                               "MomentsML": args.momentsml_training_data,
                               "LensMC": args.lensmc_training_data}

    # Set up output

    logger.info("Generating shear estimates product...")

    # For the filename, we want to set it up in a subfolder so we don't get too many files
    subfolder_number = os.getpid() % 256
    subfolder_name = "data/s" + str(subfolder_number)

    qualified_subfolder_name = os.path.join(args.workdir, subfolder_name)

    if not os.path.exists(qualified_subfolder_name):
        # Can we create it?
        try:
            os.mkdir(qualified_subfolder_name)
        except Exception as e:
            logger.error("Directory (" + qualified_subfolder_name + ") does not exist and cannot be created.")
            raise e

    # Determine the instance ID to use for the estimates file
    qualified_stacked_image_data_filename = os.path.join(
        args.workdir, get_data_filename(args.stacked_image, workdir=args.workdir))
    with fits.open(qualified_stacked_image_data_filename, mode='denywrite', memmap=False, load_images=False) as f:
        header = f[0].header
        if ppt_mv.model_hash_label in header:
            estimates_instance_id = header[ppt_mv.model_hash_label]
        elif ppt_mv.obs_id_label in header:
            estimates_instance_id = str(header[ppt_mv.obs_id_label])
        else:
            logger.warning("Cannot determine proper instance ID for filenames. Using hash of image header.")
            estimates_instance_id = hash_any(header, format="base64")

        # Fix banned characters in the instance_id, add the pid, and enforce the maximum length
        estimates_instance_id = estimates_instance_id.replace('.', '-').replace('+', '-')
        estimates_instance_id = str(os.getpid()) + "-" + estimates_instance_id
        estimates_instance_id = estimates_instance_id[0:ppt_mv.short_instance_id_maxlen - 4]

    shear_estimates_prod = products.she_measurements.create_dpd_she_measurements(
        BFD_filename=None,
        KSB_filename=get_allowed_filename("KSB-SHM", estimates_instance_id,
                                          version=SHE_CTE.__version__,
                                          subdir=subfolder_name),
        LensMC_filename=get_allowed_filename("LensMC-SHM", estimates_instance_id,
                                             version=SHE_CTE.__version__,
                                             subdir=subfolder_name),
        MomentsML_filename=get_allowed_filename("MomentsML-SHM", estimates_instance_id,
                                                version=SHE_CTE.__version__,
                                                subdir=subfolder_name),
        REGAUSS_filename=get_allowed_filename("REGAUSS-SHM", estimates_instance_id,
                                              version=SHE_CTE.__version__,
                                              subdir=subfolder_name),
        spatial_footprint=os.path.join(args.workdir, args.stacked_image))

    # Init some values that we'll determine the first time possible
    observation_id = None
    observation_time = None
    pointing_id_list = None
    tile_id_list = None

    if not dry_run:

        method_shear_estimates = {}

        # Check for methods in the pipeline options

        # extract pipeline_config from the args, this has already been read in by
        # the calling program so no need to read it here
        pipeline_config = args.pipeline_config

        # Use methods specified in the command-line first
        if args.methods is not None and len(args.methods) > 0:
            methods = args.methods
        elif AnalysisConfigKeys.ES_METHODS.value in pipeline_config:
            methods = pipeline_config[AnalysisConfigKeys.ES_METHODS.value].split()
        elif CalibrationConfigKeys.ES_METHODS.value in pipeline_config:
            methods = pipeline_config[CalibrationConfigKeys.ES_METHODS.value].split()
        else:
            # Default to using all methods
            methods = list(estimation_methods.keys())

        # Determine which method we want chains from
        if args.chains_method is not None:
            chains_method = args.chains_method
        elif AnalysisConfigKeys.ES_CHAINS_METHOD.value in pipeline_config:
            chains_method = pipeline_config[AnalysisConfigKeys.ES_CHAINS_METHOD.value]
        elif CalibrationConfigKeys.ES_CHAINS_METHOD.value in pipeline_config:
            chains_method = pipeline_config[CalibrationConfigKeys.ES_CHAINS_METHOD.value]
        else:
            # Default to using all methods
            chains_method = default_chains_method

        if not chains_method in methods:
            raise ValueError("chains_method (\"" + str(chains_method) + "\") not in methods to run (" +
                             str(methods) + ").")

        # Determine whether or not to use fast mode
        if AnalysisConfigKeys.ES_FAST_MODE.value in pipeline_config:
            fast_mode_bool = pipeline_config[AnalysisConfigKeys.ES_FAST_MODE.value].lower() in ["true", "t"]
            if fast_mode_bool:
                fast_mode = "fast"
            else:
                fast_mode = "normal"
        else:
            fast_mode = default_fast_mode

        for method in methods:

            logger.info("Estimating shear with method " + method + "...")

            load_training_data = loading_methods[method]

            estimate_shear = estimation_methods[method]

            shear_estimates_filename = shear_estimates_prod.get_method_filename(method)

            if method == chains_method:
                return_chains = True
            else:
                return_chains = False

            hdulist = fits.HDUList()

            try:

                # Check we've supplied a method for it
                if estimate_shear is None:
                    raise NotImplementedError("No shear measurement method supplied for method " + method + ".")

                if load_training_data is not None:
                    training_data_filename = training_data_filenames[method]
                    if training_data_filename == 'None':
                        training_data_filename = None
                    if training_data_filename is None and method not in ("KSB", "LensMC", "REGAUSS"):
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

                # need to add ids to iterate over

                shear_estimates_results = estimate_shear(data_stack,
                                                         training_data=training_data,
                                                         calibration_data=calibration_data,
                                                         workdir=args.workdir,
                                                         debug=args.debug,
                                                         return_chains=return_chains,
                                                         mode=fast_mode)

                if return_chains or method == "LensMC":
                    shear_estimates_table, chains_table = shear_estimates_results
                else:
                    shear_estimates_table = shear_estimates_results

                # Update the table meta with observation and tile info
                (observation_id,
                 observation_time,
                 pointing_id_list,
                 tile_id_list) = fill_measurements_table_meta(t=shear_estimates_table,
                                                              mer_final_catalog_products=mer_final_catalog_products,
                                                              vis_calibrated_frame_products=vis_calibrated_frame_products,
                                                              observation_id=observation_id,
                                                              observation_time=observation_time,
                                                              pointing_id_list=pointing_id_list,
                                                              tile_id_list=tile_id_list)

                if return_chains:

                    chains_data_filename = get_allowed_filename(method.upper() + "-CHAINS", estimates_instance_id,
                                                                version=SHE_CTE.__version__,
                                                                subdir=subfolder_name)

                    chains_table.write(os.path.join(args.workdir, chains_data_filename))

                    chains_prod = products.she_lensmc_chains.create_lensmc_chains_product(chains_data_filename)

                    # Fill in metadata for the chains product
                    chains_prod.Data.ObservationId = observation_id
                    chains_prod.Data.ObservationDateTime = observation_time
                    chains_prod.Data.PointingIdList = pointing_id_list
                    if data_stack.object_id_list_product is not None:
                        chains_prod.Data.TileList = data_stack.object_id_list_product.Data.TileList

                    write_xml_product(chains_prod, os.path.join(args.workdir, args.she_lensmc_chains))

                if not is_in_format(shear_estimates_table, sm_tf):
                    raise ValueError("Invalid implementation: Shear estimation table returned in invalid format " +
                                     "for method " + method + ".")

                hdulist.append(table_to_hdu(shear_estimates_table))

                # Cleanup loaded data for this method
                del training_data, calibration_data

            except Exception as e:

                if isinstance(e, NotImplementedError) or (isinstance(e, ValueError) and
                                                          "Invalid implementation:" in str(e)):
                    logger.warning(str(e))
                else:
                    import traceback
                    logger.warning("Failsafe exception block triggered with exception: " + str(e) + ".\n"
                                   "Traceback: " + "".join(traceback.format_tb(e.__traceback__)))

                hdulist = fits.HDUList()

                # Create an empty estimates table
                shear_estimates_table = initialisation_methods[method]()
                tf = table_formats[method]

                for r in range(len(data_stack.detections_catalogue[mfc_tf.ID])):

                    # Fill it with NaN measurements and 1e99 errors

                    shear_estimates_table.add_row({tf.ID: data_stack.detections_catalogue[mfc_tf.ID][r],
                                                   tf.g1: np.NaN,
                                                   tf.g2: np.NaN,
                                                   tf.g1_err: 1e99,
                                                   tf.g2_err: 1e99,
                                                   tf.g1g2_covar: np.NaN,
                                                   tf.fit_class: 2,
                                                   tf.fit_flags: flags.flag_unclassified_failure,
                                                   tf.ra: data_stack.detections_catalogue[mfc_tf.gal_x_world][r],
                                                   tf.dec: data_stack.detections_catalogue[mfc_tf.gal_y_world][r], })

                hdulist.append(table_to_hdu(shear_estimates_table))

                if return_chains:
                    chains_data_filename = get_allowed_filename(method.upper() + "-CHAINS", estimates_instance_id,
                                                                version=SHE_CTE.__version__,
                                                                subdir=subfolder_name)

                    # Create an empty chains table
                    chains_table = initialise_lensmc_chains_table(optional_columns=[lmcc_tf.ra, lmcc_tf.dec])

                    for r in range(len(data_stack.detections_catalogue[mfc_tf.ID])):

                        # Fill it with NaN measurements and 1e99 errors

                        chains_table.add_row({lmcc_tf.ID: data_stack.detections_catalogue[mfc_tf.ID][r],
                                              lmcc_tf.g1: [np.NaN] * len_chain,
                                              lmcc_tf.g2: [np.NaN] * len_chain,
                                              lmcc_tf.fit_flags: flags.flag_unclassified_failure,
                                              lmcc_tf.ra: [data_stack.detections_catalogue[mfc_tf.gal_x_world][r]] * len_chain,
                                              lmcc_tf.dec: [data_stack.detections_catalogue[mfc_tf.gal_y_world][r]] * len_chain, })

                    chains_table.write(os.path.join(args.workdir, chains_data_filename))

                    chains_prod = products.she_lensmc_chains.create_lensmc_chains_product(chains_data_filename)

                    # Fill in metadata for the chains product if possible
                    if observation_id is not None:
                        chains_prod.Data.ObservationId = observation_id
                        chains_prod.Data.ObservationDateTime = observation_time
                        chains_prod.Data.PointingIdList = pointing_id_list
                        if data_stack.object_id_list_product is not None:
                            chains_prod.Data.TileList = data_stack.object_id_list_product.Data.TileList

                    write_xml_product(chains_prod, os.path.join(args.workdir, args.she_lensmc_chains))

            method_shear_estimates[method] = shear_estimates_table

            # Output the shear estimates
            hdulist.writeto(os.path.join(args.workdir, shear_estimates_filename), overwrite=True)

    else:  # Dry run

        for method in methods:

            filename = shear_estimates_prod.get_method_filename(method)

            hdulist = fits.HDUList()

            shear_estimates_table = initialisation_methods[method]()

            shm_hdu = table_to_hdu(shear_estimates_table)
            hdulist.append(shm_hdu)

            hdulist.writeto(os.path.join(args.workdir, filename), overwrite=True)

    # If we're not using all methods, don't write unused ones in the product
    for method in estimation_methods:
        if method not in methods:
            shear_estimates_prod.set_method_filename(method, None)

    # Fill in metadata for the estimates product
    shear_estimates_prod.Data.ObservationId = observation_id
    shear_estimates_prod.Data.ObservationDateTime = observation_time
    shear_estimates_prod.Data.PointingIdList = pointing_id_list
    if data_stack.object_id_list_product is not None:
        shear_estimates_prod.Data.TileList = data_stack.object_id_list_product.Data.TileList

    write_xml_product(shear_estimates_prod, args.shear_estimates_product, workdir=args.workdir)

    logger.info("Finished shear estimation.")

    logger.debug("# Exiting estimate_shears_from_args() successfully.")

    return
