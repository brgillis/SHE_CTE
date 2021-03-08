""" @file ObjectIdSplit.py

    Created 14 Mar 2019

    Split point executable for splitting up processing of objects into batches.
"""
import argparse
from copy import deepcopy
import math
import os

from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_xml_product, write_xml_product,
                             get_allowed_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import read_analysis_config, AnalysisConfigKeys
from SHE_PPT.products import she_object_id_list, mer_final_catalog
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import get_arguments_string
from astropy.table import Table

import SHE_CTE
import numpy as np


__updated__ = "2021-03-08"

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


default_batch_size = 20
default_max_batches = 0

logger = getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ObjectIdSplit defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Input arguments
    parser.add_argument('--mer_final_catalog_tables', type=str)

    parser.add_argument('--data_images', type=str, default=None,
                        help='.json listfile containing filenames of data image products.')

    parser.add_argument("--pipeline_config", default=None, type=str,
                        help="Pipeline-wide configuration file.")

    # Output arguments
    parser.add_argument('--object_ids', type=str)

    # Required pipeline arguments
    parser.add_argument('--workdir', type=str,)
    parser.add_argument('--logdir', type=str,)
    parser.add_argument('--debug', action='store_true',
                        help="Set to enable debugging protocols")
    parser.add_argument('--profile', action='store_true')

    logger.debug('# Exiting SHE_CTE_ObjectIdSplit defineSpecificProgramOptions()')

    return parser


def object_id_split_from_args(args):
    """ Core function for implementing a split by object ID
    """

    logger.debug('# Entering object_id_split_from_args(args)')

    # Read in the pipeline configuration if present
    if args.pipeline_config is None or args.pipeline_config == "None":
        logger.warning("No pipeline configuration found. Using default batch size of " + str(default_batch_size))
        batch_size = default_batch_size
        max_batches = default_max_batches
        ids_to_use = None
    else:
        pipeline_config = read_analysis_config(args.pipeline_config, workdir=args.workdir)

        # Check for the batch size key
        if AnalysisConfigKeys.OID_BATCH_SIZE.value not in pipeline_config:
            logger.info("Key " + AnalysisConfigKeys.OID_BATCH_SIZE.value + " not found in pipeline config " + args.pipeline_config + ". " +
                        "Using default batch size of " + str(default_batch_size))
            batch_size = default_batch_size
        else:
            batch_size = int(pipeline_config[AnalysisConfigKeys.OID_BATCH_SIZE.value])
            if batch_size < 0:
                raise ValueError("Invalid batch size: " + str(batch_size) + ". Must be >= 0.")
            logger.info("Using batch size of: " + str(batch_size))

        # Check for the max_batches key
        if AnalysisConfigKeys.OID_MAX_BATCHES.value not in pipeline_config:
            logger.info("Key " + AnalysisConfigKeys.OID_MAX_BATCHES.value + " not found in pipeline config " + args.pipeline_config + ". " +
                        "Using default max batches of " + str(default_max_batches))
            max_batches = default_max_batches
        else:
            max_batches = int(pipeline_config[AnalysisConfigKeys.OID_MAX_BATCHES.value])
            if max_batches < 0:
                raise ValueError("Invalid max batches: " + str(max_batches) + ". Must be >= 0.")
            logger.info("Using max batches of: " + str(max_batches))

        # Check for the IDs key
        if AnalysisConfigKeys.OID_IDS.value not in pipeline_config:
            logger.info("Key " + AnalysisConfigKeys.OID_IDS.value + " not found in pipeline config " + args.pipeline_config + ". " +
                        "Using default of using all IDs")
            ids_to_use = None
        else:
            ids_to_use_str = pipeline_config[AnalysisConfigKeys.OID_IDS.value].split()

            if len(ids_to_use_str) == 0:
                ids_to_use = None
                logger.info("Using all IDs")
            else:
                ids_to_use = list(map(int, ids_to_use_str))
                logger.info("Using limited selection of IDs:" + str(ids_to_use))

    # Read in the data images

    logger.info("Reading data images...")

    data_stack = SHEFrameStack.read(exposure_listfile_filename=args.data_images,
                                    detections_listfile_filename=args.mer_final_catalog_tables,
                                    workdir=args.workdir,
                                    save_products=True,
                                    memmap=True,
                                    mode='denywrite')

    # Get the tile list, pointing list, and observation id from the input data products

    tile_list = []

    for detections_catalogue_product in data_stack.detections_catalogue_products:
        tile_index = detections_catalogue_product.Data.TileIndex
        tile_product_id = detections_catalogue_product.Header.ProductId
        tile_list.append((tile_index, tile_product_id))

    pointing_list = []
    observation_id = None

    for exposure_product in data_stack.exposure_products:
        pointing_list.append(exposure_product.Data.ObservationSequence.PointingId)
        if observation_id is None:
            observation_id = exposure_product.Data.ObservationSequence.ObservationId
        else:
            # Check it's consistent and warn if not
            new_observation_id = exposure_product.Data.ObservationSequence.ObservationId
            if not observation_id == new_observation_id:
                logger.warning(f"Inconsistent observation IDs: {observation_id} and {new_observation_id}.")

    if ids_to_use is not None:

        all_ids_array = np.array(ids_to_use)

    else:

        all_ids = set(data_stack.detections_catalogue[mfc_tf.ID].data)

        logger.info("Finished reading in IDs from mer_final_catalog.")

        # Prune IDs that aren't in the images
        good_ids = []

        for gal_id in all_ids:

            # Get a stack of the galaxy images
            gal_stamp_stack = data_stack.extract_galaxy_stack(gal_id=gal_id,
                                                              width=1,)

            # Do we have any data for this object?
            if not gal_stamp_stack.is_empty():
                good_ids.append(gal_id)

        all_ids_array = np.array(good_ids)

    num_ids = len(all_ids_array)

    # If batch size is zero, use all IDs in one batch
    if batch_size == 0:
        batch_size = num_ids

    num_batches = int(math.ceil(num_ids / batch_size))

    if max_batches > 0:
        limited_num_batches = np.min((num_batches, max_batches))
    else:
        limited_num_batches = num_batches

    logger.info("Splitting IDs into " + str(limited_num_batches) + " batches of size " + str(batch_size) + ".")

    id_split_indices = np.linspace(batch_size, (num_batches - 1) * batch_size,
                                   num_batches - 1, endpoint=True, dtype=int)

    id_arrays = np.split(all_ids_array, id_split_indices)

    # Perform some quick sanity checks
    assert len(id_arrays) == num_batches
    assert len(id_arrays[0]) == batch_size or num_batches == 1
    assert len(id_arrays[-1]) <= batch_size

    # Start outputting the ID lists for each batch

    # Keep a list of all product filenames
    id_list_product_filename_list = []

    logger.info("Writing ID lists into products.")

    for i in range(limited_num_batches):

        logger.debug("Writing ID list #" + str(i) + " to product.")

        # For the filename, we want to set it up in a subfolder so we don't get too many files
        subfolder_number = i % 256
        subfolder_name = "data/s" + str(subfolder_number)

        qualified_subfolder_name = os.path.join(args.workdir, subfolder_name)

        if not os.path.exists(qualified_subfolder_name):
            # Can we create it?
            try:
                os.mkdir(qualified_subfolder_name)
            except Exception as e:
                logger.error("Directory (" + qualified_subfolder_name + ") does not exist and cannot be created.")
                raise e

        # Get a filename for this batch and store it in the list
        batch_id_list_product_filename = get_allowed_filename(type_name="OBJ-ID-LIST",
                                                              instance_id=str(i),
                                                              extension=".xml",
                                                              version=SHE_CTE.__version__,
                                                              subdir=subfolder_name,
                                                              processing_function="SHE")
        id_list_product_filename_list.append(batch_id_list_product_filename)

        # Create the product and fill out metadata
        batch_id_list_product = she_object_id_list.create_dpd_she_object_id_list(id_list=list(id_arrays[i]))

        batch_id_list_product.Data.PointingIdList = observation_id
        batch_id_list_product.Data.ObservationId = observation_id

        num_tiles = len(tile_list)

        base_tile_object = deepcopy(batch_id_list_product.Data.TileList[0])
        batch_id_list_product.Data.TileList = [base_tile_object] * num_tiles
        tile_list_binding = batch_id_list_product.Data.TileList

        for tile_i, (tile_index, tile_product_id) in enumerate(tile_list):
            tile_list_binding[tile_i].TileIndex = tile_index
            tile_list_binding[tile_i].TileProductId = tile_product_id

        # Save the product
        write_xml_product(batch_id_list_product, batch_id_list_product_filename, workdir=args.workdir)

        logger.debug("Successfully wrote ID list #" + str(i) + " to product: " + batch_id_list_product_filename)

    # Output the listfile
    write_listfile(os.path.join(args.workdir, args.object_ids), id_list_product_filename_list)
    logger.info("Finished writing listfile of object ID list products to " + args.object_ids)

    logger.debug('# Exiting object_id_split_from_args normally')

    return


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ObjectIdSplit mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_ObjectIdSplit",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    try:

        if args.profile:
            import cProfile
            cProfile.runctx("object_id_split_from_args(args)", {},
                            {"object_id_split_from_args": object_id_split_from_args,
                             "args": args, },
                            filename="object_id_split.prof")
        else:
            object_id_split_from_args(args)
    except Exception as e:
        # logger.warning("Failsafe exception block triggered with exception: " + str(e))
        raise

    logger.debug('# Exiting SHE_CTE_ObjectIdSplit mainMethod()')

    return


def main():
    """
    @brief
        Alternate entry point for non-Elements execution.
    """

    parser = defineSpecificProgramOptions()

    args = parser.parse_args()

    mainMethod(args)

    return


if __name__ == "__main__":
    main()
