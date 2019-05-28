""" @file ObjectIdSplit.py

    Created 14 Mar 2019

    Split point executable for splitting up processing of objects into batches.
"""

__updated__ = "2019-05-27"

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

import argparse
import math
import os

import SHE_CTE
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_xml_product, write_xml_product,
                             get_allowed_filename, find_file)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import read_config, ConfigKeys
from SHE_PPT.products import object_id_list, detections
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import get_arguments_string
from astropy.table import Table
import numpy as np

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
    parser.add_argument('--detections_tables', type=str)

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
        logger.warn("No pipeline configuration found. Using default batch size of " + str(default_batch_size))
        batch_size = default_batch_size
    else:
        pipeline_config = read_config(args.pipeline_config, workdir=args.workdir)

        # Check for the batch size key
        if ConfigKeys.OID_BATCH_SIZE.value not in pipeline_config:
            logger.info("Key " + ConfigKeys.OID_BATCH_SIZE.value + " not found in pipeline config " + args.pipeline_config + ". " +
                        "Using default batch size of " + str(default_batch_size))
            batch_size = default_batch_size
        else:
            batch_size = int(pipeline_config[ConfigKeys.OID_BATCH_SIZE.value])
            if batch_size < 0:
                raise ValueError("Invalid batch size: " + str(batch_size) + ". Must be >= 0.")
            logger.info("Using batch size of: " + str(batch_size))

        # Check for the max_batches key
        if ConfigKeys.OID_MAX_BATCHES.value not in pipeline_config:
            logger.info("Key " + ConfigKeys.OID_MAX_BATCHES.value + " not found in pipeline config " + args.pipeline_config + ". " +
                        "Using default max batches of " + str(default_max_batches))
            max_batches = default_max_batches
        else:
            max_batches = int(pipeline_config[ConfigKeys.OID_MAX_BATCHES.value])
            if max_batches < 0:
                raise ValueError("Invalid max batches: " + str(max_batches) + ". Must be >= 0.")
            logger.info("Using max batches of: " + str(max_batches))

    # Read in each detections table and add the IDs in it to a global set
    all_ids = set()

    logger.info("Reading in IDs from detections tables from: " + args.detections_tables)

    detections_table_product_filenames = read_listfile(os.path.join(args.workdir, args.detections_tables))

    for tile_detections_table_product_filename in detections_table_product_filenames:

        # Read in the product and get the filename of the table

        tile_detections_table_product = read_xml_product(os.path.join(args.workdir,
                                                                      tile_detections_table_product_filename))

        if not isinstance(tile_detections_table_product, detections.dpdMerFinalCatalog):
            raise TypeError("Detections product is of invalid type: " + type(tile_detections_table_product))

        tile_detections_table_filename = tile_detections_table_product.get_data_filename()

        # Read in the table

        tile_detections_table = Table.read(os.path.join(args.workdir, tile_detections_table_filename))

        if not is_in_format(tile_detections_table, detf, verbose=True, ignore_metadata=True):
            raise TypeError("Input detections table is of invalid format.")

        # Get the ID list from it and add it to the set
        all_ids.update(tile_detections_table[detf.ID].data)

    logger.info("Finished reading in IDs from detections table.")

    # Convert IDs to a numpy array for easier handling
    all_ids_array = np.array(list(all_ids))
    num_ids = len(all_ids_array)

    # If batch size is zero, use all IDs in one batch
    if batch_size == 0:
        batch_size = num_ids

    num_batches = int(math.ceil(num_ids / batch_size))

    if max_batches > 0:
        limited_num_batches = np.min((num_batches, max_batches))

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

        # Create and save the product
        batch_id_list_product = object_id_list.create_dpd_she_object_id_list(id_list=list(id_arrays[i]))
        write_xml_product(batch_id_list_product, os.path.join(args.workdir, batch_id_list_product_filename))

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
        # logger.warn("Failsafe exception block triggered with exception: " + str(e))
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
