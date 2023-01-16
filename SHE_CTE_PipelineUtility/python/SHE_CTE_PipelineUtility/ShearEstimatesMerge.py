""" @file ShearEstimatesMerge.py

    Created 14 Mar 2019

    Merge point executable for the split of shear estimation over object ID, merging into a single output product
    per Field of View.
"""

__updated__ = "2021-08-18"

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

import multiprocessing as mp
import os
from typing import Any, Dict, Tuple, Type, Union

import numpy as np
from astropy import table

import SHE_CTE
from SHE_CTE.executor import CteLogOptions, SheCteExecutor
from SHE_PPT import products
from SHE_PPT.argument_parser import SheArgumentParser
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.constants.shear_estimation_methods import ShearEstimationMethods
from SHE_PPT.executor import ReadConfigArgs
from SHE_PPT.file_io import (get_allowed_filename, read_listfile, read_xml_product, write_xml_product)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (AnalysisConfigKeys, ConfigKeys, read_analysis_config)
from SHE_PPT.table_formats.she_lensmc_chains import tf as lmcc_tf
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_utility import is_in_format

# Set up dicts for pipeline config defaults and types
D_SEM_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    AnalysisConfigKeys.SEM_NUM_THREADS: 8,
    }

D_SEM_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    AnalysisConfigKeys.SEM_NUM_THREADS: int,
    }

D_SEM_CONFIG_CLINE_ARGS: Dict[ConfigKeys, str] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    AnalysisConfigKeys.SEM_NUM_THREADS: "number_threads",
    }

EXEC_NAME = "SHE_CTE_ShearEstimatesMerge"

logger = getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger.debug('#')
    logger.debug(f'# Entering {EXEC_NAME} defineSpecificProgramOptions()')
    logger.debug('#')

    parser = SheArgumentParser()

    # Input data
    parser.add_input_arg('--shear_estimates_product_listfile', type = str)
    parser.add_input_arg('--she_lensmc_chains_listfile', type = str)

    # Output data
    parser.add_output_arg('--merged_she_measurements', type = str)
    parser.add_output_arg('--merged_she_lensmc_chains', type = str,
                          help = 'XML data product to contain LensMC chains data.')

    # Input arguments (can't be used in pipeline)
    parser.add_option_arg('--number_threads', type = int,
                          help = 'Number of parallel threads to use.')

    logger.debug(f'# Exiting {EXEC_NAME} defineSpecificProgramOptions()')

    return parser


def read_lensmc_chains_tables(she_lensmc_chains_table_product_filename, workdir):
    try:

        logger.debug("Loading chains from file: " + she_lensmc_chains_table_product_filename)

        # Read in the product and get the filename of the table

        she_lensmc_chains_product = read_xml_product(
            os.path.join(workdir, she_lensmc_chains_table_product_filename))

        if not isinstance(she_lensmc_chains_product, products.she_lensmc_chains.dpdSheLensMcChains):
            raise TypeError(f"Chains product is of invalid type: {type(she_lensmc_chains_product)}")

    except Exception as e:

        logger.warning("Failsafe block encountered exception: " + str(e))
        return

    try:

        she_lensmc_chains_table_filename = she_lensmc_chains_product.get_filename()

        if she_lensmc_chains_table_filename is None or she_lensmc_chains_table_filename == "None":

            logger.debug("No chains avaialble from file: " + she_lensmc_chains_table_product_filename)

        else:

            she_lensmc_chains_table = table.Table.read(os.path.join(
                workdir, she_lensmc_chains_table_filename))

            if not is_in_format(she_lensmc_chains_table, lmcc_tf, verbose = True, ignore_metadata = True):
                raise TypeError("Input chains table is of invalid format.")

    except Exception as e:

        logger.warning("Failsafe block encountered exception: " + str(e))
        return

    logger.debug("Finished loading chains from file: " + she_lensmc_chains_table_product_filename)

    return she_lensmc_chains_table


def read_method_estimates_tables(she_measurements_table_product_filename, workdir):
    observation_id = None
    observation_time = None
    pointing_id_list = None
    tile_list = None

    try:

        logger.debug("Loading shear estimates from file: " + she_measurements_table_product_filename)

        # Read in the product and get the filename of the table

        she_measurements_table_product = read_xml_product(
            os.path.join(workdir, she_measurements_table_product_filename))

        if not isinstance(she_measurements_table_product, products.she_measurements.dpdSheMeasurements):
            raise TypeError(f"Shear product is of invalid type: {type(she_measurements_table_product)}")

        if observation_id is None:
            observation_id = she_measurements_table_product.Data.ObservationId
            observation_time = she_measurements_table_product.Data.ObservationDateTime
            pointing_id_list = she_measurements_table_product.Data.PointingIdList
            tile_list = she_measurements_table_product.Data.TileList

    except Exception as e:

        logger.warning("Failsafe block encountered exception: " + str(e))
        return

    # Loop over methods and read in the table

    she_measurements_tables = {}

    for method in ShearEstimationMethods:

        try:

            she_measurements_method_table_filename = she_measurements_table_product.get_method_filename(method)

            if she_measurements_method_table_filename is None or she_measurements_method_table_filename == "None":
                logger.debug(f"No shear estimates available for method {method}")
                continue

            she_measurements_method_table = table.Table.read(os.path.join(
                workdir, she_measurements_method_table_filename))

            if not is_in_format(she_measurements_method_table, sm_tf, verbose = True, ignore_metadata = True):
                raise TypeError(f"Input shear estimates table for method {method} is of invalid format.")

            # Append the table to the list of tables
            she_measurements_tables[method] = she_measurements_method_table

        except Exception as e:

            logger.warning("Failsafe block encountered exception: " + str(e))
            return

        logger.debug("Finished loading shear estimates from file: " + she_measurements_table_product_filename)

    return (she_measurements_tables, observation_id, observation_time, pointing_id_list, tile_list)


def she_measurements_merge_from_args(args):
    """ Core function for implementing a merge of shear estimates tables
    """

    logger.debug('# Entering she_measurements_merge_from_args(args)')

    # Determine how many threads we'll use
    number_threads = args.pipeline_config[AnalysisConfigKeys.SEM_NUM_THREADS]

    # If number_threads is 0 or lower, assume it means this many fewer than the cpu count
    if number_threads <= 0:
        number_threads = max(1, mp.cpu_count() + number_threads)

    # Keep a list of all shear estimates tables for each method
    she_measurements_tables = dict.fromkeys(ShearEstimationMethods)
    for method in ShearEstimationMethods:
        # Start with an empty list of the tables
        she_measurements_tables[method] = []

    logger.info("Loading shear estimates from files listed in: " + args.shear_estimates_product_listfile)

    measurements_product_filenames = read_listfile(
        os.path.join(args.workdir, args.shear_estimates_product_listfile))

    logger.info("Loading chains from files listed in: " + args.she_lensmc_chains_listfile)

    chains_product_filenames = read_listfile(
        os.path.join(args.workdir, args.she_lensmc_chains_listfile))

    

    # Read the measurements tables

    input_tuples = [(
        she_measurements_table_product_filename, args.workdir) for
        she_measurements_table_product_filename in measurements_product_filenames]

    with mp.Pool(processes = number_threads) as pool:
        pool_she_measurements_tables_and_metadata = pool.starmap(read_method_estimates_tables, input_tuples)

    (full_l_she_measurements_tables,
        full_l_observation_ids,
        full_l_observation_times,
        full_l_pointing_id_lists,
        full_l_tile_lists) = zip(*pool_she_measurements_tables_and_metadata)

    # Read the chains tables

    input_tuples = [(
        she_lensmc_chains_table_product_filename, args.workdir) for 
        she_lensmc_chains_table_product_filename in chains_product_filenames
        ]

    with mp.Pool(processes = number_threads) as pool:
        she_lensmc_chains_tables = pool.starmap(read_lensmc_chains_tables, input_tuples)


    l_she_measurements_tables = [x for x in full_l_she_measurements_tables if x is not None]
    l_observation_ids = [x for x in full_l_observation_ids if x is not None]
    l_observation_times = [x for x in full_l_observation_times if x is not None]
    l_pointing_id_lists = [x for x in full_l_pointing_id_lists if x is not None]
    l_tile_lists = [x for x in full_l_tile_lists if x is not None]

    # Check metadata is consistent
    for l in (l_observation_ids, l_observation_times, l_pointing_id_lists, l_tile_lists):
        l_arr = np.array(l)
        if not (l_arr == l_arr[0]).all():
            logger.warning("Metadata is not consistent through all batches. Will use metadata from first batch.")

    observation_id = l_observation_ids[0]
    observation_time = l_observation_times[0]
    pointing_id_list = l_pointing_id_lists[0]
    tile_list = l_tile_lists[0]

    # Sort the tables into the expected format
    for method in ShearEstimationMethods:

        for i in range(len(l_she_measurements_tables)):
            if method not in l_she_measurements_tables[i]:
                continue
            t = l_she_measurements_tables[i][method]
            if t is None or len(t) == 0:
                continue
            she_measurements_tables[method].append(t)

    logger.info("Finished loading shear estimates from files listed in: " + args.shear_estimates_product_listfile)

    logger.info("Combining shear estimates tables.")

    # Combine the shear estimates tables for each method and output them

    combined_she_measurements_tables = dict.fromkeys(ShearEstimationMethods)

    # Create the output products
    combined_she_measurements_product = products.she_measurements.create_she_measurements_product(
        spatial_footprint = os.path.join(args.workdir, measurements_product_filenames[0]))
    combined_she_lensmc_chains_product = products.she_lensmc_chains.create_lensmc_chains_product(
        spatial_footprint = os.path.join(args.workdir, measurements_product_filenames[0]))

    # Set the metadata for the measurements product
    combined_she_measurements_product.Data.ObservationId = observation_id
    combined_she_measurements_product.Data.ObservationDateTime = observation_time
    combined_she_measurements_product.Data.PointingIdList = pointing_id_list
    combined_she_measurements_product.Data.TileList = tile_list

    # Set the metadata for the chains product
    combined_she_lensmc_chains_product.Data.ObservationId = observation_id
    combined_she_lensmc_chains_product.Data.ObservationDateTime = observation_time
    combined_she_lensmc_chains_product.Data.PointingIdList = pointing_id_list
    combined_she_lensmc_chains_product.Data.TileList = tile_list

    for method in ShearEstimationMethods:

        # Skip if no data for this method
        if len(she_measurements_tables[method]) == 0:
            combined_she_measurements_product.set_method_filename(method, None)
            continue

        # Combine the tables
        combined_she_measurements_tables[method] = table.vstack(she_measurements_tables[method],
                                                                metadata_conflicts = "silent")

        # Get a filename for the table
        combined_she_measurements_table_filename = get_allowed_filename(type_name = "SHEAR-EST-" + method.name,
                                                                        instance_id = 'MERGED',
                                                                        extension = ".fits",
                                                                        version = SHE_CTE.__version__,
                                                                        subdir = "data",
                                                                        processing_function = "SHE")
        combined_she_measurements_product.set_method_filename(method, combined_she_measurements_table_filename)

        # Output the combined table
        combined_she_measurements_tables[method].write(os.path.join(args.workdir,
                                                                    combined_she_measurements_table_filename),
                                                       format = "fits")

        logger.info("Combined shear estimates for method " + method.value +
                    " output to: " + combined_she_measurements_table_filename)

    # Combine the chains tables

    # Combine the tables
    combined_she_lensmc_chains_tables = table.vstack(she_lensmc_chains_tables,
                                                     metadata_conflicts = "silent")

    # Get a filename for the table
    combined_she_lensmc_chains_table_filename = get_allowed_filename(type_name = "SHEAR-CHAIN",
                                                                     instance_id = 'MERGED',
                                                                     extension = ".fits",
                                                                     version = SHE_CTE.__version__,
                                                                     subdir = "data",
                                                                     processing_function = "SHE")
    combined_she_lensmc_chains_product.set_filename(combined_she_lensmc_chains_table_filename)

    # Output the combined table
    combined_she_lensmc_chains_tables.write(os.path.join(args.workdir,
                                                         combined_she_lensmc_chains_table_filename),
                                            format = "fits")

    logger.info("Combined shear estimates for method " + method.value +
                " output to: " + combined_she_measurements_table_filename)

    # Save the products
    write_xml_product(combined_she_measurements_product, args.merged_she_measurements, args.workdir)
    logger.info("Combined shear estimates product output to: " + args.merged_she_measurements)

    write_xml_product(combined_she_lensmc_chains_product, args.merged_she_lensmc_chains, args.workdir)
    logger.info("Combined chains product output to: " + args.merged_she_lensmc_chains)

    logger.debug('# Exiting she_measurements_merge_from_args normally')


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    executor = SheCteExecutor(run_from_args_function = she_measurements_merge_from_args,
                              log_options = CteLogOptions(executable_name = EXEC_NAME),
                              config_args = ReadConfigArgs(d_config_defaults = D_SEM_CONFIG_DEFAULTS,
                                                           d_config_types = D_SEM_CONFIG_TYPES,
                                                           d_config_cline_args = D_SEM_CONFIG_CLINE_ARGS,
                                                           s_config_keys_types = {AnalysisConfigKeys},
                                                           ), )

    executor.run(args, logger = logger, pass_args_as_dict = False)


def main():
    """
    @brief
        Alternate entry point for non-Elements execution.
    """

    parser = defineSpecificProgramOptions()

    args = parser.parse_args()

    mainMethod(args)


if __name__ == "__main__":
    main()
