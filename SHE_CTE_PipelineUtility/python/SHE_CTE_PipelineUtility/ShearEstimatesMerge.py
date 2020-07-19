""" @file ShearEstimatesMerge.py

    Created 14 Mar 2019

    Merge point executable for the split of shear estimation over object ID, merging into a single output product
    per Field of View.
"""

__updated__ = "2020-07-10"

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
import os

from SHE_PPT import products
from SHE_PPT.file_io import (read_listfile,
                             read_xml_product, write_xml_product,
                             get_allowed_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import read_config, ConfigKeys
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import get_arguments_string
from astropy import table

import SHE_CTE
import multiprocessing as mp

logger = getLogger(__name__)

methods = ("KSB", "REGAUSS", "LensMC", "BFD", "MomentsML")

default_number_threads = 8


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ShearEstimatesMerge defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Input data
    parser.add_argument('--input_she_measurements_listfile', type=str)

    parser.add_argument("--pipeline_config", default=None, type=str,
                        help="Pipeline-wide configuration file.")

    # Output data
    parser.add_argument('--output_she_measurements', type=str)

    # Input arguments
    parser.add_argument('--number_threads', type=int, default=None,
                        help='Number of parallel threads to use.')

    # Required pipeline arguments
    parser.add_argument('--workdir', type=str,)
    parser.add_argument('--logdir', type=str,)
    parser.add_argument('--debug', action='store_true',
                        help="Set to enable debugging protocols")
    parser.add_argument('--profile', action='store_true')

    logger.debug('# Exiting SHE_CTE_ShearEstimatesMerge defineSpecificProgramOptions()')

    return parser


def read_method_estimates_tables(she_measurements_table_product_filename, workdir):

    try:

        logger.debug("Loading shear estimates from file: " + she_measurements_table_product_filename)

        # Read in the product and get the filename of the table

        she_measurements_table_product = read_xml_product(
            os.path.join(workdir, she_measurements_table_product_filename))

        if not isinstance(she_measurements_table_product, products.she_measurements.dpdShearMeasurement):
            raise TypeError("Shear product is of invalid type: " + type(she_measurements_table_product))

    except Exception as e:

        logger.warning("Failsafe block encountered exception: " + str(e))
        return

    # Loop over methods and read in the table

    she_measurements_tables = {}

    for method in methods:

        try:

            she_measurements_method_table_filename = she_measurements_table_product.get_method_filename(method)

            if she_measurements_method_table_filename is None or she_measurements_method_table_filename == "None":
                logger.debug("No shear estimates available for method: " + method)
                continue

            she_measurements_method_table = table.Table.read(os.path.join(
                workdir, she_measurements_method_table_filename))

            if not is_in_format(she_measurements_method_table, sm_tf, verbose=True, ignore_metadata=True):
                raise TypeError("Input shear estimates table for method {} is of invalid format.".format(method))

            # Append the table to the list of tables
            she_measurements_tables[method] = she_measurements_method_table

        except Exception as e:

            logger.warning("Failsafe block encountered exception: " + str(e))
            return

        logger.debug("Finished loading shear estimates from file: " + she_measurements_table_product_filename)

    return she_measurements_tables


def she_measurements_merge_from_args(args):
    """ Core function for implementing a merge of shear estimates tables
    """

    logger.debug('# Entering she_measurements_merge_from_args(args)')

    # Read in the pipeline config
    try:
        pipeline_config = read_config(args.pipeline_config, workdir=args.workdir)
        if pipeline_config is None:
            pipeline_config = {}
    except Exception as e:
        logger.warning("Failsafe exception block triggered when trying to read pipeline config. " +
                    "Exception was: " + str(e))
        pipeline_config = {}

    # Determine how many threads we'll use
    if args.number_threads is not None:
        number_threads = args.number_threads
    elif ConfigKeys.SEM_NUM_THREADS.value in pipeline_config:
        number_threads = pipeline_config[ConfigKeys.MB_NUM_THREADS.value]
        if number_threads.lower() == "none":
            number_threads = default_number_threads
    else:
        number_threads = default_number_threads

    # If number_threads is 0 or lower, assume it means this many fewer than the cpu count
    if number_threads <= 0:
        number_threads = max(1, mp.cpu_count() + number_threads)

    # Keep a list of all shear estimates tables for each method
    she_measurements_tables = dict.fromkeys(methods)
    for method in methods:
        # Start with an empty list of the tables
        she_measurements_tables[method] = []

    logger.info("Loading shear estimates from files listed in: " + args.input_she_measurements_listfile)

    she_measurements_table_product_filenames = read_listfile(
        os.path.join(args.workdir, args.input_she_measurements_listfile))

    # If using just one thread, don't bother with multiprocessing to read tables

    if number_threads == 1:

        full_l_she_measurements_tables = [read_method_estimates_tables(
            f, args.workdir) for f in she_measurements_table_product_filenames]

        l_she_measurements_tables = [t for t in full_l_she_measurements_tables if t is not None]

    else:

        pool = mp.Pool(processes=number_threads)
        pool_she_measurements_tables = [pool.apply_async(read_method_estimates_tables, args=(
            she_measurements_table_product_filename, args.workdir)) for she_measurements_table_product_filename in she_measurements_table_product_filenames]

        pool.close()
        pool.join()

        l_she_measurements_tables = [a.get() for a in pool_she_measurements_tables if a.get() is not None]

    # Sort the tables into the expected format
    for method in methods:

        for i in range(len(l_she_measurements_tables)):
            if method not in l_she_measurements_tables[i]:
                continue
            t = l_she_measurements_tables[i][method]
            if t is None or len(t) == 0:
                continue
            she_measurements_tables[method].append(t)

    logger.info("Finished loading shear estimates from files listed in: " + args.input_she_measurements_listfile)

    logger.info("Combining shear estimates tables.")

    # Combine the shear estimates tables for each method and output them

    combined_she_measurements_tables = dict.fromkeys(methods)

    # Create the output product
    combined_she_measurements_product = products.she_measurements.create_she_measurements_product(
        spatial_footprint=os.path.join(args.workdir, she_measurements_table_product_filenames[0]))

    for method in methods:

        # Skip if no data for this method
        if len(she_measurements_tables[method]) == 0:
            combined_she_measurements_product.set_method_filename(method, "None")
            continue

        # Combine the tables
        combined_she_measurements_tables[method] = table.vstack(she_measurements_tables[method])

        # Get a filename for the table
        combined_she_measurements_table_filename = get_allowed_filename(type_name="SHEAR-EST-" + method.upper(),
                                                                       instance_id='MERGED',
                                                                       extension=".fits",
                                                                       version=SHE_CTE.__version__,
                                                                       subdir="data",
                                                                       processing_function="SHE")
        combined_she_measurements_product.set_method_filename(method, combined_she_measurements_table_filename)

        # Output the table
        combined_she_measurements_tables[method].write(os.path.join(args.workdir,
                                                                   combined_she_measurements_table_filename),
                                                      format="fits")

        logger.info("Combined shear estimates for method " + method +
                    " output to: " + combined_she_measurements_table_filename)

    # Save the product
    write_xml_product(combined_she_measurements_product, args.output_she_measurements, args.workdir)

    logger.info("Combined shear estimates product output to: " + args.output_she_measurements)

    logger.debug('# Exiting she_measurements_merge_from_args normally')

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
    logger.debug('# Entering SHE_CTE_ShearEstimatesMerge mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_ShearEstimatesMerge",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    try:

        if args.profile:
            import cProfile
            cProfile.runctx("she_measurements_merge_from_args(args)", {},
                            {"she_measurements_merge_from_args": she_measurements_merge_from_args,
                             "args": args, },
                            filename="she_measurements_merge.prof")
        else:
            she_measurements_merge_from_args(args)
    except Exception as e:
        # logger.warning("Failsafe exception block triggered with exception: " + str(e))
        raise

    logger.debug('# Exiting SHE_CTE_ShearEstimatesMerge mainMethod()')

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
