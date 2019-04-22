""" @file ShearEstimatesMerge.py

    Created 14 Mar 2019

    Merge point executable for the split of shear estimation over object ID, merging into a single output product
    per Field of View.
"""

__updated__ = "2019-04-22"

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

import SHE_CTE
from SHE_PPT import products
from SHE_PPT.file_io import (read_listfile,
                             read_xml_product, write_xml_product,
                             get_allowed_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.shear_estimates import tf as setf
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import get_arguments_string
from astropy import table


logger = getLogger(__name__)

methods = ("KSB", "REGAUSS", "LensMC", "BFD", "MomentsML")


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

    # Input arguments
    parser.add_argument('--input_shear_estimates_listfile', type=str)

    # Output arguments
    parser.add_argument('--output_shear_estimates', type=str)

    # Required pipeline arguments
    parser.add_argument('--workdir', type=str,)
    parser.add_argument('--logdir', type=str,)
    parser.add_argument('--debug', action='store_true',
                        help="Set to enable debugging protocols")
    parser.add_argument('--profile', action='store_true')

    logger.debug('# Exiting SHE_CTE_ShearEstimatesMerge defineSpecificProgramOptions()')

    return parser


def shear_estimates_merge_from_args(args):
    """ Core function for implementing a merge of shear estimates tables
    """

    logger.debug('# Entering shear_estimates_merge_from_args(args)')

    # Keep a list of all shear estimates tables for each method
    shear_estimates_tables = dict.fromkeys(methods)
    for method in methods:
        # Start with an empty list of the tables
        shear_estimates_tables[method] = []

    logger.info("Loading shear estimates from files listed in: " + args.input_shear_estimates_listfile)

    shear_estimates_table_product_filenames = read_listfile(
        os.path.join(args.workdir, args.input_shear_estimates_listfile))

    # We'll use a failsafe exception block in case of corrupt data, and combine what we can. But raise
    # an exception if we can't get anything
    at_least_one_good_estimates_table = False

    for shear_estimates_table_product_filename in shear_estimates_table_product_filenames:

        try:

            logger.debug("Loading shear estimates from file: " + shear_estimates_table_product_filename)

            # Read in the product and get the filename of the table

            shear_estimates_table_product = read_xml_product(
                os.path.join(args.workdir, shear_estimates_table_product_filename))

            if not isinstance(shear_estimates_table_product, products.shear_estimates.dpdShearMeasurement):
                raise TypeError("Shear product is of invalid type: " + type(shear_estimates_table_product))

        except Exception as e:

            logger.warn("Failsafe block encountered exception: " + str(e))
            continue

        # Loop over methods and read in the table

        for method in methods:

            try:

                shear_estimates_method_table_filename = shear_estimates_table_product.get_method_filename(method)

                if shear_estimates_method_table_filename is None or shear_estimates_method_table_filename == "None":
                    logger.debug("No shear estimates available for method: " + method)
                    continue

                shear_estimates_method_table = table.Table.read(os.path.join(
                    args.workdir, shear_estimates_method_table_filename))

                if not is_in_format(shear_estimates_method_table, setf, verbose=True, ignore_metadata=True):
                    raise TypeError("Input shear estimates table for method {} is of invalid format.".format(method))

                # Append the table to the list of tables
                shear_estimates_tables[method].append(shear_estimates_method_table)

                if len(shear_estimates_method_table) > 0:
                    at_least_one_good_estimates_table = True

            except Exception as e:

                logger.warn("Failsafe block encountered exception: " + str(e))
                continue

            logger.debug("Finished loading shear estimates from file: " + shear_estimates_table_product_filename)

    if not at_least_one_good_estimates_table:
        raise RuntimeError("Not able to load any shear estimates tables successfully.")

    logger.info("Finished loading shear estimates from files listed in: " + args.input_shear_estimates_listfile)

    logger.info("Combining shear estimates tables.")

    # Combine the shear estimates tables for each method and output them

    combined_shear_estimates_tables = dict.fromkeys(methods)

    # Create the output product
    combined_shear_estimates_product = products.shear_estimates.create_shear_estimates_product()

    for method in methods:

        # Skip if no data for this method
        if len(shear_estimates_tables[method]) == 0:
            combined_shear_estimates_product.set_method_filename(method, "None")
            continue

        # Combine the tables
        combined_shear_estimates_tables[method] = table.vstack(shear_estimates_tables[method])

        # Get a filename for the table
        combined_shear_estimates_table_filename = get_allowed_filename(type_name="SHEAR-EST-" + method.upper(),
                                                                       instance_id='MERGED',
                                                                       extension=".fits",
                                                                       version=SHE_CTE.__version__,
                                                                       subdir="data",
                                                                       processing_function="SHE")
        combined_shear_estimates_product.set_method_filename(method, combined_shear_estimates_table_filename)

        # Output the table
        combined_shear_estimates_tables[method].write(os.path.join(args.workdir,
                                                                   combined_shear_estimates_table_filename),
                                                      format="fits")

        logger.info("Combined shear estimates for method " + method +
                    " output to: " + combined_shear_estimates_table_filename)

    # Save the product
    write_xml_product(combined_shear_estimates_product, os.path.join(args.workdir,
                                                                     args.output_shear_estimates))

    logger.info("Combined shear estimates product output to: " + args.output_shear_estimates)

    logger.debug('# Exiting shear_estimates_merge_from_args normally')

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

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE 0.7 SHE_CTE_ShearEstimatesMerge",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    try:

        if args.profile:
            import cProfile
            cProfile.runctx("shear_estimates_merge_from_args(args)", {},
                            {"shear_estimates_merge_from_args": shear_estimates_merge_from_args,
                             "args": args, },
                            filename="shear_estimates_merge.prof")
        else:
            shear_estimates_merge_from_args(args)
    except Exception as e:
        # logger.warn("Failsafe exception block triggered with exception: " + str(e))
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
