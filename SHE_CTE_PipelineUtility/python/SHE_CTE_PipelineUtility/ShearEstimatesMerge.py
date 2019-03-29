""" @file ObjectIdMerge.py

    Created 14 Mar 2019

    Merge point executable for merging up batches into objects.
"""

__updated__ = "2019-03-29"

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
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_xml_product, write_xml_product,
                             get_allowed_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.shear_estimates import tf as setf
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import get_arguments_string
from astropy import table
from astropy.table import Table


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
    logger.debug('# Entering SHE_CTE_ObjectIdMerge defineSpecificProgramOptions()')
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

    logger.debug('# Exiting SHE_CTE_ObjectIdMerge defineSpecificProgramOptions()')

    return parser


def object_id_merge_from_args(args):
    """ Core function for implementing a merge of object ID
    """

    logger.debug('# Entering object_id_merge_from_args(args)')

    # Keep a list of all shear estimates tables for each method
    shear_estimates_tables = dict.fromkeys(methods)
    for method in methods:
        # Start with an empty list of the tables
        shear_estimates_tables[method] = []

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

            if not isinstance(shear_estimates_table_product, shear_estimates.dpdShearMeasurement):
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

                shear_estimates_method_table = Table.read(os.path.join(
                    args.workdir, shear_estimates_method_table_filename))

                if not is_in_format(shear_estimates_method_table, setf, verbose=True, ignore_metadata=True):
                    raise TypeError("Input shear estimates table for method {} is of invalid format.".format(method))

                # Get the ID list from it and add it to the set
                shear_estimates_tables[method].append(shear_estimates_method_table)

                if len(shear_estimates_method_table) > 0:
                    at_least_one_good_estimates_table = True

            except Exception as e:

                logger.warn("Failsafe block encountered exception: " + str(e))
                continue

    if not at_least_one_good_estimates_table:
        raise RuntimeError("Not able to load any shear estimates tables successfully.")

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
                                                                       release="00.07",
                                                                       subdir="data",
                                                                       processing_function="SHE")
        combined_shear_estimates_product.set_method_filename(method, combined_shear_estimates_table_filename)

        # Output the table
        combined_shear_estimates_tables[method].write(os.path.join(args.workdir,
                                                                   combined_shear_estimates_table_filename),
                                                      format="fits")

    # Get a filename for the output product
    combined_shear_estimates_product_filename = get_allowed_filename(type_name="SHEAR-EST",
                                                                     instance_id='MERGED',
                                                                     extension=".xml",
                                                                     release="00.07",
                                                                     subdir="data",
                                                                     processing_function="SHE")

    # Save the product
    write_xml_product(combined_shear_estimates_product, os.path.join(args.workdir,
                                                                     combined_shear_estimates_product_filename))

    logger.debug('# Exiting object_id_merge_from_args normally')

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
    logger.debug('# Entering SHE_CTE_ObjectIdMerge mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE 0.7 SHE_CTE_ObjectIdMerge",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    try:

        if args.profile:
            import cProfile
            cProfile.runctx("object_id_merge_from_args(args)", {},
                            {"object_id_merge_from_args": object_id_merge_from_args,
                             "args": args, },
                            filename="object_id_merge.prof")
        else:
            object_id_merge_from_args(args)
    except Exception as e:
        # logger.warn("Failsafe exception block triggered with exception: " + str(e))
        raise

    logger.debug('# Exiting SHE_CTE_ObjectIdMerge mainMethod()')

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
