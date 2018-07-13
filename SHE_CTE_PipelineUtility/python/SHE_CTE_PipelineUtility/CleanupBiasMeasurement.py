""" @file CleanupBiasMeasurement.py

    Created 5 July 2018

    Main program for cleaning up intermediate files created for the bias measurement pipeline.
"""

__updated__ = "2018-07-13"

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
from SHE_PPT.file_io import read_listfile, read_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string

products.calibrated_frame.init()
products.stacked_frame.init()
products.psf_image.init()
products.mosaic.init()
products.stack_mosaic.init()
products.detections.init()
products.details.init()
products.shear_estimates.init()
products.simulation_config.init()


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_CleanupBiasMeasurement defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Input arguments
    parser.add_argument('--simulation_config', type=str)
    parser.add_argument('--data_images', type=str)
    parser.add_argument('--stacked_data_image', type=str)
    parser.add_argument('--psf_images_and_tables', type=str)
    parser.add_argument('--segmentation_images', type=str)
    parser.add_argument('--stacked_segmentation_image', type=str)
    parser.add_argument('--detections_tables', type=str)
    parser.add_argument('--details_table', type=str)
    parser.add_argument('--shear_estimates', type=str)

    parser.add_argument('--shear_bias_statistics_in', type=str)  # Needed to ensure it waits until ready

    # Output arguments
    parser.add_argument('--shear_bias_statistics_out', type=str)

    # Required pipeline arguments
    parser.add_argument('--workdir', type=str,)
    parser.add_argument('--logdir', type=str,)
    parser.add_argument('--debug', action='store_true',
                        help="Set to enable debugging protocols")

    logger.debug('# Exiting SHE_Pipeline_Run defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_CleanupBiasMeasurement mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE 0.5 SHE_CTE_CleanupBiasMeasurement",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    def remove_file(qualified_filename):
        if not os.path.exists(qualified_filename):
            logger.warn("Expected file '" + qualified_filename + "' does not exist.")
            return
        os.remove(qualified_filename)
        return

    def remove_product(qualified_filename):

        # Load the product
        if not os.path.exists(qualified_filename):
            logger.warn("Expected file '" + qualified_filename + "' does not exist.")
            return
        p = read_xml_product(qualified_filename)

        # Remove all files this points to
        if hasattr(p, "get_all_filenames"):
            data_filenames = p.get_all_filenames()
            for data_filename in data_filenames:
                if data_filename is not None:
                    remove_file(os.path.join(args.workdir, data_filename))
        else:
            logger.error("Product " + qualified_filename + " has no 'get_all_filenames' method.")
            if not args.debug:
                raise RuntimeError("Product " + qualified_filename + " has no 'get_all_filenames' method.")

        # Remove the product itself
        remove_file(qualified_filename)

    # Clean up all files pointed to by listfiles and the listfiles themselves
    for listfile_name in (args.data_images,
                          args.psf_images_and_tables,
                          args.segmentation_images,
                          args.detections_tables):

        qualified_listfile_name = os.path.join(args.workdir, listfile_name)

        if not os.path.exists(qualified_listfile_name):
            logger.warn("Expected file '" + qualified_listfile_name + "' does not exist.")
            continue

        filenames = read_listfile(qualified_listfile_name)

        for filename in filenames:
            remove_product(os.path.join(args.workdir, filename))

        remove_file(qualified_listfile_name)

    # Now remove all products
    for product_filename in (args.simulation_config,
                             args.stacked_data_image,
                             args.stacked_segmentation_image,
                             args.details_table,
                             args.shear_estimates):
        remove_product(os.path.join(args.workdir, product_filename))

    # Move the statistics product to the new name
    qualified_stats_in_filename = os.path.join(args.workdir, args.shear_bias_statistics_in)
    qualified_stats_out_filename = os.path.join(args.workdir, args.shear_bias_statistics_out)

    os.rename(qualified_stats_in_filename, qualified_stats_out_filename)

    logger.debug('# Exiting SHE_CTE_CleanupBiasMeasurement mainMethod()')

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
