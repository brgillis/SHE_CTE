""" @file CleanupBiasMeasurement.py

    Created 5 July 2018

    Main program for cleaning up intermediate files created for the bias measurement pipeline.
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

import os
import shutil
from typing import Any, Dict, Tuple, Type, Union

from SHE_CTE.executor import CteLogOptions, SheCteExecutor
from SHE_PPT.argument_parser import SheArgumentParser
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.executor import ReadConfigArgs
from SHE_PPT.file_io import read_listfile, read_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (CalibrationConfigKeys, ConfigKeys)

# Set up dicts for pipeline config defaults and types
D_CBM_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    CalibrationConfigKeys.CBM_CLEANUP: True,
    }

D_CBM_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    CalibrationConfigKeys.CBM_CLEANUP: bool,
    }

D_CBM_CONFIG_CLINE_ARGS: Dict[ConfigKeys, str] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    }

EXEC_NAME = "SHE_CTE_CleanupBiasMeasurement"

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

    # Input arguments
    parser.add_input_arg('--simulation_config', type = str)
    parser.add_input_arg('--data_images', type = str)
    parser.add_input_arg('--stacked_data_image', type = str)
    parser.add_input_arg('--psf_images_and_tables', type = str)
    parser.add_input_arg('--segmentation_images', type = str)
    parser.add_input_arg('--stacked_segmentation_image', type = str)
    parser.add_input_arg('--detections_tables', type = str)
    parser.add_input_arg('--details_table', type = str)
    parser.add_input_arg('--shear_estimates', type = str)

    parser.add_input_arg('--shear_bias_statistics_in', type = str)  # Needed to ensure it waits until ready

    # Output arguments
    parser.add_output_arg('--shear_bias_statistics_out', type = str)

    logger.debug(f'# Exiting {EXEC_NAME} defineSpecificProgramOptions()')

    return parser


def cleanup_bias_measurement_from_args(args):
    """ Function which performs the heavy lifting of actually cleaning up unneeded intermediate files
        in the bias measurement pipeline.
    """

    qualified_stats_in_filename = os.path.join(args.workdir, args.shear_bias_statistics_in)
    qualified_stats_out_filename = os.path.join(args.workdir, args.shear_bias_statistics_out)

    # Read in the pipeline config, which tells us whether to clean up or not
    if args.pipeline_config is None:
        logger.warning("No pipeline configuration found. Being safe and not cleaning up.")
        return
    pipeline_config = args.pipeline_config

    # Check for the cleanup key
    clean_up = pipeline_config[CalibrationConfigKeys.CBM_CLEANUP]
    if not clean_up:
        logger.info("Config is set to not clean up, so ending execution.")

        # Copy the statistics product to the new name
        shutil.copy(qualified_stats_in_filename, qualified_stats_out_filename)
        return

    # Move the statistics product to the new name
    os.rename(qualified_stats_in_filename, qualified_stats_out_filename)

    # If we get here, we are cleaning up

    def remove_file(qualified_filename):
        if qualified_filename[-1] == "/":
            logger.warning("Attempted to remove directory " + qualified_filename)
            return 1
        if not os.path.exists(qualified_filename):
            logger.warning(f"Expected file '{qualified_filename}' does not exist.")
            return 1
        os.remove(qualified_filename)
        return 0

    def remove_product(qualified_filename):

        if qualified_filename[-1] == "/":
            logger.warning("Attempted to remove directory " + qualified_filename)
            return

        # Load the product
        if not os.path.exists(qualified_filename):
            logger.warning("Expected file '" + qualified_filename + "' does not exist.")
            return
        p = read_xml_product(qualified_filename)

        # Remove all files this points to
        if hasattr(p, "get_all_filenames"):
            data_filenames = p.get_all_filenames()
            for data_filename in data_filenames:
                if (data_filename is not None and data_filename != "default_filename.fits" and
                        data_filename != "" and data_filename != "None"):
                    if remove_file(os.path.join(args.workdir, data_filename)):
                        logger.warning("...in removal of product " + qualified_filename)
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
            logger.warning("Expected file '" + qualified_listfile_name + "' does not exist.")
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


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    executor = SheCteExecutor(run_from_args_function = cleanup_bias_measurement_from_args,
                              log_options = CteLogOptions(executable_name = EXEC_NAME),
                              config_args = ReadConfigArgs(d_config_defaults = D_CBM_CONFIG_DEFAULTS,
                                                           d_config_types = D_CBM_CONFIG_TYPES,
                                                           d_config_cline_args = D_CBM_CONFIG_CLINE_ARGS,
                                                           s_config_keys_types = {CalibrationConfigKeys},
                                                           ))

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
