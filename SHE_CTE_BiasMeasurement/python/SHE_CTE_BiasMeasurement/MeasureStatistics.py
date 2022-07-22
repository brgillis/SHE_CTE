""" @file MeasureStatistics.py

    Created 7 Apr 2017


    Executable for measuring necessary statistics on a set of shear
    measurements.
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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA

from typing import Any, Dict, Tuple, Type, Union

from SHE_CTE.executor import CteLogOptions, SheCteExecutor
from SHE_PPT.argument_parser import SheArgumentParser
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.executor import ReadConfigArgs
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (CalibrationConfigKeys, ConfigKeys)
from .measure_statistics import measure_statistics_from_args

# Set up dicts for pipeline config defaults and types
D_MS_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    CalibrationConfigKeys.MS_ARCHIVE_DIR   : None,
    CalibrationConfigKeys.MS_WEBDAV_ARCHIVE: False,
    CalibrationConfigKeys.MS_WEBDAV_DIR    : "/mnt/webdav",
    }

D_MS_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    CalibrationConfigKeys.MS_ARCHIVE_DIR   : str,
    CalibrationConfigKeys.MS_WEBDAV_ARCHIVE: bool,
    CalibrationConfigKeys.MS_WEBDAV_DIR    : str,
    }

D_MS_CONFIG_CLINE_ARGS: Dict[ConfigKeys, str] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    CalibrationConfigKeys.MS_ARCHIVE_DIR   : "archive_dir",
    CalibrationConfigKeys.MS_WEBDAV_ARCHIVE: "webdav_dir",
    CalibrationConfigKeys.MS_WEBDAV_DIR    : "webdav_archive",
    }

EXEC_NAME = 'SHE_CTE_EstimateShears'

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
    parser.add_input_arg('--details_table', type = str,
                         help = "Details table data product")
    parser.add_input_arg('--shear_estimates', type = str,
                         help = "Shear estimates data product")

    # Output data
    parser.add_output_arg('--she_bias_statistics', type = str,
                          help = 'Desired name of the output shear bias statistics data product')

    # Archive directory - only default value can be used in pipeline
    parser.add_option_arg('--archive_dir', type = str)

    parser.add_option_arg('--webdav_dir', type = str,
                          help = "Path of the WebDAV mount.")

    parser.add_option_arg('--webdav_archive', action = "store_true",
                          help = "If set, will mount/demount webdav for archiving, and workspace will be relative to " +
                                 "the webdav mount.")

    logger.debug(f'# Exiting {EXEC_NAME} mainMethod()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to measure bias.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    executor = SheCteExecutor(run_from_args_function = measure_statistics_from_args,
                              log_options = CteLogOptions(executable_name = EXEC_NAME),
                              config_args = ReadConfigArgs(d_config_defaults = D_MS_CONFIG_DEFAULTS,
                                                           d_config_types = D_MS_CONFIG_TYPES,
                                                           d_config_cline_args = D_MS_CONFIG_CLINE_ARGS,
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
