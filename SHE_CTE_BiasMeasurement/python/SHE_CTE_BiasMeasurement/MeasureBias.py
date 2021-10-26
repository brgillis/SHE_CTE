""" @file MeasureBias.py

    Created 7 Apr 2017

    Main program for measuring bias of shear estimates.
"""

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
from typing import Any, Dict, Tuple, Type, Union

import SHE_CTE
from EL_PythonUtils.utilities import get_arguments_string
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (CalibrationConfigKeys, ConfigKeys, GlobalConfigKeys, read_calibration_config)
from .measure_bias import measure_bias_from_args

# Set up dicts for pipeline config defaults and types
D_MB_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    CalibrationConfigKeys.MB_ARCHIVE_DIR   : None,
    CalibrationConfigKeys.MB_NUM_THREADS   : 8,
    CalibrationConfigKeys.MB_WEBDAV_ARCHIVE: False,
    CalibrationConfigKeys.MB_WEBDAV_DIR    : None,
    }

D_MB_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    CalibrationConfigKeys.MB_ARCHIVE_DIR   : str,
    CalibrationConfigKeys.MB_NUM_THREADS   : int,
    CalibrationConfigKeys.MB_WEBDAV_ARCHIVE: bool,
    CalibrationConfigKeys.MB_WEBDAV_DIR    : str,
    }

D_MB_CONFIG_CLINE_ARGS: Dict[ConfigKeys, str] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    CalibrationConfigKeys.MB_ARCHIVE_DIR   : "archive_dir",
    CalibrationConfigKeys.MB_NUM_THREADS   : "number_threads",
    CalibrationConfigKeys.MB_WEBDAV_ARCHIVE: "webdav_dir",
    CalibrationConfigKeys.MB_WEBDAV_DIR    : "webdav_archive",
    }


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureBias defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile', action = 'store_true',
                        help = 'Store profiling data for execution.')

    # Input data
    parser.add_argument('--she_bias_statistics', type = str,
                        help = 'Listfile pointing to shear bias statistics objects.')
    parser.add_argument('--bootstrap_seed', type = int, default = 0,
                        help = 'Seed for bootstrapping of errors')

    parser.add_argument("--pipeline_config", default = None, type = str,
                        help = "Pipeline-wide configuration file.")

    # Output data
    parser.add_argument('--she_bias_measurements', type = str,
                        help = 'Desired name of the output shear bias statistics data product')
    parser.add_argument('--details_table_head', type = str, default = None,
                        help = 'Desired head for the filenames of the output details tables')

    # Input arguments
    parser.add_argument('--number_threads', type = int, default = None,
                        help = 'Number of parallel threads to use.')
    parser.add_argument('--recovery_bias_statistics_filename', type = str, default = "she_bias_statistics.xml",
                        help = 'Expected name of bias statistics files for when operating in recovery mode')
    parser.add_argument('--recovery_bias_measurements_filename', type = str, default = "she_bias_measurements.xml",
                        help = 'Expected name of bias measurements files for when operating in recovery mode')
    parser.add_argument('--store_measurements_only', action = 'store_true',
                        help = 'If set, the resulting bias measurements file will contain only the measurements. ' +
                               'Otherwise, it will also store data for all individual galaxy measurements.')
    parser.add_argument('--use_bias_only', action = 'store_true',
                        help = 'If set will calculate using only bias measurements and not full statistics.')

    # Archive directory
    parser.add_argument('--archive_dir', type = str, default = None)

    parser.add_argument('--webdav_dir', type = str, default = "/mnt/webdav",
                        help = "Path of the WebDAV mount.")

    parser.add_argument('--webdav_archive', action = "store_true",
                        help = "If set, will mount/demount webdav for archiving, and workspace will be relative to " +
                               "the webdav mount.")

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type = str, default = ".")
    parser.add_argument('--logdir', type = str, default = ".")

    logger.debug('# Exiting SHE_CTE_MeasureBias defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to measure bias.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureBias mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd = "E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_MeasureBias",
                                    store_true = ["profile", "debug", "webdav_archive", "store_measurements_only"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    # load the pipeline config in
    args.pipeline_config = read_calibration_config(args.pipeline_config,
                                                   workdir = args.workdir,
                                                   d_defaults = D_MB_CONFIG_DEFAULTS,
                                                   d_cline_args = D_MB_CONFIG_CLINE_ARGS,
                                                   parsed_args = args,
                                                   d_types = D_MB_CONFIG_TYPES)

    # check if profiling is to be enabled from the pipeline config
    profiling = args.pipeline_config[GlobalConfigKeys.PIP_PROFILE]

    if profiling:

        import cProfile
        logger.info("Profiling enabled")

        filename = os.path.join(args.workdir, args.logdir, "measure_bias.prof")
        logger.info("Writing profiling data to %s", filename)

        cProfile.runctx("measure_bias_from_args(args)", {},
                        {"measure_bias_from_args": measure_bias_from_args,
                         "args"                  : args}, filename = filename)
    else:
        logger.info("Profiling disabled")
        measure_bias_from_args(args)

    logger.debug('# Exiting SHE_CTE_MeasureBias mainMethod()')


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
