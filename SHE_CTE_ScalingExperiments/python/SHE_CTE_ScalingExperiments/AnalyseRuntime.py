""" @file AnalyseRuntime.py

    Created 27 Aug 2021
    

    Executable that analyses the runtime.
"""

__updated__ = "2021-08-26"

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

import argparse
import os

from EL_PythonUtils.utilities import get_arguments_string
# from SHE_PPT.constants.config import D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES, D_GLOBAL_CONFIG_CLINE_ARGS
from SHE_PPT.logging import getLogger
# from SHE_PPT.pipeline_utility import (read_calibration_config, CalibrationConfigKeys, ConfigKeys, GlobalConfigKeys)

import SHE_CTE

from .analyse_runtime import analyse_runtime_from_args



def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_AnalyseRuntime defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # parser.add_argument('--profile', action='store_true',
    #                     help='Store profiling data for execution.')

    # Input data
    parser.add_argument('--pipeline_config', type=str,
                        help="ScalingExperiments config file", required=True)

    parser.add_argument('--stamp_timing_listfile', type=str,
                        help="listfile for the timing information for stamp extraction", required=True)

    parser.add_argument('--split_timing_listfile', type=str,
                        help="listfile for the timing information for splitting the input fits", required=True)

    # Output data
    parser.add_argument('--results', type=str,
                        help='json containing the timing results and links to plots', required=True)

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('# Exiting SHE_CTE_AnalyseRuntime defineSpecificProgramOptions')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_AnalyseRuntime mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_AnalyseRuntime",
                                    )#store_true=["profile", "debug", "webdav_archive"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    analyse_runtime_from_args(args)

    logger.debug('# Exiting SHE_CTE_AnalyseRuntime mainMethod()')


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