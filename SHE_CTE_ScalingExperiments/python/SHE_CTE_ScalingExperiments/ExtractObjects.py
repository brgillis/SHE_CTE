""" @file ExtractObjects.py

    Created 25 Aug 2021
    

    Executable for Extracting the objects from the MER tables that we wish to fit.
"""

__updated__ = "2021-08-25"

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

from .extract_objects import extract_objects_from_args



def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ExtractObjects defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # parser.add_argument('--profile', action='store_true',
    #                     help='Store profiling data for execution.')

    # Input data
    parser.add_argument('--pipeline_config', type=str,
                        help="ScalingExperiments config file", required=True)

    parser.add_argument('--stacked_image', type=str,
                        help="stacked image from VIS", required=True)

    parser.add_argument('--exposures', type=str,
                        help="The json pointing to the exposures from VIS,", required=True)

    parser.add_argument('--catalogue_listfile', type=str,
                        help="Listfile pointing to the MER catalogue files", required=True)

    

    # Output data
    parser.add_argument('--output_objects_list', type=str,
                        help='json containing the objects', required=True)

    parser.add_argument('--combined_catalogue', type=str,
                        help='All the MER catalogues joined together', required=True)

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('# Exiting SHE_CTE_ExtractObjects defineSpecificProgramOptions')

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
    logger.debug('# Entering SHE_CTE_ExtractObjects mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_ExtractObjects",
                                    )#store_true=["profile", "debug", "webdav_archive"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    extract_objects_from_args(args)

    logger.debug('# Exiting SHE_CTE_ExtractObjects mainMethod()')


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