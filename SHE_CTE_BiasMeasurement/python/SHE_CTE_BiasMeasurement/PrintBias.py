""" @file PrintBias.py

    Created 16 July 2018

    Main program for printing out bias of shear estimates
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

import argparse
import os

from EL_PythonUtils.utilities import get_arguments_string
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES, D_GLOBAL_CONFIG_CLINE_ARGS
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import read_config, CalibrationConfigKeys, GlobalConfigKeys

import SHE_CTE
from SHE_CTE_BiasMeasurement.print_bias import print_bias_from_product_filename


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_PrintBias defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')

    # Input data
    parser.add_argument('--she_bias_measurements', type=str,
                        help='Desired name of the output shear bias statistics data product')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('# Exiting SHE_CTE_PrintBias defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_PrintBias mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_PrintBias",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    # load the pipeline config in
    args.pipeline_config = read_config(args.pipeline_config,
                                       workdir=args.workdir,
                                       config_keys=CalibrationConfigKeys,
                                       defaults=D_GLOBAL_CONFIG_DEFAULTS,
                                       d_cline_args=D_GLOBAL_CONFIG_CLINE_ARGS,
                                       parsed_args=args,
                                       d_types=D_GLOBAL_CONFIG_TYPES)

    # check if profiling is to be enabled from the pipeline config
    profiling = args.pipeline_config[GlobalConfigKeys.PIP_PROFILE]

    if profiling:
        import cProfile
        logger.info("Profiling enabled")

        filename = os.path.join(args.workdir, args.logdir, "print_bias.prof")
        logger.info("Writing profiling data to %s", filename)

        cProfile.runctx("print_bias_from_product_filename(bias_measurements, workdir)", {},
                        {"print_bias_from_product_filename": print_bias_from_product_filename,
                         "bias_measurements": args.she_bias_measurements,
                         "workdir": args.workdir, },
                        filename=filename)
    else:
        logger.info("Profiling disabled")
        print_bias_from_product_filename(args.she_bias_measurements, args.workdir)

    logger.debug('# Exiting SHE_CTE_PrintBias mainMethod()')


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
