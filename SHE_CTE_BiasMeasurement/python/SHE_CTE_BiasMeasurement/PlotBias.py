""" @file PlotBias.py

    Created 19 July 2018

    Main program for plotting bias sensitivity
"""

__updated__ = "2018-08-31"

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

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_BiasMeasurement.plot_bias_measurements import plot_bias_measurements_from_args
from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string
import numpy as np


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_PlotBias defineSpecificProgramOptions()')
    logger.debug('#')

    # Set up the command-line arguments
    parser = argparse.ArgumentParser()

    # Input Data
    parser.add_argument('--methods', nargs='*', default=["KSB", "REGAUSS", "MomentsML", "LensMC"],
                        help='Methods to plot bias measurements for.')
    parser.add_argument('--root_data_folder', default="../",
                        help="Folder containing folders which contain bias measurements. Must be either absolute " +
                             "or relative to workdir.")
    parser.add_argument('--data_folder_head', default="workspace_",
                        help="Head of folder names containing bias measurements, minus the E??P??S?? portion.")
    parser.add_argument('--output_file_name_head', default="sensitivity_testing",
                        help="Desired head of filenames for created files.")
    parser.add_argument('--output_format', default="png",
                        help="File format (and extension) of created images.")
    parser.add_argument('--hide', action="store_true",
                        help="If set, will not display plots when they are generated.")

    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default="/home/user/Work/workspace")
    parser.add_argument('--logdir', type=str, default="logs")

    logger.debug('# Exiting SHE_CTE_PlotBias defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_PlotBias mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE 0.6.9 SHE_CTE_PlotBias",
                                    store_true=["profile", "debug", "hide"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    if args.profile:
        import cProfile
        cProfile.runctx("plot_bias_measurements_from_args(args)", {},
                        {"plot_bias_measurements_from_args": plot_bias_measurements_from_args,
                         "args": args, },
                        filename="plot_bias_measurements.prof")
    else:
        plot_bias_measurements_from_args(args)

    logger.debug('# Exiting SHE_CTE_PlotBias mainMethod()')

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
