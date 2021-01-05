""" @file DrawSensTypes.py

    Created 5 January 2021

    Main function to draw a plot explaining different types of sensitivity testing.
"""

__updated__ = "2021-01-05"

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
from os.path import join

from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string

import SHE_CTE
import matplotlib.pyplot as pyplot
import numpy as np


# TODO: Define constants here


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_DrawSensTypes defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')

    # Output data
    parser.add_argument('--sens_types_plot', type=str, default="sens_types.pdf",
                        help='Desired name of the output sensitivity types plot.')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('# Exiting SHE_CTE_DrawSensTypes defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_DrawSensTypes mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_DrawSensTypes",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    if args.profile:
        import cProfile
        cProfile.runctx("draw_sens_types_from_args(args)", {},
                        {"draw_sens_types_from_args": draw_sens_types_from_args,
                         "args": args}, filename="measure_bias.prof")
    else:
        draw_sens_types_from_args(args)

    logger.debug('# Exiting SHE_CTE_DrawSensTypes mainMethod()')

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


def draw_sens_types_from_args(args):
    """ @TODO main docstring
    """

    return
