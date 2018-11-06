""" @file MeasureStatistics.py

    Created 7 Apr 2017
    

    Executable for measuring necessary statistics on a set of shear
    measurements.
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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA

import argparse

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_BiasMeasurement.measure_statistics import measure_statistics_from_args
from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureStatistics defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')

    # Input data
    parser.add_argument('--details_table', type=str,
                        help="Details table data product")
    parser.add_argument('--shear_estimates', type=str,
                        help="Shear estimates data product")

    parser.add_argument("--pipeline_config", default=None, type=str,
                        help="Pipeline-wide configuration file.")

    # Output data
    parser.add_argument('--shear_bias_statistics', type=str,
                        help='Desired name of the output shear bias statistics data product')

    # Archive directory - only default value can be used in pipeline
    parser.add_argument('--archive_dir', type=str, default=None)
    
    parser.add_argument('--webdav_dir', type=str, default="/mnt/webdav",
                        help="Path of the WebDAV mount.")

    parser.add_argument('--webdav_archive', action="store_true",
                        help="If set, will mount/demount webdav for archiving, and workspace will be relative to " +
                        "the webdav mount.")

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('# Exiting SHE_CTE_MeasureStatistics mainMethod()')

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
    logger.debug('# Entering SHE_CTE_EstimateShears mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE 0.6.9 SHE_CTE_MeasureStatistics",
                                    store_true=["profile", "debug", "webdav_archive"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    if args.profile:
        import cProfile
        cProfile.runctx("measure_statistics_from_args(args)", {},
                        {"measure_statistics_from_args": measure_statistics_from_args,
                         "args": args}, filename="measure_statistics.prof")
    else:
        measure_statistics_from_args(args)

    logger.debug('# Exiting MeasureStatistics mainMethod()')

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
