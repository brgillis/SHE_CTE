""" @file MeasureStatistics.py

    Created 7 Apr 2017

    Executable for measuring necessary statistics on a set of shear
    measurements.

    ---------------------------------------------------------------------

    Copyright (C) 2017 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
from SHE_GST_IceBRGpy.logging import getLogger

from SHE_CTE_BiasMeasurement import magic_values as mv

def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureStatistics defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile',action='store_true',
                        help='Store profiling data for execution.')
    
    # Input data
    parser.add_argument('input_cat', type=str,
                        help="Filename of the input file catalogue, which should be in " +
                        "ascii.ecsv format.")
    
    parser.add_argument('--Pe', type=str, default="p0",
                        help='P(e) to use. Must be one of m2, m1, p0, p1, or p2.')
    
    # Output data
    parser.add_argument('--output_cat',type=str,
                        help='Desired name of the output file catalogue.')

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

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_EstimateShears mainMethod()')
    logger.debug('#')
        
    if args.profile:
        import cProfile
        cProfile.runctx("measure_statistics_from_args(vars(args))",{},
                        {"measure_statistics_from_args":measure_statistics_from_args,
                         "args":args},filename="measure_statistics.prof")
    else:
        measure_statistics_from_args(vars(args))

    logger.debug('# Exiting MeasureBias mainMethod()')

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