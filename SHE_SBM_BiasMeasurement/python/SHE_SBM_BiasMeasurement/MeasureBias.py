""" @file MeasureBias.py

    Created 7 Apr 2017

    Main program for measuring bias of shear estimates.

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
from icebrgpy.logging import getLogger

from SHE_SBM_BiasMeasurement import magic_values as mv
from SHE_SBM_BiasMeasurement.measure_bias import measure_bias_from_args

def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile',action='store_true',
                        help='Store profiling data for execution.')
    
    # Input data
    parser.add_argument('input_dir',type=str,
                        help='Path of the directory containing shear measurements and details tables.')
    parser.add_argument('--required_input_pattern', type=str, default=None,
                        help='Required pattern in file names for measurements or details tables.')
    parser.add_argument('--input_depth', type=int, default=0,
                        help='Maximum subfolder depth to search for input files. 0 = input_dir only.')
    
    parser.add_argument('--e', type=str, default="p0",
                        help='P(e) to use. Must be one of m2, m1, p0, p1, or p2.')
    
    # Output data
    parser.add_argument('--output_file_name',type=str,default=mv.default_output_filename,
                        help='Desired name of the output table containing bias measurements.')
    parser.add_argument('--output_format',type=str,default=mv.default_output_format,
                        help='Desired format of the output table.')

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
    logger.debug('# Entering SHE_SBM_EstimateShears mainMethod()')
    logger.debug('#')
        
    if args.profile:
        import cProfile
        cProfile.runctx("measure_bias_from_args(vars(args))",{},
                        {"measure_bias_from_args":measure_bias_from_args,
                         "args":args},filename="measure_bias.prof")
    else:
        measure_bias_from_args(vars(args))

    logger.debug('# Exiting MeasureBias mainMethod()')

    return