""" @file EstimateShears.py

    Created 27 Mar 2017

    Main program for estimating shears on simulation data.

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

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_CTE_ShearEstimation.estimate_shears import estimate_shears_from_args

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
    parser.add_argument('detections_table_file_name',type=str,
                        help='Name of the data table containing positions of detected galaxies')
    parser.add_argument('galaxies_image_file_name',type=str,
                        help='Name of the simulated image file to process.')
    parser.add_argument('psf_image_file_name',type=str,
                        help='Name of the file containing the PSF image.')
    
    # Info on noise
    parser.add_argument('--gain', type=float,
                        help="Gain of the galaxy image, in e-/ADU")
    parser.add_argument('--read_noise', type=float,
                        help="Read noise of the galaxy image, in e-/pixel")
    parser.add_argument('--subtracted_sky_level', type=float,
                        help="Sky level that's been subtracted from the galaxy image, in ADU/arcsec^2")
    parser.add_argument('--p_of_e_table_file_name',type=str,
                        help='Name of the file containing the P(e) table.')
    
    # Method to use
    parser.add_argument('--method',type=str,default='KSB',
                        help='Shape measurement method to use. Allowed: \'KSB\' (default)')
    
    # Output data
    parser.add_argument('--output_file_name',type=str,default=None,
                        help='Desired name of the output fits file containing shape measurements.')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to estimate shears.

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
        cProfile.runctx("estimate_shears_from_args(vars(args))",{},
                        {"estimate_shears_from_args":estimate_shears_from_args,
                         "args":args},filename="measure_shapes.prof")
    else:
        estimate_shears_from_args(vars(args))

    logger.debug('# Exiting EstimateShear mainMethod()')

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