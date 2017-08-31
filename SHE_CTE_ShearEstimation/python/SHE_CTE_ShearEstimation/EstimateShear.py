""" @file EstimateShears.py

    Created 27 Mar 2017

    Main program for estimating shears on simulation data.

    ---------------------------------------------------------------------

    Copyright (C) 2012-2020 Euclid Science Ground Segment      
       
    This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General    
    Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)    
    any later version.    
       
    This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied    
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more    
    details.    
       
    You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to    
    the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
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

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_EstimateShear defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile',action='store_true',
                        help='Store profiling data for execution.')
    
    # Input arguments
    
    parser.add_argument('--data_image',type='str',
                        help='FITS file containing the observed/simulated image.')
    
    parser.add_argument('--psf_image',type='str',
                        help='FITS file containing the PSF image.')
    
    parser.add_argument('--noise_image',type='str',
                        help='FITS file containing the noise map.')
    
    parser.add_argument('--mask_image',type='str',
                        help='FITS file containing the mask.')
    
    parser.add_argument('--segmentation_image',type='str',
                        help='FITS file containing the segmentation map.')
    
    parser.add_argument('--detections_table',type='str',
                        help='FITS file containing the detections table.')
    
    # Output arguments
    
    parser.add_argument('--shear_measurements_table',type='str',
                        help='FITS file to contain the shear measurements table.')
    
    parser.add_argument('--mcmc_chain_table',type='str',
                        help='FITS file to contain the MCMC chains table.')

    logger.debug('# Exiting SHE_CTE_EstimateShear defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_EstimateShear mainMethod()')
    logger.debug('#')
        
    if args.profile:
        import cProfile
        cProfile.runctx("estimate_shears_from_args(vars(args))",{},
                        {"estimate_shears_from_args":estimate_shears_from_args,
                         "args":args},filename="measure_shapes.prof")
    else:
        estimate_shears_from_args(vars(args))

    logger.debug('# Exiting SHE_CTE_EstimateShear mainMethod()')

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