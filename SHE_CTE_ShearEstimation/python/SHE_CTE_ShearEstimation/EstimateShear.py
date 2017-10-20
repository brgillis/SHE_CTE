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
from SHE_CTE.magic_values import force_dry_run
from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_CTE_ShearEstimation.estimate_shears import estimate_shears_from_args

from SHE_GST_IceBRGpy.logging import getLogger

def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_EstimateShear defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile',action='store_true',
                        help='Store profiling data for execution.')
    parser.add_argument('--dry_run',action='store_true',
                        help='Dry run (no data processed).')
    
    # Required input arguments
    
    parser.add_argument('--data_images',type=str,
                        help='.json listfile containing filenames of data images.')
    
    parser.add_argument('--psf_images_and_tables',type=str,
                        help='.json listfile containing filenames of psf images.')
    
    parser.add_argument('--segmentation_images',type=str,
                        help='.json listfile containing filenames of segmentation map images.')
    
    parser.add_argument('--detections_tables',type=str,
                        help='.json listfile containing filenames of detections tables.')
    
    parser.add_argument('--galaxy_population_priors_table',type=str,
                        help='Filename of galaxy population priors table (fits binary table).')
    
    parser.add_argument('--calibration_parameters_product',type=str,
                        help='Filename of calibration parameters product (XML data product).')
    
    parser.add_argument('--calibration_parameters_listfile',type=str,
                        help='.json listfile containing filenames of calibration parameters files.')
    
    # Optional input arguments (cannot be used in pipeline)
    
    parser.add_argument('--methods',type=str, nargs='*', default=[],
                        help='Which shear estimation methods to apply. If not specified, all will be run.')
    
    # Output arguments
    
    parser.add_argument('--shear_estimates_product',type=str,
                        help='XML data product to contain file links to the shear estimates tables.')
    parser.add_argument('--shear_estimates_listfile',type=str,
                        help='.json listfile to contain filenames of shear estimates tables.')
    
    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir',type=str,default=".")
    parser.add_argument('--logdir',type=str,default=".")

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
    
    dry_run = args.dry_run or force_dry_run
        
    if args.profile:
        import cProfile
        cProfile.runctx("estimate_shears_from_args(args)",{},
                        {"estimate_shears_from_args":estimate_shears_from_args,
                         "args":args,
                         dry_run:dry_run},
                        filename="measure_shapes.prof")
    else:
        estimate_shears_from_args(args,dry_run)

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