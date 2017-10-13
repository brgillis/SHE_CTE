""" @file MakeMockAnalysisData.py

    Created 12 Oct 2017

    Executable for making mock data for the analysis pipeline.

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

from ElementsKernel.Logging import getLogger

from SHE_CTE_Pipeline import magic_values as mv
from SHE_CTE.magic_values import force_dry_run

if force_dry_run:
    from SHE_CTE_Pipeline.mock_analysis_data_dry import make_mock_analysis_data
else:
    from SHE_CTE_Pipeline.mock_analysis_data import make_mock_analysis_data

def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program, using all possible configurations.

    @return
        An  ArgumentParser.
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MakeMockAnalysisData defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Option for profiling
    parser.add_argument('--profile',action='store_true',
                        help='Store profiling data for execution.')
    
    # Setup options (Cannot be used in pipeline)
    parser.add_argument('--num_exposures', type=int, default=4,
                        help='Number of exposures.')
    parser.add_argument('--num_detectors', type=int, default=36,
                        help='Number of detectors per exposure.')
    
    # Output filenames
    parser.add_argument('--data_images', type=str, default="mock_data_images.json",
                        help="Desired name of listfile containing mock data images.")
    
    parser.add_argument('--psf_calibration_products', type=str, default="mock_psf_calibration_products.json",
                        help="Desired name of listfile containing mock PSF calibration products.")
    
    parser.add_argument('--segmentation_images', type=str, default="mock_segmentation_images.json",
                        help="Desired name of listfile containing mock segmentation images.")
    
    parser.add_argument('--detections_tables', type=str, default="mock_detections_tables.json",
                        help="Desired name of listfile containing mock detections tables.")
    
    parser.add_argument('--astrometry_products', type=str, default="mock_astrometry_products.json",
                        help="Desired name of listfile containing mock astrometry products.")
    
    parser.add_argument('--aocs_time_series_products', type=str, default="mock_aocs_time_series_products.json",
                        help="Desired name of listfile containing mock aocs time series products.")
    
    parser.add_argument('--mission_time_products', type=str, default="mock_mission_time_products.json",
                        help="Desired name of listfile containing mock mission time products.")
    
    parser.add_argument('--galaxy_population_priors_table', type=str, default="mock_galaxy_population_priors_table.fits",
                        help="Desired name of fits file containing mock galaxy_population_priors.")
    
    parser.add_argument('--calibration_parameters_product', type=str, default="mock_calibration_parameters_product.bin",
                        help="Desired name of mock calibration parameters data product.")
    
    parser.add_argument('--calibration_parameters_listfile', type=str, default="mock_calibration_parameters_listfile.json",
                        help="Desired name of listfile containing files associated with the mock calibration parameters product.")
    
    parser.add_argument('--shear_validation_statistics_table', type=str, default="mock_shear_validation_statistics_table.fits",
                        help="Desired name of fits file containing mock shear validation statistics.")
    

    logger.debug('Exiting SHE_CTE_MakeMockAnalysisData defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to generate galaxy images.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MakeMockAnalysisData mainMethod()')
    logger.debug('#')
        
    if args.profile:
        import cProfile
        cProfile.runctx("make_mock_analysis_data(args)",{},
                        {"make_mock_analysis_data":make_mock_analysis_data,
                         "args":args,},
                        filename="make_mock_analysis_data.prof")
    else:
        make_mock_analysis_data(args)

    logger.debug('Exiting SHE_CTE_MakeMockAnalysisData mainMethod()')

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