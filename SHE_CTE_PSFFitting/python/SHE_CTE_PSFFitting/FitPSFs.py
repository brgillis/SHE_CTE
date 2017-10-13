""" @file FitPSFs.py

    Created 12 Oct 2017

    Executable for fitting PSFs and generating images of them.

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

from ElementsKernel.Loggingimport getLogger

from SHE_CTE_PSFFitting import magic_values as mv
from SHE_CTE.magic_values import force_dry_run

if force_dry_run:
    from SHE_CTE_PSFFitting.fit_psfs_dry import fit_psfs
else:
    from SHE_PSM import fit_psfs

def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program, using all possible configurations.

    @return
        An  ArgumentParser.
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_FitPSFs defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Option for profiling
    parser.add_argument('--profile',action='store_true',
                        help='Store profiling data for execution.')
    
    # Input filenames
    parser.add_argument('data_images',type=str,
                        help='Filename for list of data images (.json listfile)')
    parser.add_argument('detections_tables',type=str,
                        help='Filename for list of detections tables (.json listfile)')
    parser.add_argument('astrometry_products',type=str,
                        help='Filename for list of astrometry data products (.json listfile)')
    parser.add_argument('aocs_time_series_products',type=str,
                        help='Filename for list of AOCS data series data products (.json listfile)')
    parser.add_argument('mission_time_products',type=str,
                        help='Filename for list of mission time data products (.json listfile)')
    parser.add_argument('psf_calibration_products',type=str,
                        help='Filename for list of PSF calibration data products (.json listfile)')
    
    # Output filenames
    parser.add_argument('psf_images_and_tables',type=str,
                        help='Desired filename for output list of PSF images and tables (.json listfile).')
    

    logger.debug('Exiting SHE_CTE_FitPSFs defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_FitPSFs mainMethod()')
    logger.debug('#')
        
    if args.profile:
        import cProfile
        cProfile.runctx("fit_psfs(args)",{},
                        {"fit_psfs":fit_psfs,
                         "args":args,},
                        filename="fit_psfs.prof")
    else:
        fit_psfs(args)

    logger.debug('Exiting SHE_CTE_FitPSFs mainMethod()')

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