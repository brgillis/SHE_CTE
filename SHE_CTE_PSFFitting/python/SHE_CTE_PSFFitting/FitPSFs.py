""" @file SimulateImages.py

    Created 21 Aug 2017

    Executable for generating simulated images.

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

from SHE_CTE_PSFFitting import magic_values as mv
from SHE_GST_IceBRGpy.logging import getLogger

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