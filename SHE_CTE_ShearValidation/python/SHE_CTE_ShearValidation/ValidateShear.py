""" @file ValidateShear.py

    Created 12 Oct 2017

    Executable for validating shear measurements and creating a combined catalogue.

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

from SHE_CTE_ValidateShear import magic_values as mv
from SHE_CTE.magic_values import force_dry_run

if force_dry_run:
    from SHE_CTE_ValidateShear.validate_shear_dry import validate_shear
else:
    from SHE_CTE_ValidateShear.validate_shear import validate_shear

def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program, using all possible configurations.

    @return
        An  ArgumentParser.
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ValidateShear defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Option for profiling
    parser.add_argument('--profile',action='store_true',
                        help='Store profiling data for execution.')
    
    # Input filenames
    parser.add_argument('shear_estimates_product',type=str,
                        help='Filename for shear estimates data product (XML data product)')
    parser.add_argument('shear_estimates_listfile',type=str,
                        help='Filename for listfile associated with shear estimates data product (.json listfile)')
    parser.add_argument('shear_validation_statistics_table',type=str,
                        help='Filename for table of shear validation statistics.')
    
    # Output filenames
    parser.add_argument('validated_shear_estimates_table',type=str,
                        help='Desired filename for output shear estimates table (multi-HDU fits file).')
    

    logger.debug('Exiting SHE_CTE_ValidateShear defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_ValidateShear mainMethod()')
    logger.debug('#')
        
    if args.profile:
        import cProfile
        cProfile.runctx("validate_shear(args)",{},
                        {"validate_shear":validate_shear,
                         "args":args,},
                        filename="validate_shear.prof")
    else:
        validate_shear(args)

    logger.debug('Exiting SHE_CTE_ValidateShear mainMethod()')

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