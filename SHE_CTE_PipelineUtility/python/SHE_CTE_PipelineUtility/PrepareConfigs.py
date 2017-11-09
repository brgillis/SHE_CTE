""" @file PrepareConfigs.py

    Created 25 Aug 2017

    Main program for preparing configuration files for parallel runs.
"""

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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

import argparse
from SHE_GST_IceBRGpy.logging import getLogger

from SHE_CTE_ShearEstimation import magic_values as mv

def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_PrepareConfigs defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()
    
    # Extra configuration files
    parser.add_argument('--config_files', nargs='*', default=[],
                        help='Extra configuration files. Each will overwrite any values specified in previous ' +
                             'files, but NOT the one specified with the --config-file option.')

    # Add in each allowed option, with a null default
    for option in allowed_options:
        default = allowed_options[option][0]
        option_type = allowed_options[option][1]
        parser.add_argument("--" + option, default=default, type=option_type)

    # Add allowed fixed params
    for allowed_fixed_param in allowed_fixed_params:
        parser.add_argument("--" + allowed_fixed_param, type=float)

    # Add allowed survey settings, with both level and setting possibilities
    for allowed_survey_setting in allowed_survey_settings:

        generation_level = allowed_survey_setting + "_level"
        parser.add_argument("--" + generation_level, type=str)

        settings = allowed_survey_setting + "_setting"
        parser.add_argument("--" + settings, type=str)
        
    # Where to save a list of the files generated
    parser.add_argument("--output_listfile", default="config_files.json", type=str,
                        help="Filename to which to output list of generated config files.")
        
    # Add option to vary noise seed instead
    parser.add_argument("--vary_noise", action="store_true",
                        help='If set, will vary noise seed instead of model seed for images.')
    
    # Add tag for file labeling
    parser.add_argument("--tag", type=str, default="def",
                        help='Tag to add to generated file names.')
    
    # Where to start the model seed
    parser.add_argument("--model_seed_start", type=int, default=1,
                        help="Where to start the model seed from.")
    
    # Where to start the noise seed
    parser.add_argument("--noise_seed_start", type=int, default=0,
                        help="Where to start the noise seed from. If zero, will follow model seed")

    logger.debug('# Exiting SHE_CTE_PrepareConfigs defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_PrepareConfigs mainMethod()')
    logger.debug('#')
        
    prepare_configs_from_args(args)

    logger.debug('# Exiting SHE_CTE_PrepareConfigs mainMethod()')

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