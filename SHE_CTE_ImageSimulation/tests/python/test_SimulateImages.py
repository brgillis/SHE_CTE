""" @file SimulateImages.py

    Created 21 Aug 2017

    Unit tests for SimulateImages executable.

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

import pytest

from SHE_CTE_ImageSimulation.SimulateImages import defineSpecificProgramOptions
from SHE_GST_GalaxyImageGeneration.config.config_default import (allowed_options,
                                                            allowed_fixed_params,
                                                            allowed_survey_settings)


class TestCase:
    """
        @brief Unit tests for SimulateImages executable.

    """


    def testArgs(self):
        """
            @brief Tests that defineSpecificProgramOptions() gives the expected default options.
        """
        
        parser = defineSpecificProgramOptions()
    
        args = parser.parse_args()
        
        arg_lib = vars(args)
        
        # Test the top-level arguments are the defaults
        assert arg_lib["profile"] == False
        assert len(arg_lib["config_files"]) == 0
        
        # Test options are all defaults
        for option in allowed_options:
            assert arg_lib[option] == allowed_options[option][0]

        # Test fixed params are all None
        for allowed_fixed_param in allowed_fixed_params:
            assert arg_lib[allowed_fixed_param] == None
    
        # Add allowed survey settings, with both level and setting possibilities
        for allowed_survey_setting in allowed_survey_settings:
    
            generation_level = allowed_survey_setting + "_level"
            assert arg_lib[generation_level] == None
    
            setting = allowed_survey_setting + "_setting"
            assert arg_lib[setting] == None
        
        return
