""" @file cti_residual_validation_test.py

    Created 15 Feb 2018

    Unit tests for CTI residual validation.
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

import numpy as np
import pytest

from SHE_CTE_ShearValidation.cti_residual_validation import (validate_cti_ellipticity_residual_bin,
                                                             validate_cti_ellipticity_residuals,
                                                             readout_split,
                                                             lower_fail_flag_offset,
                                                             upper_fail_flag_offset)


class TestCTIResidualValidation:
    """
        Unit tests for CTI residual validation.
    """
    
    @classmethod
    def setup_class(cls):
        
        # Seed for consistent results
        np.random.seed(1234)
        
        cls.nvals = 1000
        
        # Positions
        cls.r_vals = np.random.random(cls.nvals)*8272.
        
        # Errors on shear measurements
        cls.g1_err = 0.2 + np.random.random(cls.nvals)*0.1
        cls.g2_err = 0.2 + np.random.random(cls.nvals)*0.1
        
        # Mock measurements with no correlation
        cls.good_g1_vals = np.random.randn(cls.nvals)*cls.g1_err
        cls.good_g2_vals = np.random.randn(cls.nvals)*cls.g2_err
        
        low_vals = cls.r_vals < readout_split
        
        # Mock measurements with a correlation with readout distance in both halves
        cls.badbad_g1_vals = np.where(low_vals,
                                      cls.good_g1_vals + 0.5/4000.*cls.r_vals,
                                      cls.good_g1_vals + 0.5/4000.*(8272-cls.r_vals),)
        
        cls.badbad_g2_vals = np.where(low_vals,
                                      cls.good_g2_vals + -1.0/4000.*cls.r_vals,
                                      cls.good_g2_vals + -1.0/4000.*(8272-cls.r_vals),)
        
        # Mock measurements with a correlation with readout distance in bottom half only
        cls.badgood_g1_vals = np.where(low_vals,
                                       cls.good_g1_vals + 0.5/4000.*cls.r_vals,
                                       cls.good_g1_vals,)
        
        cls.badgood_g2_vals = np.where(low_vals,
                                       cls.good_g2_vals + -1.0/4000.*cls.r_vals,
                                       cls.good_g2_vals,)
        
        # Mock measurements with a correlation with readout distance in top half only
        cls.goodbad_g1_vals = np.where(low_vals,
                                       cls.good_g1_vals,
                                       cls.good_g1_vals + 0.5/4000.*(8272-cls.r_vals),)
        
        cls.goodbad_g2_vals = np.where(low_vals,
                                       cls.good_g2_vals,
                                       cls.good_g2_vals + -1.0/4000.*(8272-cls.r_vals),)
        
        return
        
    @classmethod
    def teardown_class(cls):
        
        del (cls.r_vals, cls.g1_err, cls.g2_err, cls.good_g1_vals, cls.good_g2_vals, cls.badbad_g1_vals, cls.badbad_g2_vals,
             cls.badgood_g1_vals, cls.badgood_g2_vals, cls.goodbad_g1_vals, cls.goodbad_g2_vals,)
        
        return
    
    def test_validate_cti_ellipticity_residual_bin(self):
        
        # Check it gives flag of 0 for good data
        assert validate_cti_ellipticity_residual_bin(self.r_vals, self.good_g1_vals) == 0
        assert validate_cti_ellipticity_residual_bin(self.r_vals, self.good_g2_vals) == 0
        
        # Check it gives proper flag for fully bad data 
        assert (validate_cti_ellipticity_residual_bin(self.r_vals, self.badbad_g1_vals) ==
                lower_fail_flag_offset + upper_fail_flag_offset)
        assert (validate_cti_ellipticity_residual_bin(self.r_vals, self.badbad_g2_vals) ==
                lower_fail_flag_offset + upper_fail_flag_offset)
        
        # Check it gives proper flag for bad data in the bottom half only
        assert (validate_cti_ellipticity_residual_bin(self.r_vals, self.badgood_g1_vals) ==
                lower_fail_flag_offset)
        assert (validate_cti_ellipticity_residual_bin(self.r_vals, self.badgood_g2_vals) ==
                lower_fail_flag_offset)
        
        # Check it gives proper flag for bad data in the top half only
        assert (validate_cti_ellipticity_residual_bin(self.r_vals, self.goodbad_g1_vals) ==
                upper_fail_flag_offset)
        assert (validate_cti_ellipticity_residual_bin(self.r_vals, self.goodbad_g2_vals) ==
                upper_fail_flag_offset)
        
        return
