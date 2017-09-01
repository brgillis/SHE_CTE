""" @file test_estimate_shear.py

    Created 1 Sep 2017

    Unit tests for the control shear estimation methods.

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

import pytest
import numpy as np

from SHE_PPT import SHEImage
from SHE_PPT.magic_values import scale_label

from SHE_CTE_ShearEstimation.estimate_shear import (get_resampled_image,)

class TestCase:
    """


    """


    def test_get_resampled_image(self):
        
        ss_factor = 5
        
        # Make a mock subsampled image
        ss_data = np.zeros((3*ss_factor,3*ss_factor))
        for i in range(3):
            for j in range(3):
                ss_data[ ss_factor*i:ss_factor*i+ss_factor, ss_factor*j:ss_factor*j+ss_factor ] += i + 3*j
                
        ss_data /= ss_data.sum()
                
        ss_image = SHEImage(ss_data)
        ss_image.header[scale_label] = 1./ss_factor
        
        rb_image = get_resampled_image(ss_image,1.)
        
        rb_data = rb_image.data
        
        rb_data /= rb_data.sum()
        
        
        assert True
