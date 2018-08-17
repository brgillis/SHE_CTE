#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

"""
File: tests/python/meas_stamp_test.py

Created on: 09/08/17
Author: user
"""

import py.test
import numpy as np
import astropy.table

from SHE_MomentsML import meas_stamp
from SHE_MomentsML import utils_table

from SHE_PPT.she_image import SHEImage


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    return lambda x,y: height*np.exp(
        -(((float(center_x)-x)/width_x)**2+((float(center_y)-y)/width_y)**2)/2)


class Test_meas_stamp(object):
   
    
    
    def test_galsim_adamom(self):
        """Make a simple stamp, and measure it"""
        
        
         # We prepare a dummy astropy Table
        cat = astropy.table.Table({"ID":[1]})
        meas_stamp.prep_galsim_adamom(cat, prefix="")
        src = cat[0]
        
        # Prepare a coordinate grid
        X, Y = np.mgrid[0:100, 0:200]+0.5 # So that the arrays contain the coordinates of the pixel *centers*: first pixel goes from 0.0 to 1.0
        assert X[0,0] == 0.5
        assert Y.shape == (100, 200)
        assert Y[0,199] == 199.5
        
        # Draw a Gaussian
        stamp = gaussian(10.0, 50.0, 100.0, 10.0, 20.0)(X, Y) + np.random.randn(*X.shape)
        # We create a SHEImage object out of this numpy array:
        stamp = SHEImage(stamp)
        # And call:
        meas_stamp.galsim_adamom(stamp, src, prefix="")
        
        assert np.hypot(src["x"]-50.0, src["y"]-100.0) < 1.0
        assert src["g1"] < 0.1
            
        #stamp.write_to_fits("/home/user/Work/Projects/SHE_MomentsML/test.fits", clobber=True)
        #print(src)
        #assert(0)
        
    def test_skystats(self): 
        """Similar"""
         
        cat = astropy.table.Table({"ID":[1]})
        meas_stamp.prep_skystats(cat)
        src = cat[0]
        
        stamp = 10000.0 + 5.0 * np.random.randn(50, 50)
        assert stamp.shape == (50, 50)
        stamp = SHEImage(stamp)
        meas_stamp.skystats(stamp, src)
        
        assert src["skystats_smad"] < 7.0
        assert src["skystats_smad"] > 3.0
       
        #print(src)
        #assert False
       
       
        

