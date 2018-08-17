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
File: tests/python/utils_io_test.py

Created on: 09/07/17
Author: Malte Tewes
"""

import os
import py.test
import SHE_MomentsML.utils_io
import numpy as np

class Test_utils_io(object):
    
    
    @classmethod
    def setup_class(cls):
        
        # Filenames for testing, will be deleted by teardown_class
        cls.testfilepath_pkl = "test_MomentsML_utils_io.pkl"
        cls.testfilepath_pklgz = "test_MomentsML_utils_io.pkl.gz"
        

    @classmethod
    def teardown_class(cls):
        
        # Delete all potentially created files:
        for testfilepath in [cls.testfilepath_pkl, cls.testfilepath_pklgz]:
            if os.path.exists(testfilepath):
                os.remove(testfilepath)
 
    
    def test_run(self):
        
        a = np.random.randn(10, 20)
        SHE_MomentsML.utils_io.write_pickle(a, self.testfilepath_pkl)
        SHE_MomentsML.utils_io.write_pickle(a, self.testfilepath_pklgz)
        
        a_pkl = SHE_MomentsML.utils_io.read_pickle(self.testfilepath_pkl)
        a_pklgz = SHE_MomentsML.utils_io.read_pickle(self.testfilepath_pklgz)
        
        assert np.array_equal(a, a_pkl)
        assert np.array_equal(a, a_pklgz)
        
        
        
        
        
        
