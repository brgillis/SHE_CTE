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
File: tests/python/params_test.py

Created on: 09/07/17
Author: Malte Tewes
"""

import py.test
import SHE_MomentsML.params
import os

class Test_params(object):
    
    def setup_class(cls):    
        # Filenames for testing, will be deleted by teardown_class
        cls.testfilepath = "test_MomentsML_params.pkl"
        
    @classmethod
    def teardown_class(cls):
        # Delete all potentially created files:
        for testfilepath in [cls.testfilepath]:
            if os.path.exists(testfilepath):
                os.remove(testfilepath)
 

    def test(self):
        params = SHE_MomentsML.params.MomentsMLParams()
        params.write_pickle(self.testfilepath)
        
        p2 = SHE_MomentsML.params.MomentsMLParams.read_pickle(self.testfilepath)
        
        
        
        
        
        
        
