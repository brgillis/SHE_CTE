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
File: tests/python/utils_table_test.py

Created on: 09/07/17
Author: Malte Tewes
"""

import astropy.table

import py.test
import SHE_MomentsML.utils_table

class Test_utils_table(object):
    
    def test_get_info(self):
        
        cat = astropy.table.Table({"a":[1,2,3], "b":[5,6,7]})
        
        txt = SHE_MomentsML.utils_table.get_info(cat)
        
        assert len(txt) > 10
        
        
        
        
        
        
        
