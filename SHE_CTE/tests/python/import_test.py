""" @file import_test.py

    Created 20 August 2021

    This module tests importing all modules in each package, to make sure nothing obvious goes wrong, even if the code
    is otherwise untested. 
"""

__updated__ = "2021-08-20"

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

from SHE_CTE import *


class TestImports():

    def test_import(self):
        """ Dummy test - the actual testing is done in the wildcard imports above, which can only be done at the
            module level.
        """
        pass
