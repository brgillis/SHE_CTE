""" @file galsim_estimate_shear_interface_test.py

    Created 6 August 2018

    Unit tests for the control shear estimation methods.
"""

__updated__ = "2018-08-06"

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

from os.path import join

from SHE_PPT.file_io import read_pickled_product, find_file
from SHE_PPT.table_formats.shear_estimates import tf as setf
import pytest

from SHE_CTE_ShearEstimation.galsim_estimate_shear import (KSB_estimate_shear, REGAUSS_estimate_shear)
import numpy as np


she_frame_location = "AUX/SHE_PPT/test_she_frame_stack_simple.bin"
ksb_training_location = "AUX/SHE_PPT/test_KSB_training_data.bin"
regauss_training_location = "AUX/SHE_PPT/test_REGAUSS_training_data.bin"

class TestCase:
    """


    """

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.workdir = tmpdir.strpath
        self.logdir = join(tmpdir.strpath, "logs")

        return

    def test_ksb(self):
        """Test that the interface for the KSB method works properly.
        """
        
        # Read in the test data
        she_frame = read_pickled_product(find_file(she_frame_location))
        ksb_training_data= read_pickled_product(find_file(ksb_training_location))
        
        ksb_cat = KSB_estimate_shear(she_frame,
                                     training_data=ksb_training_data,
                                     calibration_data=None,
                                     workdir=self.workdir)
        
        # Check that we have valid data
        for row in ksb_cat:
            for colname in (setf.g1, setf.g2, setf.g1_err, setf.g2_err):
                g = row[colname]
                if not (g>-1 and g<1):
                    raise Exception("Bad value for " + colname + ": " + str(g))
                
        return

    def test_regauss(self):
        """Test that the interface for the REGAUSS method works properly.
        """
        
        # Read in the test data
        she_frame = read_pickled_product(find_file(she_frame_location))
        regauss_training_data= read_pickled_product(find_file(regauss_training_location))
        
        regauss_cat = REGAUSS_estimate_shear(she_frame,
                                             training_data=regauss_training_data,
                                             calibration_data=None,
                                             workdir=self.workdir)
        
        # Check that we have valid data
        for row in regauss_cat:
            for colname in (setf.g1, setf.g2, setf.g1_err, setf.g2_err):
                g = row[colname]
                if not (g>-1 and g<1):
                    raise Exception("Bad value for " + colname + ": " + str(g))
                
        return
