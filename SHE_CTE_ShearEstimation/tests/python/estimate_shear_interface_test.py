""" @file estimate_shear_interface_test.py

    Created 6 August 2018

    Unit tests for the control shear estimation methods.
"""

__updated__ = "2018-11-16"

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

from copy import deepcopy
from os.path import join
import pytest
import time

from SHE_CTE_ShearEstimation.bfd_measure_moments import bfd_measure_moments
from SHE_CTE_ShearEstimation.galsim_estimate_shear import (KSB_estimate_shear, REGAUSS_estimate_shear)
import SHE_LensMC.SHE_measure_shear
from SHE_MomentsML.estimate_shear import estimate_shear as ML_estimate_shear
from SHE_PPT.file_io import read_pickled_product, find_file
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.shear_estimates import tf as setf
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
        original_she_frame = deepcopy(she_frame)
        ksb_training_data = read_pickled_product(find_file(ksb_training_location))

        ksb_cat = KSB_estimate_shear(she_frame,
                                     training_data=ksb_training_data,
                                     calibration_data=None,
                                     workdir=self.workdir)

        # Check that we have valid data
        for row in ksb_cat:
            for colname in (setf.g1, setf.g2, setf.g1_err, setf.g2_err):
                g = row[colname]
                if not (g > -1 and g < 1):
                    raise Exception("Bad value for " + colname + ": " + str(g))

        # Check that the input data isn't changed by the method
        assert she_frame == original_she_frame

        return

    def test_regauss(self):
        """Test that the interface for the REGAUSS method works properly.
        """

        # Read in the test data
        she_frame = read_pickled_product(find_file(she_frame_location))
        original_she_frame = deepcopy(she_frame)
        regauss_training_data = read_pickled_product(find_file(regauss_training_location))

        regauss_cat = REGAUSS_estimate_shear(she_frame,
                                             training_data=regauss_training_data,
                                             calibration_data=None,
                                             workdir=self.workdir)

        # Check that we have valid data
        for row in regauss_cat:
            for colname in (setf.g1, setf.g2, setf.g1_err, setf.g2_err):
                g = row[colname]
                if not (g > -1 and g < 1):
                    raise Exception("Bad value for " + colname + ": " + str(g))

        # Check that the input data isn't changed by the method
        assert she_frame == original_she_frame

        return

    @pytest.mark.skip(reason="Duplicated from SHE_LensMC")
    def test_lensmc(self):
        """Test that the interface for the LensMC method works properly.
        """

        # define logger object
        logger = getLogger(__name__)

        logger.debug('Load auxiliary data.')

        # access mock data in PPT
        she_frame = read_pickled_product(find_file(she_frame_location))
        training_data = read_pickled_product(find_file('AUX/SHE_LensMC/test_LensMC_training_data.bin'))

        original_she_frame = deepcopy(she_frame)

        # make sure we seed the random number generator to produce consistent results
        seed = 67794

        logger.debug('Call fit_frame_stack().')

        t0 = time.time()

        # call the SHE wrapper of LensMC
        shear_estimates = SHE_LensMC.SHE_measure_shear.fit_frame_stack(
            she_frame, training_data=training_data, seed=seed)

        logger.info('Runtime = {:.3f}'.format(time.time() - t0))

        for row in shear_estimates:
            logger.info('ID = {}, chi2 = {:.3f}'.format(row[setf.ID], row[setf.chi2]))

        logger.debug('Check results.')

        # check results
        for row in shear_estimates:
            for colname in (setf.g1, setf.g2):
                g = row[colname]
                assert (g > -1 and g < 1)
            assert row[setf.g1_err] > 0
            assert row[setf.g2_err] > 0
            assert row[setf.re] > 0
            assert row[setf.flux] > 0
            assert row[setf.snr] > 0
            assert row[setf.chi2] > 0
            assert row[setf.dof] > 0

        # Check that the input data isn't changed by the method
        assert she_frame == original_she_frame

        return

    @pytest.mark.skip(reason="Errors are known issue to be addressed")
    def test_momentsml(self):
        """Test that the interface for the MomentsML method works properly.
        """

        # Read in the test data
        she_frame = read_pickled_product(find_file(she_frame_location))
        original_she_frame = deepcopy(she_frame)

        momentsml_cat = ML_estimate_shear(she_frame,
                                          training_data=None,
                                          calibration_data=None,
                                          workdir=self.workdir)

        # Check that we have valid data
        for row in momentsml_cat:
            for colname in (setf.g1, setf.g2, setf.g1_err, setf.g2_err):
                g = row[colname]
                if not (g > -1 and g < 1):
                    raise Exception("Bad value for " + colname + ": " + str(g))

        # Check that the input data isn't changed by the method
        assert she_frame == original_she_frame

        return

    def test_bfd(self):
        """Test that the interface for the BFD method works properly.
        """

        # Read in the test data
        she_frame = read_pickled_product(find_file(she_frame_location))
        original_she_frame = deepcopy(she_frame)

        bfd_cat = bfd_measure_moments(she_frame,
                                      training_data=None,
                                      calibration_data=None,
                                      workdir=self.workdir)

        # TODO: Add test that values in the catalogue are reasonable

        # Check that the input data isn't changed by the method
        assert she_frame == original_she_frame

        return
