""" @file measure_statistics_test.py

    Created 22 June 2017

    Unit tests for measuring shear bias statistics.
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

from numpy.testing import assert_almost_equal

import pytest

from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics
from SHE_PPT import table_formats
from SHE_PPT.math import BiasMeasurements, linregress_with_errors
from SHE_PPT.table_formats.details import tf as datf
from SHE_PPT.table_formats.shear_estimates import tf as setf
import numpy as np


class TestMeasureStatistics:
    """Tests for measuring bias statistics.
    """

    @classmethod
    def setup_class(cls):

        # Work done in setup() method since we need to use the tmpdir fixture

        return

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.workdir = tmpdir.strpath

        return

    def test_measure_statistics(self):

        # Start by setting up some mock data for the test
        details = table_formats.details.initialise_details_table()
        shear_estimates = table_formats.shear_estimates.initialise_shear_estimates_table()

        ex_m1 = 0
        ex_c1 = 0

        g1_true = np.array([-0.02, -0.01, 0.00, 0.01, 0.02])
        g1_est = (1 + ex_m1) * g1_true + ex_c1 + np.array([0.25, -0.25, 0, -0.25, 0.25])
        g1_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])

        ex_m2 = 0
        ex_c2 = 0

        g2_true = np.array([0.00, 0.01, 0.02, 0.03, 0.04])
        g2_est = (1 + ex_m2) * g2_true + ex_c2 + np.array([-0.25, 0.25, 0, 0.25, -0.25])
        g2_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])

        for i in range(len(g1_true)):
            details.add_row(vals={datf.ID: i})
            shear_estimates.add_row(vals={setf.ID: i})

        details[datf.g1] = g1_true
        details[datf.g2] = g2_true

        shear_estimates[setf.g1] = g1_est
        shear_estimates[setf.g2] = g2_est

        shear_estimates[setf.g1_err] = g1_err
        shear_estimates[setf.g2_err] = g2_err

        # Now, try using the calculate_shear_bias_statistics function and check the results
        g1_bias_stats, g2_bias_stats = calculate_shear_bias_statistics(shear_estimates, details)

        g1_bias = BiasMeasurements(LinregressResults(g1_bias_stats))
        g2_bias = BiasMeasurements(LinregressResults(g2_bias_stats))

        assert_almost_equal(g1_bias.m, ex_m1)
        assert_almost_equal(g1_bias.c, ex_c1)

        assert_almost_equal(g2_bias.m, ex_m2)
        assert_almost_equal(g2_bias.c, ex_c2)

        return
