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
from os.path import join

import pytest

from SHE_CTE_BiasMeasurement.measure_statistics import measure_statistics_from_args
from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics
from SHE_PPT import products
from SHE_PPT import table_formats
from SHE_PPT.file_io import write_xml_product
from SHE_PPT.math import BiasMeasurements, LinregressResults, linregress_with_errors
from SHE_PPT.table_formats.details import tf as datf
from SHE_PPT.table_formats.shear_estimates import tf as setf
import numpy as np


products.details.init()
products.shear_estimates.init()
products.shear_bias_stats.init()


class Args(object):

    def __init__(self):
        self.workdir = None
        self.details_table = None
        self.shear_estimates = None
        self.shear_bias_statistics = None


class TestMeasureStatistics:
    """Tests for measuring bias statistics.
    """

    @classmethod
    def setup_class(cls):

        # Set up some mock data for the test
        cls.details = table_formats.details.initialise_details_table()
        cls.shear_estimates = table_formats.shear_estimates.initialise_shear_estimates_table()

        cls.ex_m1 = 0
        cls.ex_c1 = 0

        g1_true = np.array([-0.02, -0.01, 0.00, 0.01, 0.02])
        g1_est = (1 + cls.ex_m1) * g1_true + cls.ex_c1 + np.array([0.25, -0.25, 0, -0.25, 0.25])
        g1_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])

        cls.ex_m2 = 0
        cls.ex_c2 = 0

        cls.g2_true = np.array([0.00, 0.01, 0.02, 0.03, 0.04])
        cls.g2_est = (1 + cls.ex_m2) * g2_true + cls.ex_c2 + np.array([-0.25, 0.25, 0, 0.25, -0.25])
        cls.g2_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])

        for i in range(len(g1_true)):
            cls.details.add_row(vals={datf.ID: i})
            cls.shear_estimates.add_row(vals={setf.ID: i})

        cls.details[datf.g1] = g1_true
        cls.details[datf.g2] = g2_true

        cls.shear_estimates[setf.g1] = g1_est
        cls.shear_estimates[setf.g2] = g2_est

        cls.shear_estimates[setf.g1_err] = g1_err
        cls.shear_estimates[setf.g2_err] = g2_err

        return

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.workdir = tmpdir.strpath

        return

    def test_calculate_shear_bias_statistics(self):
        """Try using the calculate_shear_bias_statistics function and check the results.
        """
        g1_bias_stats, g2_bias_stats = calculate_shear_bias_statistics(self.shear_estimates, self.details)

        g1_bias = BiasMeasurements(LinregressResults(g1_bias_stats))
        g2_bias = BiasMeasurements(LinregressResults(g2_bias_stats))

        assert_almost_equal(g1_bias.m, self.ex_m1)
        assert_almost_equal(g1_bias.c, self.ex_c1)

        assert_almost_equal(g2_bias.m, self.ex_m2)
        assert_almost_equal(g2_bias.c, self.ex_c2)

        return

    def test_measure_statistics_from_args(self):
        """Try using the measure_statistics_from_args function and check the results.
        """

        # Set up the arguments object
        args = Args()
        args.workdir = self.workdir
        args.details_table = "test_details_table.xml"
        args.shear_estimates = "test_shear_estimates.xml"
        args.shear_bias_statistics = "test_shear_statistics.xml"

        # Set up the files to be read in

        details_filename = "test_details_table.fits"
        details_product = products.details.create_details_product(details_filename)
        write_xml_product(details_product, join(args.workdir, args.details_table))
        self.details.write(join(args.workdir, args.details_filename), format="fits")

        shear_estimates_filename = "test_shear_estimates.fits"
        shear_estimates_product = products.shear_estimates.create_shear_estimates_product(
            KSB_filename=shear_estimates_filename)
        write_xml_product(shear_estimates_product, join(args.workdir, args.shear_estimates))
        self.shear_estimates.write(join(args.workdir, args.shear_estimates_filename), format="fits")

        # Call the function
        measure_statistics_from_args(args)

        # Read in and check the results
        shear_bias_statistics_product = read_xml_product(join(args.workdir, args.shear_bias_statistics))

        g1_bias_stats, g2_bias_stats = shear_bias_statistics_product.get_KSB_statistics()

        g1_bias = BiasMeasurements(LinregressResults(g1_bias_stats))
        g2_bias = BiasMeasurements(LinregressResults(g2_bias_stats))

        assert_almost_equal(g1_bias.m, self.ex_m1)
        assert_almost_equal(g1_bias.c, self.ex_c1)

        assert_almost_equal(g2_bias.m, self.ex_m2)
        assert_almost_equal(g2_bias.c, self.ex_c2)
