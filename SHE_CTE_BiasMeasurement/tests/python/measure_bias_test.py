""" @file measure_bias_test.py

    Created 25 June 2018

    Unit tests for measuring shear bias statistics.
"""
from numpy.testing import assert_almost_equal
from os.path import join

import pytest

from SHE_CTE_BiasMeasurement.measure_bias import measure_bias_from_args
from SHE_PPT import products
from SHE_PPT.file_io import write_xml_product, read_xml_product,\
    get_allowed_filename, write_listfile
from SHE_PPT.math import LinregressStatistics
import numpy as np


__updated__ = "2018-06-25"

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


products.shear_bias_stats.init()
products.shear_bias_measurements.init()


class Args(object):

    def __init__(self):
        self.workdir = None
        self.logdir = None
        self.profile = False
        self.shear_bias_statistics = None
        self.shear_bias_measurement = None


class TestMeasureStatistics:
    """Tests for measuring bias statistics.
    """

    @classmethod
    def setup_class(cls):

        # Set up some mock data for the test

        cls.ex_m1 = 0
        cls.ex_c1 = 0

        cls.g1_bias_statistics = LinregressStatistics()
        cls.g1_bias_statistics.w = 80.0
        cls.g1_bias_statistics.xm = 0.0
        cls.g1_bias_statistics.x2m = 0.00019999999494757503
        cls.g1_bias_statistics.ym = 0.0
        cls.g1_bias_statistics.xym = 0.0001999999955296518
        g1_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])

        cls.ex_m2 = 0
        cls.ex_c2 = 0

        cls.g2_bias_statistics = LinregressStatistics()
        cls.g2_bias_statistics.w = 80.0
        cls.g2_bias_statistics.xm = 0.019999999552965164
        cls.g2_bias_statistics.x2m = 0.0005999999862979166
        cls.g2_bias_statistics.ym = 0.02000000000000001
        cls.g2_bias_statistics.xym = 0.0005999999865889553

        return

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.workdir = tmpdir.strpath
        self.logdir = join(tmpdir.strpath, "logs")

        return

    def test_measure_bias_from_args(self):
        """Try using the measure_bias_from_args function and check the results.
        """

        # Set up the arguments object
        args = Args()
        args.workdir = self.workdir
        args.logdir = self.logdir
        args.shear_bias_statistics = "test_shear_statistics.json"
        args.shear_bias_measurements = "test_shear_bias_measurements.xml"

        # Set up the files to be read in

        shear_bias_statistics_prod = products.shear_bias_stats.create_shear_bias_statistics_product()
        shear_bias_statistics_prod.set_KSB_statistics(self.g1_bias_statistics, self.g2_bias_statistics)

        shear_bias_statistics_filenames = []
        for i in range(10):
            filename = get_allowed_filename("test_shear_bias_stats", str(i), extension=".xml")
            write_xml_product(shear_bias_statistics_prod, join(args.workdir, filename))
            shear_bias_statistics_filenames.append(filename)
        write_listfile(join(args.workdir, args.shear_bias_statistics), shear_bias_statistics_filenames)

        # Call the function
        measure_bias_from_args(args)

        # Read in and check the results
        shear_bias_measurements_product = read_xml_product(join(args.workdir, args.shear_bias_measurements))

        g1_bias_measurements, g2_bias_measurements = shear_bias_measurements_product.get_KSB_bias_measurements()

        assert_almost_equal(g1_bias.m, self.ex_m1)
        assert_almost_equal(g1_bias.c, self.ex_c1)

        assert_almost_equal(g2_bias.m, self.ex_m2)
        assert_almost_equal(g2_bias.c, self.ex_c2)
