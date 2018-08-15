""" @file measure_bias_test.py

    Created 25 June 2018

    Unit tests for measuring shear bias statistics.
"""

__updated__ = "2018-08-10"

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

from SHE_PPT import products
from SHE_PPT.file_io import write_xml_product, read_xml_product,\
    get_allowed_filename, write_listfile
from SHE_PPT.math import LinregressStatistics

from SHE_CTE_BiasMeasurement.measure_bias import measure_bias_from_args
import numpy as np


class Args(object):

    def __init__(self):
        self.workdir = None
        self.logdir = None
        self.profile = False
        self.shear_bias_statistics = None
        self.shear_bias_measurement = None
        self.archive_dir = None
        self.webdav_dir = None
        self.webdav_archive = False


class TestMeasureStatistics:
    """Tests for measuring bias statistics.
    """

    @classmethod
    def setup_class(cls):

        # Set up some mock data for the test

        cls.ex_m1_0 = 0
        cls.ex_c1_0 = 0

        cls.g1_0_bias_statistics = LinregressStatistics()
        cls.g1_0_bias_statistics.w = 80.0
        cls.g1_0_bias_statistics.xm = 0.0
        cls.g1_0_bias_statistics.x2m = 0.00019999999494757503
        cls.g1_0_bias_statistics.ym = 0.0
        cls.g1_0_bias_statistics.xym = 0.0001999999955296518

        cls.ex_m1_1 = 0.02
        cls.ex_c1_1 = -0.05

        cls.g1_1_bias_statistics = LinregressStatistics()
        cls.g1_1_bias_statistics.w = 160.0
        cls.g1_1_bias_statistics.xm = 0.0
        cls.g1_1_bias_statistics.x2m = 0.00019999999494757503
        cls.g1_1_bias_statistics.ym = -0.05000000000000001
        cls.g1_1_bias_statistics.xym = 0.00020399999544024476

        cls.ex_m2_0 = 0
        cls.ex_c2_0 = 0

        cls.g2_0_bias_statistics = LinregressStatistics()
        cls.g2_0_bias_statistics.w = 80.0
        cls.g2_0_bias_statistics.xm = 0.019999999552965164
        cls.g2_0_bias_statistics.x2m = 0.0005999999862979166
        cls.g2_0_bias_statistics.ym = 0.02000000000000001
        cls.g2_0_bias_statistics.xym = 0.0005999999865889553

        cls.ex_m2_1 = 0.1
        cls.ex_c2_1 = 0.0

        cls.g2_1_bias_statistics = LinregressStatistics()
        cls.g2_1_bias_statistics.w = 160.0
        cls.g2_1_bias_statistics.xm = 0.019999999552965164
        cls.g2_1_bias_statistics.x2m = 0.0005999999862979166
        cls.g2_1_bias_statistics.ym = 0.022000000000000013
        cls.g2_1_bias_statistics.xym = 0.0006599999852478506

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

        shear_bias_statistics_prod_0 = products.shear_bias_stats.create_shear_bias_statistics_product()
        shear_bias_statistics_prod_0.set_KSB_statistics(self.g1_0_bias_statistics, self.g2_0_bias_statistics)
        shear_bias_statistics_prod_1 = products.shear_bias_stats.create_shear_bias_statistics_product()
        shear_bias_statistics_prod_1.set_KSB_statistics(self.g1_1_bias_statistics, self.g2_1_bias_statistics)

        shear_bias_statistics_filenames_0 = []
        shear_bias_statistics_filenames_1 = []
        shear_bias_statistics_filenames_01 = []

        for i in range(2):

            filename_0 = get_allowed_filename("bias-stats-0", str(i), extension=".xml", subdir=None)
            filename_1 = get_allowed_filename("bias-stats-1", str(i), extension=".xml", subdir=None)

            write_xml_product(shear_bias_statistics_prod_0, join(args.workdir, filename_0))
            write_xml_product(shear_bias_statistics_prod_1, join(args.workdir, filename_1))

            shear_bias_statistics_filenames_0.append(filename_0)
            shear_bias_statistics_filenames_1.append([filename_1, ])  # We should also be able to read a length-1 tuple
            if i % 2 == 0:
                shear_bias_statistics_filenames_01.append(filename_0)
            else:
                shear_bias_statistics_filenames_01.append(filename_1)

        # Test for bias 0

        write_listfile(join(args.workdir, args.shear_bias_statistics), shear_bias_statistics_filenames_0)

        # Call the function
        measure_bias_from_args(args)

        # Read in and check the results
        shear_bias_measurements_product = read_xml_product(join(args.workdir, args.shear_bias_measurements))

        g1_bias, g2_bias = shear_bias_measurements_product.get_KSB_bias_measurements()

        assert_almost_equal(g1_bias.m, self.ex_m1_0)
        assert_almost_equal(g1_bias.c, self.ex_c1_0)

        assert_almost_equal(g2_bias.m, self.ex_m2_0)
        assert_almost_equal(g2_bias.c, self.ex_c2_0)

        # Test for bias 1

        write_listfile(join(args.workdir, args.shear_bias_statistics), shear_bias_statistics_filenames_1)

        # Call the function
        measure_bias_from_args(args)

        # Read in and check the results
        shear_bias_measurements_product = read_xml_product(join(args.workdir, args.shear_bias_measurements))

        g1_bias, g2_bias = shear_bias_measurements_product.get_KSB_bias_measurements()

        assert_almost_equal(g1_bias.m, self.ex_m1_1)
        assert_almost_equal(g1_bias.c, self.ex_c1_1)

        assert_almost_equal(g2_bias.m, self.ex_m2_1)
        assert_almost_equal(g2_bias.c, self.ex_c2_1)

        # Test for bias 01

        write_listfile(join(args.workdir, args.shear_bias_statistics), shear_bias_statistics_filenames_01)

        # Call the function
        measure_bias_from_args(args)

        # Read in and check the results
        shear_bias_measurements_product = read_xml_product(join(args.workdir, args.shear_bias_measurements))

        g1_bias, g2_bias = shear_bias_measurements_product.get_KSB_bias_measurements()

        assert_almost_equal(g1_bias.m, (self.ex_m1_0 + 2 * self.ex_m1_1) / 3)
        assert_almost_equal(g1_bias.c, (self.ex_c1_0 + 2 * self.ex_c1_1) / 3)

        assert_almost_equal(g2_bias.m, (self.ex_m2_0 + 2 * self.ex_m2_1) / 3)
        assert_almost_equal(g2_bias.c, (self.ex_c2_0 + 2 * self.ex_c2_1) / 3)

        return
