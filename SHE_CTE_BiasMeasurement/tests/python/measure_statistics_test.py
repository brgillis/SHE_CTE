""" @file measure_statistics_test.py

    Created 22 June 2018

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

from copy import deepcopy
from numpy.testing import assert_almost_equal
from os.path import join
import pytest

from SHE_PPT import products
from SHE_PPT import table_formats
from SHE_PPT.file_io import write_xml_product, read_xml_product
from SHE_PPT.math import BiasMeasurements, LinregressResults, linregress_with_errors
from SHE_PPT.table_formats.details import tf as datf
from SHE_PPT.table_formats.shear_estimates import tf as setf

from SHE_CTE_BiasMeasurement.measure_statistics import measure_statistics_from_args
from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics
import numpy as np


class Args(object):

    def __init__(self):
        self.workdir = None
        self.logdir = None
        self.profile = False
        self.details_table = None
        self.shear_estimates = None
        self.shear_bias_statistics = None
        self.archive_dir = None
        self.webdav_dir = None
        self.webdav_archive = False


class TestMeasureStatistics:
    """Tests for measuring bias statistics.
    """

    @classmethod
    def setup_class(cls):

        # Set up some mock data for the test
        cls.details = table_formats.details.initialise_details_table()
        cls.details_group = table_formats.details.initialise_details_table()
        cls.shear_estimates = table_formats.shear_estimates.initialise_shear_estimates_table()
        cls.shear_estimates_group = table_formats.shear_estimates.initialise_shear_estimates_table()

        cls.len_group = 2  # Number of galaxies per group

        # G1 data

        cls.ex_m1 = 0.02
        cls.ex_c1 = -0.05

        g1_true = np.array([-0.02, -0.01, 0.00, 0.01, 0.02])
        g1_est = (1 + cls.ex_m1) * g1_true + cls.ex_c1 + np.array([0.25, -0.25, 0, -0.25, 0.25])
        g1_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])
        cls.ex_w1 = np.sum(g1_err**-2)

        g1_true_group = np.concatenate((g1_true, g1_true))
        g1_est_group = np.concatenate((g1_est, g1_est))
        g1_err_group = np.concatenate((g1_err, g1_err))

        g1_est_group += np.array([0.01] * 5 + [-0.01] * 5)

        # G2 data

        cls.ex_m2 = 0.1
        cls.ex_c2 = 0.0

        g2_true = np.array([0.00, 0.01, 0.02, 0.03, 0.04])
        g2_est = (1 + cls.ex_m2) * g2_true + cls.ex_c2 + np.array([-0.25, 0.25, 0, 0.25, -0.25])
        g2_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])
        cls.ex_w2 = np.sum(g2_err**-2)

        g2_true_group = np.concatenate((g2_true, g2_true))
        g2_est_group = np.concatenate((g2_est, g2_est))
        g2_err_group = np.concatenate((g2_err, g2_err))

        g2_est_group += np.array([-0.2] * 5 + [0.2] * 5)

        for i in range(len(g1_true)):
            cls.details.add_row(vals={datf.ID: i,
                                      datf.group_ID: i})
            cls.shear_estimates.add_row(vals={setf.ID: i})
        for j in range(cls.len_group):
            for i in range(len(g1_true)):
                cls.details_group.add_row(vals={datf.ID: i * cls.len_group + j,
                                                datf.group_ID: i})
                cls.shear_estimates_group.add_row(vals={setf.ID: i * cls.len_group + j})

        # Save details data
        cls.details[datf.g1] = g1_true
        cls.details[datf.g2] = g2_true

        # Save details_group data
        cls.details_group[datf.g1] = g1_true_group
        cls.details_group[datf.g2] = g2_true_group

        # Save shear_estimates data
        cls.shear_estimates[setf.g1] = g1_est
        cls.shear_estimates[setf.g2] = g2_est
        cls.shear_estimates[setf.g1_err] = g1_err
        cls.shear_estimates[setf.g2_err] = g2_err

        # Save shear_estimates_group data
        cls.shear_estimates_group[setf.g1] = g1_est_group
        cls.shear_estimates_group[setf.g2] = g2_est_group
        cls.shear_estimates_group[setf.g1_err] = g1_err_group
        cls.shear_estimates_group[setf.g2_err] = g2_err_group

        return

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.workdir = tmpdir.strpath
        self.logdir = join(tmpdir.strpath, "logs")

        return

    def test_calculate_shear_bias_statistics(self):
        """Try using the calculate_shear_bias_statistics function and check the results.
        """
        g1_bias_stats, g2_bias_stats = calculate_shear_bias_statistics(self.shear_estimates, self.details)

        g1_bias = BiasMeasurements(LinregressResults(g1_bias_stats))
        g2_bias = BiasMeasurements(LinregressResults(g2_bias_stats))

        assert_almost_equal(g1_bias_stats.w, self.ex_w1)
        assert_almost_equal(g2_bias_stats.w, self.ex_w2)

        assert_almost_equal(g1_bias.m, self.ex_m1)
        assert_almost_equal(g1_bias.c, self.ex_c1)

        assert_almost_equal(g2_bias.m, self.ex_m2)
        assert_almost_equal(g2_bias.c, self.ex_c2)

        return

    def test_bad_shear_bias_input(self):
        """Tests that the calculate_shear_bias_statistics method is resilient to bad measurements.
        """

        # Set up data with some bad input

        somebad_shear_estimates = table_formats.shear_estimates.initialise_shear_estimates_table()
        allbad_shear_estimates = table_formats.shear_estimates.initialise_shear_estimates_table()

        for i in range(len(self.shear_estimates)):
            somebad_shear_estimates.add_row(vals={setf.ID: i})
            allbad_shear_estimates.add_row(vals={setf.ID: i})

        somebad_shear_estimates[setf.g1] = self.shear_estimates[setf.g1]
        somebad_shear_estimates[setf.g2] = self.shear_estimates[setf.g2]
        somebad_shear_estimates[setf.g1_err] = self.shear_estimates[setf.g1_err]
        somebad_shear_estimates[setf.g2_err] = self.shear_estimates[setf.g2_err]

        allbad_shear_estimates[setf.g1] = self.shear_estimates[setf.g1]
        allbad_shear_estimates[setf.g2] = self.shear_estimates[setf.g2]
        allbad_shear_estimates[setf.g1_err] = self.shear_estimates[setf.g1_err]
        allbad_shear_estimates[setf.g2_err] = self.shear_estimates[setf.g2_err]

        for i in [0, 3]:
            row = somebad_shear_estimates[i]
            row[setf.g1], row[setf.g2] = np.nan, np.nan
            row[setf.g1_err], row[setf.g2_err] = np.inf, np.inf
        for i in range(len(self.shear_estimates)):
            row = allbad_shear_estimates[i]
            row[setf.g1], row[setf.g2] = np.nan, np.nan
            row[setf.g1_err], row[setf.g2_err] = np.inf, np.inf

        # Check that we don't crash by calling calculate_shear_bias_statistics
        g1_somebad_bias_stats, g2_somebad_bias_stats = calculate_shear_bias_statistics(
            somebad_shear_estimates, self.details)
        g1_allbad_bias_stats, g2_allbad_bias_stats = calculate_shear_bias_statistics(
            allbad_shear_estimates, self.details)

        # Check the statistics all have the expected weight

        assert_almost_equal(g1_somebad_bias_stats.w, 0.6 * self.ex_w1)
        assert_almost_equal(g2_somebad_bias_stats.w, 0.6 * self.ex_w2)

        assert not np.isnan(g1_somebad_bias_stats.xm)
        assert not np.isnan(g1_somebad_bias_stats.x2m)
        assert not np.isnan(g1_somebad_bias_stats.ym)
        assert not np.isnan(g1_somebad_bias_stats.xym)

        assert_almost_equal(g1_allbad_bias_stats.w, 0.)
        assert_almost_equal(g2_allbad_bias_stats.w, 0.)

        return

    def test_bad_shear_bias_input_group(self):
        """Tests that the calculate_shear_bias_statistics method is resilient to bad measurements.
        """

        # Set up data with some bad input

        somebad_shear_estimates = table_formats.shear_estimates.initialise_shear_estimates_table()
        allbad_shear_estimates = table_formats.shear_estimates.initialise_shear_estimates_table()

        for j in range(self.len_group):
            for i in range(len(self.shear_estimates)):
                somebad_shear_estimates.add_row(vals={setf.ID: i * self.len_group + j})
                allbad_shear_estimates.add_row(vals={setf.ID: i * self.len_group + j})

        somebad_shear_estimates[setf.g1] = self.shear_estimates_group[setf.g1]
        somebad_shear_estimates[setf.g2] = self.shear_estimates_group[setf.g2]
        somebad_shear_estimates[setf.g1_err] = self.shear_estimates_group[setf.g1_err]
        somebad_shear_estimates[setf.g2_err] = self.shear_estimates_group[setf.g2_err]

        allbad_shear_estimates[setf.g1] = self.shear_estimates_group[setf.g1]
        allbad_shear_estimates[setf.g2] = self.shear_estimates_group[setf.g2]
        allbad_shear_estimates[setf.g1_err] = self.shear_estimates_group[setf.g1_err]
        allbad_shear_estimates[setf.g2_err] = self.shear_estimates_group[setf.g2_err]

        for i in [0, 3, 7, 8]:
            row = somebad_shear_estimates[i]
            row[setf.g1], row[setf.g2] = np.nan, np.nan
            row[setf.g1_err], row[setf.g2_err] = np.inf, np.inf
        for i in range(len(self.shear_estimates_group)):
            row = allbad_shear_estimates[i]
            row[setf.g1], row[setf.g2] = np.nan, np.nan
            row[setf.g1_err], row[setf.g2_err] = np.inf, np.inf

        # Check that we don't crash by calling calculate_shear_bias_statistics
        g1_somebad_bias_stats, g2_somebad_bias_stats = calculate_shear_bias_statistics(
            somebad_shear_estimates, self.details_group)
        g1_allbad_bias_stats, g2_allbad_bias_stats = calculate_shear_bias_statistics(
            allbad_shear_estimates, self.details_group)

        # Check the statistics all have the expected weight

        assert_almost_equal(g1_somebad_bias_stats.w, self.len_group * 0.6 * self.ex_w1)
        assert_almost_equal(g2_somebad_bias_stats.w, self.len_group * 0.6 * self.ex_w2)

        assert not np.isnan(g1_somebad_bias_stats.xm)
        assert not np.isnan(g1_somebad_bias_stats.x2m)
        assert not np.isnan(g1_somebad_bias_stats.ym)
        assert not np.isnan(g1_somebad_bias_stats.xym)

        assert_almost_equal(g1_allbad_bias_stats.w, 0.)
        assert_almost_equal(g2_allbad_bias_stats.w, 0.)

        return

    def test_calculate_shear_bias_statistics_group(self):
        """Try using the calculate_shear_bias_statistics function and check the results on grouped data.
        """
        g1_bias_stats, g2_bias_stats = calculate_shear_bias_statistics(self.shear_estimates_group, self.details_group)

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
        args.logdir = self.logdir
        args.details_table = "test_details_table.xml"
        args.shear_estimates = "test_shear_estimates.xml"
        args.shear_bias_statistics = "test_shear_statistics.xml"

        # Set up the files to be read in

        details_filename = "test_details_table.fits"
        details_product = products.details.create_details_product(details_filename)
        write_xml_product(details_product, args.details_table, workdir=args.workdir)
        self.details.write(join(args.workdir, details_filename), format="fits")

        shear_estimates_filename = "test_shear_estimates.fits"
        shear_estimates_product = products.shear_estimates.create_shear_estimates_product(
            KSB_filename=shear_estimates_filename)
        write_xml_product(shear_estimates_product, args.shear_estimates, workdir=args.workdir)
        self.shear_estimates.write(join(args.workdir, shear_estimates_filename), format="fits")

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
