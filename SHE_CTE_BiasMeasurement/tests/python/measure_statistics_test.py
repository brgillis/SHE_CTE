""" @file measure_statistics_test.py

    Created 22 June 2018

    Unit tests for measuring shear bias statistics.
"""

__updated__ = "2021-08-18"

# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details_table.
#
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

import os
from os.path import join

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from SHE_CTE_BiasMeasurement.MeasureStatistics import D_MS_CONFIG_DEFAULTS, D_MS_CONFIG_TYPES
from SHE_CTE_BiasMeasurement.measure_statistics import measure_statistics_from_args
from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics
from SHE_PPT import products, table_formats
from SHE_PPT.file_io import read_xml_product, write_xml_product
from SHE_PPT.math import BiasMeasurements, LinregressResults
from SHE_PPT.pipeline_utility import read_calibration_config
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksbm_tf
from SHE_PPT.table_formats.she_simulated_catalog import tf as datf


class Args(object):

    def __init__(self):
        self.workdir = None
        self.logdir = None
        self.profile = False
        self.details_table = None
        self.shear_estimates = None
        self.she_bias_statistics = None
        self.archive_dir = None
        self.webdav_dir = None
        self.webdav_archive = False
        self.number_threads = 1
        self.pipeline_config = read_calibration_config(None,
                                                       d_defaults = D_MS_CONFIG_DEFAULTS,
                                                       d_types = D_MS_CONFIG_TYPES)
        self.store_measurements_only = False


class TestMeasureStatistics:
    """Tests for measuring bias statistics.
    """

    @classmethod
    def setup_class(cls):

        # Set up some mock data for the test
        cls.details_table = table_formats.she_simulated_catalog.tf.init_table()
        cls.details_table_group = table_formats.she_simulated_catalog.tf.init_table()
        cls.she_measurements = table_formats.she_ksb_measurements.tf.init_table()
        cls.she_measurements_group = table_formats.she_ksb_measurements.tf.init_table()

        cls.len_group = 2  # Number of galaxies per group

        # G1 data

        cls.ex_m1 = 0.02
        cls.ex_c1 = -0.05

        g1_true = np.array([-0.02, -0.01, 0.00, 0.01, 0.02])
        g1_est = (1 + cls.ex_m1) * g1_true + cls.ex_c1 + np.array([0.25, -0.25, 0, -0.25, 0.25])
        g1_err = np.array([0.25, 0.25, 0.25, 0.25, 0.25])
        cls.ex_w1 = np.sum(g1_err ** -2)

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
        cls.ex_w2 = np.sum(g2_err ** -2)

        g2_true_group = np.concatenate((g2_true, g2_true))
        g2_est_group = np.concatenate((g2_est, g2_est))
        g2_err_group = np.concatenate((g2_err, g2_err))

        g2_est_group += np.array([-0.2] * 5 + [0.2] * 5)

        for i in range(len(g1_true)):
            cls.details_table.add_row(vals = {datf.ID      : i,
                                              datf.group_ID: i})
            cls.she_measurements.add_row(vals = {ksbm_tf.ID: i})
        for j in range(cls.len_group):
            for i in range(len(g1_true)):
                cls.details_table_group.add_row(vals = {datf.ID      : i * cls.len_group + j,
                                                        datf.group_ID: i})
                cls.she_measurements_group.add_row(vals = {ksbm_tf.ID: i * cls.len_group + j})

        # Save details_table data
        cls.details_table[datf.g1] = g1_true
        cls.details_table[datf.g2] = g2_true

        # Save details_table_group data
        cls.details_table_group[datf.g1] = g1_true_group
        cls.details_table_group[datf.g2] = g2_true_group

        # Save she_measurements data
        cls.she_measurements[ksbm_tf.g1] = g1_est
        cls.she_measurements[ksbm_tf.g2] = g2_est
        cls.she_measurements[ksbm_tf.g1_err] = g1_err
        cls.she_measurements[ksbm_tf.g2_err] = g2_err

        # Save she_measurements_group data
        cls.she_measurements_group[ksbm_tf.g1] = g1_est_group
        cls.she_measurements_group[ksbm_tf.g2] = g2_est_group
        cls.she_measurements_group[ksbm_tf.g1_err] = g1_err_group
        cls.she_measurements_group[ksbm_tf.g2_err] = g2_err_group

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse = True)
    def setup(self, tmpdir):
        self.workdir = tmpdir.strpath
        self.logdir = join(tmpdir.strpath, "logs")

    def test_calculate_shear_bias_statistics(self):
        """Try using the calculate_shear_bias_statistics function and check the results.
        """
        g1_bias_stats, g2_bias_stats = calculate_shear_bias_statistics(self.she_measurements, self.details_table)

        g1_bias = BiasMeasurements(LinregressResults(g1_bias_stats))
        g2_bias = BiasMeasurements(LinregressResults(g2_bias_stats))

        assert_almost_equal(g1_bias_stats.w, self.ex_w1, decimal = 4)
        assert_almost_equal(g2_bias_stats.w, self.ex_w2, decimal = 4)

        assert_almost_equal(g1_bias.m, self.ex_m1, decimal = 4)
        assert_almost_equal(g1_bias.c, self.ex_c1, decimal = 4)

        assert_almost_equal(g2_bias.m, self.ex_m2, decimal = 4)
        assert_almost_equal(g2_bias.c, self.ex_c2, decimal = 4)

    def test_bad_shear_bias_input(self):
        """Tests that the calculate_shear_bias_statistics method is resilient to bad measurements.
        """

        # Set up data with some bad input

        somebad_she_measurements = table_formats.she_ksb_measurements.tf.init_table()
        allbad_she_measurements = table_formats.she_ksb_measurements.tf.init_table()

        for i in range(len(self.she_measurements)):
            somebad_she_measurements.add_row(vals = {ksbm_tf.ID: i})
            allbad_she_measurements.add_row(vals = {ksbm_tf.ID: i})

        somebad_she_measurements[ksbm_tf.g1] = self.she_measurements[ksbm_tf.g1]
        somebad_she_measurements[ksbm_tf.g2] = self.she_measurements[ksbm_tf.g2]
        somebad_she_measurements[ksbm_tf.g1_err] = self.she_measurements[ksbm_tf.g1_err]
        somebad_she_measurements[ksbm_tf.g2_err] = self.she_measurements[ksbm_tf.g2_err]

        allbad_she_measurements[ksbm_tf.g1] = self.she_measurements[ksbm_tf.g1]
        allbad_she_measurements[ksbm_tf.g2] = self.she_measurements[ksbm_tf.g2]
        allbad_she_measurements[ksbm_tf.g1_err] = self.she_measurements[ksbm_tf.g1_err]
        allbad_she_measurements[ksbm_tf.g2_err] = self.she_measurements[ksbm_tf.g2_err]

        for i in [0, 3]:
            row = somebad_she_measurements[i]
            row[ksbm_tf.g1], row[ksbm_tf.g2] = np.nan, np.nan
            row[ksbm_tf.g1_err], row[ksbm_tf.g2_err] = np.inf, np.inf
        for i in range(len(self.she_measurements)):
            row = allbad_she_measurements[i]
            row[ksbm_tf.g1], row[ksbm_tf.g2] = np.nan, np.nan
            row[ksbm_tf.g1_err], row[ksbm_tf.g2_err] = np.inf, np.inf

        # Check that we don't crash by calling calculate_shear_bias_statistics
        g1_somebad_bias_stats, g2_somebad_bias_stats = calculate_shear_bias_statistics(
            somebad_she_measurements, self.details_table)
        g1_allbad_bias_stats, g2_allbad_bias_stats = calculate_shear_bias_statistics(
            allbad_she_measurements, self.details_table)

        # Check the statistics all have the expected weight

        assert_almost_equal(g1_somebad_bias_stats.w, 0.6 * self.ex_w1, decimal = 4)
        assert_almost_equal(g2_somebad_bias_stats.w, 0.6 * self.ex_w2, decimal = 4)

        assert not np.isnan(g1_somebad_bias_stats.xm)
        assert not np.isnan(g1_somebad_bias_stats.x2m)
        assert not np.isnan(g1_somebad_bias_stats.ym)
        assert not np.isnan(g1_somebad_bias_stats.xym)

        assert_almost_equal(g1_allbad_bias_stats.w, 0., decimal = 4)
        assert_almost_equal(g2_allbad_bias_stats.w, 0., decimal = 4)

    def test_bad_shear_bias_input_group(self):
        """Tests that the calculate_shear_bias_statistics method is resilient to bad measurements.
        """

        # Set up data with some bad input

        somebad_she_measurements = table_formats.she_ksb_measurements.tf.init_table()
        allbad_she_measurements = table_formats.she_ksb_measurements.tf.init_table()

        for j in range(self.len_group):
            for i in range(len(self.she_measurements)):
                somebad_she_measurements.add_row(vals = {ksbm_tf.ID: i * self.len_group + j})
                allbad_she_measurements.add_row(vals = {ksbm_tf.ID: i * self.len_group + j})

        somebad_she_measurements[ksbm_tf.g1] = self.she_measurements_group[ksbm_tf.g1]
        somebad_she_measurements[ksbm_tf.g2] = self.she_measurements_group[ksbm_tf.g2]
        somebad_she_measurements[ksbm_tf.g1_err] = self.she_measurements_group[ksbm_tf.g1_err]
        somebad_she_measurements[ksbm_tf.g2_err] = self.she_measurements_group[ksbm_tf.g2_err]

        allbad_she_measurements[ksbm_tf.g1] = self.she_measurements_group[ksbm_tf.g1]
        allbad_she_measurements[ksbm_tf.g2] = self.she_measurements_group[ksbm_tf.g2]
        allbad_she_measurements[ksbm_tf.g1_err] = self.she_measurements_group[ksbm_tf.g1_err]
        allbad_she_measurements[ksbm_tf.g2_err] = self.she_measurements_group[ksbm_tf.g2_err]

        for i in [0, 3, 7, 8]:
            row = somebad_she_measurements[i]
            row[ksbm_tf.g1], row[ksbm_tf.g2] = np.nan, np.nan
            row[ksbm_tf.g1_err], row[ksbm_tf.g2_err] = np.inf, np.inf
        for i in range(len(self.she_measurements_group)):
            row = allbad_she_measurements[i]
            row[ksbm_tf.g1], row[ksbm_tf.g2] = np.nan, np.nan
            row[ksbm_tf.g1_err], row[ksbm_tf.g2_err] = np.inf, np.inf

        # Check that we don't crash by calling calculate_shear_bias_statistics
        g1_somebad_bias_stats, g2_somebad_bias_stats = calculate_shear_bias_statistics(
            somebad_she_measurements, self.details_table_group)
        g1_allbad_bias_stats, g2_allbad_bias_stats = calculate_shear_bias_statistics(
            allbad_she_measurements, self.details_table_group)

        # Check the statistics all have the expected weight

        assert_almost_equal(g1_somebad_bias_stats.w, self.len_group * 0.6 * self.ex_w1, decimal = 4)
        assert_almost_equal(g2_somebad_bias_stats.w, self.len_group * 0.6 * self.ex_w2, decimal = 4)

        assert not np.isnan(g1_somebad_bias_stats.xm)
        assert not np.isnan(g1_somebad_bias_stats.x2m)
        assert not np.isnan(g1_somebad_bias_stats.ym)
        assert not np.isnan(g1_somebad_bias_stats.xym)

        assert_almost_equal(g1_allbad_bias_stats.w, 0., decimal = 4)
        assert_almost_equal(g2_allbad_bias_stats.w, 0., decimal = 4)

    def test_calculate_shear_bias_statistics_group(self):
        """Try using the calculate_shear_bias_statistics function and check the results on grouped data.
        """
        g1_bias_stats, g2_bias_stats = calculate_shear_bias_statistics(
            self.she_measurements_group, self.details_table_group)

        g1_bias = BiasMeasurements(LinregressResults(g1_bias_stats))
        g2_bias = BiasMeasurements(LinregressResults(g2_bias_stats))

        assert_almost_equal(g1_bias.m, self.ex_m1, decimal = 4)
        assert_almost_equal(g1_bias.c, self.ex_c1, decimal = 4)

        assert_almost_equal(g2_bias.m, self.ex_m2, decimal = 4)
        assert_almost_equal(g2_bias.c, self.ex_c2, decimal = 4)

    def test_measure_statistics_from_args(self):
        """Try using the measure_statistics_from_args function and check the results.
        """

        # Set up the arguments object
        args = Args()
        args.workdir = self.workdir
        args.logdir = self.logdir
        args.details_table = "test_details_table.xml"
        args.shear_estimates = "test_she_measurements.xml"
        args.she_bias_statistics = "test_shear_statistics.xml"

        # Set up the files to be read in

        os.makedirs(os.path.join(args.workdir, "data"))

        details_table_filename = "test_details_table.fits"
        details_table_product = products.she_simulated_catalog.create_dpd_she_simulated_catalog(details_table_filename)
        write_xml_product(details_table_product, args.details_table, workdir = args.workdir)
        self.details_table.write(join(args.workdir, "data/" + details_table_filename), format = "fits")

        she_measurements_filename = "test_she_measurements.fits"
        she_measurements_product = products.she_measurements.create_she_measurements_product(
            KSB_filename = she_measurements_filename)
        write_xml_product(she_measurements_product, args.shear_estimates, workdir = args.workdir)
        self.she_measurements.write(join(args.workdir, "data/" + she_measurements_filename), format = "fits")

        # Call the function
        measure_statistics_from_args(args)

        # Read in and check the results
        she_bias_statistics_product = read_xml_product(join(args.workdir, args.she_bias_statistics))

        g1_bias_stats, g2_bias_stats = she_bias_statistics_product.get_KSB_bias_statistics(workdir = self.workdir)

        g1_bias = BiasMeasurements(LinregressResults(g1_bias_stats))
        g2_bias = BiasMeasurements(LinregressResults(g2_bias_stats))

        assert_almost_equal(g1_bias.m, self.ex_m1, decimal = 4)
        assert_almost_equal(g1_bias.c, self.ex_c1, decimal = 4)

        assert_almost_equal(g2_bias.m, self.ex_m2, decimal = 4)
        assert_almost_equal(g2_bias.c, self.ex_c2, decimal = 4)
