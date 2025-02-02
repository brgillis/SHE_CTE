""" @file statistics_calculation.py

    Created 22 June 2018

    Function to calculate bias statistics from tables of shear measurements and input details.
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
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

from SHE_PPT.logging import getLogger
from SHE_PPT.math import get_linregress_statistics, LinregressStatistics
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksbm_tf
from SHE_PPT.table_formats.she_lensmc_measurements import tf as lmcm_tf
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_formats.she_momentsml_measurements import tf as mmlm_tf
from SHE_PPT.table_formats.she_regauss_measurements import tf as regm_tf
from SHE_PPT.table_formats.she_simulated_catalog import tf as simc_tf
from astropy import table
from numpy.testing import assert_allclose

import numpy as np

table_formats = {"she.ksbMeasurements": ksbm_tf,
                 "she.regaussMeasurements": regm_tf,
                 "she.momentsmlMeasurements": mmlm_tf,
                 "she.lensmcMeasurements": lmcm_tf, }


def compress_details_and_measurements(combined_table):
    """
    @brief
        Compress measurements when shape noise cancellation was used, to combine
        measurements made on the same input values

    @param combined_table <astropy.table.Table>

    @return compressed_table <astropy.table.Table>
    """

    groups = {}

    # First go through and see what we need to combine
    for row in combined_table:
        group_id = row[simc_tf.group_ID]
        if not group_id in groups:
            groups[group_id] = [row]
        else:
            groups[group_id].append(row)

    # Check if we actually need to combine anything
    if len(combined_table) == len(groups):
        # No need to combine; just return original table
        return combined_table

    tf = table_formats[combined_table.meta[sm_tf.m.fits_def]]

    # Create a new table for compressed data
    compressed_table = table.Table(names=[simc_tf.group_ID, simc_tf.g1, simc_tf.g2,
                                          tf.g1, tf.g2, tf.g1_err, tf.g2_err],
                                   dtype=[simc_tf.dtypes[simc_tf.group_ID], simc_tf.dtypes[simc_tf.g1], simc_tf.dtypes[simc_tf.g2],
                                          tf.dtypes[tf.g1], tf.dtypes[tf.g2],
                                          tf.dtypes[tf.g1_err], tf.dtypes[tf.g2_err]])

    for group_id in groups:

        # Get a list of only good measurements
        all_rows = groups[group_id]
        good_rows = []

        for row in all_rows:
            g1, g2, g1_err, g2_err = row[tf.g1], row[tf.g2], row[tf.g1_err], row[tf.g2_err]
            if (g1 > -2 and g1 < 2 and g2 > -2 and g2 < 2 and
                    g1_err > 0 and g1_err < 1e99 and g2_err > 0 and g2_err < 1e99):
                good_rows.append(row)

        # Sort the data into numpy arrays
        num_good_rows = len(good_rows)
        data = {}
        for name in [simc_tf.g1, simc_tf.g2, tf.g1, tf.g2, tf.g1_err, tf.g2_err]:
            data[name] = np.zeros(num_good_rows)
            for i in range(num_good_rows):
                data[name][i] = good_rows[i][name]

        # Check we have a non-zero number of good values
        if num_good_rows == 0:
            continue

        # Check all real values are close. If not, we shouldn't be grouping
        assert_allclose(data[simc_tf.g1], data[simc_tf.g1][0])
        assert_allclose(data[simc_tf.g2], data[simc_tf.g2][0])

        # Calculate unweighted means of the measurements but retain full weight
        # Must be unweighted to retain benefits of shape noise cancellation
        g1_weight = data[tf.g1_err] ** -2
        g2_weight = data[tf.g2_err] ** -2

        total_g1_weight = g1_weight.sum()
        total_g2_weight = g2_weight.sum()

        if total_g1_weight <= 0 or total_g2_weight <= 0:
            raise ValueError("Bad weights in combining shear measurements.")

        combined_g1_err = total_g1_weight ** -0.5
        combined_g2_err = total_g2_weight ** -0.5

        mean_g1 = np.mean(data[tf.g1])
        mean_g2 = np.mean(data[tf.g2])

        # Add these to the compressed table
        compressed_table.add_row(vals={simc_tf.group_ID: group_id,
                                       simc_tf.g1: data[simc_tf.g1][0],
                                       simc_tf.g2: data[simc_tf.g2][0],
                                       tf.g1: mean_g1,
                                       tf.g2: mean_g2,
                                       tf.g1_err: combined_g1_err,
                                       tf.g2_err: combined_g2_err})

    return compressed_table


def calculate_shear_bias_statistics(estimates_table, details_table):
    """Calculates shear bias statistics from the provided shear estimates table and input details table.

    Parameters
    ----------
    estimates_table : astropy.table.Table
        Table containing shear estimates
    details_table : astropy.table.Table
        Table containing true galaxy details
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureStatistics calculate_shear_bias_statistics()')
    logger.debug('#')

    # If there are no rows in the estimates table, exit early will an empty statistics object
    if len(estimates_table) == 0:
        g1_stats = LinregressStatistics()
        g2_stats = LinregressStatistics()
        for stats in g1_stats, g2_stats:
            stats.w = 0
            stats.xm = 0
            stats.x2m = 0
            stats.ym = 0
            stats.xym = 0
        return g1_stats, g2_stats

    tf = table_formats[estimates_table.meta[sm_tf.m.fits_def]]

    # Create a combined table, joined on galaxy ID
    if tf.ID != simc_tf.ID:
        details_table.rename_column(simc_tf.ID, tf.ID)
    combined_table = table.join(estimates_table, details_table, keys=tf.ID)
    if tf.ID != simc_tf.ID:
        details_table.rename_column(tf.ID, simc_tf.ID)
    combined_table.meta[tf.m.fits_def] = estimates_table.meta[sm_tf.m.fits_def]

    # Compress the table on group ID to properly handle shape noise cancellation
    compressed_table = compress_details_and_measurements(combined_table)

    # Get stats for both g1 and g2
    bias_stats = []
    for g_est_colname, g_err_colname, g_true_colname in ((tf.g1, tf.g1_err, simc_tf.g1,),
                                                         (tf.g2, tf.g2_err, simc_tf.g2,),):
        lx = compressed_table[g_true_colname].data
        ly = compressed_table[g_est_colname].data
        ly_err = compressed_table[g_err_colname].data

        bias_stats.append(get_linregress_statistics(lx, ly, ly_err))

    logger.debug('# Exiting SHE_CTE_MeasureStatistics calculate_shear_bias_statistics()')

    return tuple(bias_stats)
