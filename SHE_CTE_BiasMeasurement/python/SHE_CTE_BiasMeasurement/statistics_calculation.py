""" @file statistics_calculation.py

    Created 22 June 2018

    Function to calculate bias statistics from tables of shear measurements and input details.
"""
from numpy.testing.utils import assert_allclose

__updated__ = "2018-06-29"

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

from astropy import table

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_PPT.logging import getLogger
from SHE_PPT.math import get_linregress_statistics
from SHE_PPT.table_formats.details import tf as datf
from SHE_PPT.table_formats.shear_estimates import tf as setf
import numpy as np


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
        group_id = row[datf.group_ID]
        if not group_id in groups:
            groups[group_id] = [row]
        else:
            groups[group_id].append(row)

    # Check if we actually need to combine anything
    if len(combined_table) == len(groups):
        # No need to combine; just return original table
        return combined_table

    # Create a new table for compressed data
    compressed_table = table.Table(names=[datf.group_ID, datf.g1, datf.g2,
                                          setf.g1, setf.g2, setf.g1_err, setf.g2_err],
                                   dtype=[datf.dtypes[datf.group_ID], datf.dtypes[datf.g1], datf.dtypes[datf.g2],
                                          setf.dtypes[setf.g1], setf.dtypes[setf.g2],
                                          setf.dtypes[setf.g1_err], setf.dtypes[setf.g2_err]])

    for group_id in groups:

        # Get a list of only good measurements
        all_rows = groups[group_id]
        good_rows = []

        for row in all_rows:
            g1, g2, g1_err, g2_err = row[setf.g1], row[setf.g2], row[setf.g1_err], row[setf.g2_err]
            if (g1 > -2 and g1 < 2 and g2 > -2 and g2 < 2 and
                    g1_err > 0 and g1_err < 1e99 and g2_err > 0 and g2_err < 1e99):
                good_rows.append(row)

        # Sort the data into numpy arrays
        num_good_rows = len(good_rows)
        data = {}
        for name in [datf.g1, datf.g2, setf.g1, setf.g2, setf.g1_err, setf.g2_err]:
            data[name] = np.zeros(num_good_rows)
            for i in range(num_good_rows):
                data[name][i] = good_rows[i][name]

        # Check we have a non-zero number of good values
        if num_good_rows == 0:
            continue

        # Check all real values are close. If not, we shouldn't be grouping
        assert_allclose(data[datf.g1], data[datf.g1][0])
        assert_allclose(data[datf.g2], data[datf.g2][0])

        # Calculate weighted means of the measurements
        g1_weight = data[setf.g1_err]**-2
        g2_weight = data[setf.g2_err]**-2

        total_g1_weight = g1_weight.sum()
        total_g2_weight = g2_weight.sum()

        if total_g1_weight <= 0 or total_g2_weight <= 0:
            raise ValueError("Bad weights in combining shear measurements.")

        mean_g1 = np.sum(data[setf.g1] * g1_weight) / total_g1_weight
        mean_g2 = np.sum(data[setf.g2] * g2_weight) / total_g2_weight

        mean_g1_err = np.sqrt(1 / np.sum(g1_weight))
        mean_g2_err = np.sqrt(1 / np.sum(g2_weight))

        # Add these to the compressed table
        compressed_table.add_row(vals={datf.group_ID: group_id,
                                       datf.g1: data[datf.g1][0],
                                       datf.g2: data[datf.g2][0],
                                       setf.g1: mean_g1,
                                       setf.g2: mean_g2,
                                       setf.g1_err: mean_g1_err,
                                       setf.g2_err: mean_g2_err})

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

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureStatistics calculate_shear_bias_statistics()')
    logger.debug('#')

    # Create a combined table, joined on galaxy ID
    if setf.ID != datf.ID:
        details_table.rename_column(datf.ID, setf.ID)
    combined_table = table.join(estimates_table, details_table, keys=setf.ID)
    if setf.ID != datf.ID:
        details_table.rename_column(setf.ID, datf.ID)

    # Compress the table on group ID to properly handle shape noise cancellation
    compressed_table = compress_details_and_measurements(combined_table)

    # Get stats for both g1 and g2
    bias_stats = []
    for g_est_colname, g_err_colname, g_true_colname in ((setf.g1, setf.g1_err, datf.g1,),
                                                         (setf.g2, setf.g2_err, datf.g2,),):
        lx = combined_table[g_true_colname].data
        ly = combined_table[g_est_colname].data
        ly_err = combined_table[g_err_colname].data

        bias_stats.append(get_linregress_statistics(lx, ly, ly_err))

    logger.debug('# Exiting SHE_CTE_MeasureStatistics calculate_shear_bias_statistics()')

    return tuple(bias_stats)
