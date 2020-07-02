""" @file statistics_calculation.py

    Created 22 June 2018

    Function to calculate bias statistics from tables of shear measurements and input details.
"""

__updated__ = "2020-07-02"

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
from SHE_PPT.math import BFDSumStatistics, get_bfd_sum_statistics
from astropy import table

from SHE_BFD_CalculateMoments.bfd import BfdPqrs
from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_PPT.table_formats.she_bfd_moments import tf as bfdm_tf
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_formats.she_simulated_catalog import tf as simc_tf


def calculate_bfd_shear_bias_statistics(estimates_table, details_table):
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
    logger.debug('# Entering SHE_CTE_MeasureStatistics calculate_bfd_shear_bias_statistics()')
    logger.debug('#')

    # If there are no rows in the estimates table, exit early will an empty statistics object
    if len(estimates_table) == 0:
        stats = BFDSumStatistics()
        stats.b1 = 0
        stats.b2 = 0
        stats.b3 = 0
        stats.b4 = 0
        stats.A11 = 0
        stats.A12 = 0
        stats.A13 = 0
        stats.A14 = 0
        stats.A22 = 0
        stats.A23 = 0
        stats.A24 = 0
        stats.A33 = 0
        stats.A34 = 0
        stats.A44 = 0

    # Create a combined table, joined on galaxy ID
    if sm_tf.ID != simc_tf.ID:
        details_table.rename_column(simc_tf.ID, sm_tf.ID)
    combined_table = table.join(estimates_table, details_table, keys=sm_tf.ID, metadata_conflicts='silent')
    if sm_tf.ID != simc_tf.ID:
        details_table.rename_column(sm_tf.ID, simc_tf.ID)

    # Get stats for file
    # start a bfd_pqrs instance

    pqr = BfdPqrs(pqr=combined_table[bfdm_tf.bfd_pqr])
    sums = pqr.get_sums(g1_true=combined_table[simc_tf.g1], g2_true=combined_table[simc_tf.g2])
    bfd_stats = get_bfd_sum_statistics(sums)

    logger.debug('# Exiting SHE_CTE_MeasureStatistics calculate_shear_bias_statistics()')

    return bfd_stats
