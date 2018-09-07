""" @file statistics_calculation.py

    Created 22 June 2018

    Function to calculate bias statistics from tables of shear measurements and input details.
"""
from SHE_PPT.logging import getLogger
from SHE_PPT.math import get_linregress_statistics, LinregressStatistics, BFDSumStatistics, get_bfd_sum_statistics
from SHE_PPT.table_formats.details import tf as datf
from SHE_PPT.table_formats.shear_estimates import tf as setf
from SHE_PPT.table_formats.bfd_moments import tf as bfdf
from astropy import table
from numpy.testing.utils import assert_allclose

from SHE_BFD_CalculateMoments.bfd import bfd_pqrs
from SHE_CTE_BiasMeasurement import magic_values as mv
import numpy as np


__updated__ = "2018-07-02"

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
    if setf.ID != datf.ID:
        details_table.rename_column(datf.ID, setf.ID)
    combined_table = table.join(estimates_table, details_table, keys=setf.ID,metadata_conflicts='silent')
    if setf.ID != datf.ID:
        details_table.rename_column(setf.ID, datf.ID)

    # Get stats for file
    # start a bfd_pqrs instance

    pqr = bfd_pqrs(pqr=combined_table[bfdf.bfd_pqr])
    sums = pqr.get_sums(g1_true=combined_table[datf.g1], g2_true=combined_table[datf.g2])
    bfd_stats = get_bfd_sum_statistics(sums)

    logger.debug('# Exiting SHE_CTE_MeasureStatistics calculate_shear_bias_statistics()')

    return bfd_stats,None # MuST return as tuple to fit with other methods
