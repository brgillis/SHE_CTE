""" @file statistics_calculation.py

    Created 22 June 2018

    Function to calculate bias statistics from tables of shear measurements and input details.
"""

__updated__ = "2018-06-22"

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
