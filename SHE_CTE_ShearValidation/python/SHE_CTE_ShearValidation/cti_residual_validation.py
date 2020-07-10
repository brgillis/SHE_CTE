""" @file cti_residual_validation.py

    Created 12 Feb 2018

    Function for performing validation of CTI residuals
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

from copy import deepcopy

from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_utility import is_in_format
import astropy.table
import astropy.wcs
from scipy.stats import linregress

import numpy as np

reg_pix = 1 # Readout register - believed to be y, but might have to fix later

readout_split = 4134.5 # Maximum pixel position for which values are read out downward (in stacked frame)

g1_fail_flag_offset = 1
g2_fail_flag_offset = 4

lower_fail_flag_offset = 1
upper_fail_flag_offset = 2

slope_fail_sigma = 5
intercept_fail_sigma = 5


def validate_cti_ellipticity_residual_bin(r, g):

    result_flag = 0

    # Check for both lower and upper half of detectors
    for flag, limits in ((lower_fail_flag_offset, (-1e99, readout_split)),
                         (upper_fail_flag_offset, (readout_split, 1e99)),
                         ):

        # Get rows within this half
        good_rows = np.logical_and(r >= limits[0], r < limits[1])
        r_binned = r[good_rows]
        g_binned = g[good_rows]

        # Perform a linear regression with errors
        slope, _intercept, _r, _p, slope_err = linregress(r_binned,
                                                          g_binned)

        # Check if we fall outside acceptable sigma for slope not being 0
        if abs(slope / slope_err) > slope_fail_sigma:
            result_flag += flag # Add flag for this half if it's bad

    return result_flag


def report_val_results(validation_results,
                       metaname,
                       shear_estimates_table):

    shear_estimates_table.meta[metaname] = repr(tuple(validation_results))


def validate_cti_ellipticity_residuals(shear_estimates_table,
                                       stacked_frame_wcs):

    # Check format of input
    if not isinstance(shear_estimates_table, astropy.table.Table):
        raise TypeError("shear_estimates_table must be an astropy.table.Table object")
    if not is_in_format(shear_estimates_table, sm_tf, strict=False):
        raise TypeError("shear_estimates_table is not in correct format.")

    if not isinstance(stacked_frame_wcs, astropy.wcs.WCS):
        raise TypeError("stacked_frame_wcs must be an astropy.wcs.WCS object")

    all_validated = True

    # Loop over different columns we're testing
    for colname, metaname, bins in ((sm_tf.snr, "CE_SN_VAL", (-1e99, 1e99)), # TODO make these magic values and store somewhere
                                    (sm_tf.sky_bg, "CE_SB_VAL", (-1e99, 1e99)),
                                    (sm_tf.re, "CE_RE_VAL", (-1e99, 1e99)),
                                    (sm_tf.color, "CE_CO_VAL", (-1e99, 1e99)),
                                    (sm_tf.epoch, "CE_EP_VAL", (-1e99, 1e99)),
                                    ):

        num_bins = len(bins) - 1

        col = shear_estimates_table[colname]

        # Set up the results array for this column
        validation_results = np.zeros(num_bins, dtype=int)

        # Loop over the bins for this column
        for bin_i in range(num_bins):
            bin_min = bins[bin_i]
            bin_max = bins[bin_i + 1]

            good_rows = np.logical_and(col >= bin_min, col < bin_max)
            r = shear_estimates_table[reg_pix][good_rows]

            # Test g1 and g2 for residuals, and flag appropriately

            validation_results[colname][bin_i] += (g1_fail_flag_offset *
                                                   validate_cti_ellipticity_residual_bin(r=r,
                                                                                         g=shear_estimates_table[sm_tf.g1][good_rows]))
            validation_results[colname][bin_i] += (g2_fail_flag_offset *
                                                   validate_cti_ellipticity_residual_bin(r=r,
                                                                                         g=shear_estimates_table[sm_tf.g2][good_rows]))
            # Report results in the header of the original table
            report_val_results(validation_results=validation_results,
                               metaname=metaname,
                               shear_estimates_table=shear_estimates_table)

            # Record if there is any failure
            if np.sum(validation_results) > 0:
                all_validated = False

    return all_validated
