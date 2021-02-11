""" @file reconciliation_functions.py

    Created 04 August 2020

    Functions to handle different ways of reconciling different shear estimates.
"""

__updated__ = "2021-02-11"

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
from SHE_PPT.utility import run_only_once
from astropy.io.ascii.core import masked

import numpy as np

logger = getLogger(__name__)


def reconcile_best(measurements_to_reconcile_table,
                   output_row,
                   sem_tf,
                   *args, **kwargs):
    """ Reconciliation method which selects the best (highest-weight) measurement.

        Parameters
        ----------
        measurements_to_reconcile_table : astropy.table.Table
            The table containing rows of different measurements of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled measurement
        sem_tf : Table Format
            The table format of the measurements table
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled measurement

        Return
        ------
        None
    """

    measurements_to_reconcile_table.sort(sem_tf.weight)

    # Find the row with the highest weight, ignoring NaN and Inf weights
    best_index = None
    best_weight = -1e99

    for i in range(len(measurements_to_reconcile_table)):
        weight = measurements_to_reconcile_table[i][sem_tf.weight]
        if not (np.isnan(weight) or np.isinf(weight)) and weight > best_weight:
            best_index = i
            best_weight = weight

    # If no good rows are available, just use the first
    if best_index is None:
        best_index = 0

    best_row = measurements_to_reconcile_table[best_index]

    # Update the output row
    for colname in output_row.colnames:
        output_row[colname] = best_row[colname]

    return


def reconcile_shape_weight(measurements_to_reconcile_table,
                           output_row,
                           sem_tf,
                           *args, **kwargs):
    """ Reconciliation method which combines measurements based on supplied shape
        measurement weights.


        Parameters
        ----------
        measurements_to_reconcile_table : astropy.table.Table
            The table containing rows of different measurements of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled measurement
        sem_tf : Table Format
            The table format of the measurements table
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled measurement

        Return
        ------
        None
    """

    weights = measurements_to_reconcile_table[sem_tf.shape_weight]

    return reconcile_weight(measurements_to_reconcile_table=measurements_to_reconcile_table,
                            output_row=output_row,
                            sem_tf=sem_tf,
                            weights=weights,
                            *args, **kwargs)


def reconcile_invvar(measurements_to_reconcile_table,
                     output_row,
                     sem_tf,
                     *args, **kwargs):
    """ Reconciliation method which combines measurements based on inverse-shape-variance
        weighting.


        Parameters
        ----------
        measurements_to_reconcile_table : astropy.table.Table
            The table containing rows of different measurements of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled measurement
        sem_tf : Table Format
            The table format of the measurements table
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled measurement

        Return
        ------
        None
    """

    weights = 1 / measurements_to_reconcile_table[sem_tf.e_var]

    return reconcile_weight(measurements_to_reconcile_table=measurements_to_reconcile_table,
                            output_row=output_row,
                            sem_tf=sem_tf,
                            weights=weights,
                            *args, **kwargs)


# Data used for weight reconciliation
props_with_independent_errors = ("g1", "g2", "g1_uncal", "g2_uncal", "ra", "dec", "re", "flux",
                                 "bulge_frac", "snr", "sersic", "unmasked_fraction")
props_to_sum = ("nexp",)
props_to_bitwise_or = ("fit_flags", "val_flags",)
props_to_copy = ("ID", "fit_class")
props_to_nan = ("chi2", "dof", "g1g2_covar", "g1g2_uncal_covar")
props_to_ignore = ("rec_flags")


@run_only_once
def warn_missing_props(missing_props):
    warning_string = "The following properties have no defined method to combine them in the reconcile_weight function:"
    for prop in missing_props:
        warning_string += " " + prop + ","
    logger.warning(warning_string[:-1])


def reconcile_weight(measurements_to_reconcile_table,
                     output_row,
                     sem_tf,
                     weights,
                     *args, **kwargs):
    """ Reconciliation method which combines measurements based on given weights.
        Parameters
        ----------
        measurements_to_reconcile_table : astropy.table.Table
            The table containing rows of different measurements of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled measurement
        sem_tf : Table Format
            The table format of the measurements table
        weights : Iterable <float>
            The weights to use for each row of the table
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled measurement

        Return
        ------
        None
    """

    # For any instances of NaN, set the weight to zero
    weights = np.where(np.logical_or(np.isnan(weights), np.isinf(weights)), 0, weights)
    # Also for any negative weights (we do this in a separate step to avoid a warning)
    weights = np.where(weights < 0, 0, weights)

    # Make a mask of objects with zero weight
    m = weights <= 0
    masked_weights = np.ma.masked_array(weights, m)

    tot_weight = weights.sum()
    highest_weight_index = np.argmax(weights)

    if tot_weight <= 0:
        logger.warning("Total weight for object " + str(measurements_to_reconcile_table[sem_tf.ID]) + " with:"
                       "weight = " + str(measurements_to_reconcile_table[sem_tf.weight]) + " is not positive. " +
                       "Will output copy of 'best' row.")
        return reconcile_best(measurements_to_reconcile_table,
                              output_row,
                              em_tf,
                              weights,
                              *args, **kwargs)

    new_props = {}

    # Combine each value we can combine simply (where all errors are independent)
    for prop in props_with_independent_errors:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        masked_column = np.ma.masked_array(measurements_to_reconcile_table[colname], m)
        new_props[colname] = (masked_column * masked_weights).sum() / tot_weight

        # If this property has an error, calculate that too
        prop_err = prop + "_err"
        if not prop_err in vars(sem_tf):
            continue
        colname_err = getattr(sem_tf, prop_err)
        masked_column_err = np.ma.masked_array(measurements_to_reconcile_table[colname_err], m)
        new_props[colname_err] = np.sqrt((np.power(masked_column_err * masked_weights, 2)).sum()) / tot_weight

    # Combine properties we sum up
    for prop in props_to_sum:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        masked_column = np.ma.masked_array(measurements_to_reconcile_table[colname], m)
        new_props[colname] = masked_column.sum()

    # Combine properties bitwise or
    for prop in props_to_bitwise_or:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        masked_column = np.ma.masked_array(measurements_to_reconcile_table[colname], m)
        new_props[colname] = np.bitwise_or.reduce(masked_column[~m], axis=0)

    # Copy the highest-weight property when we just copy
    for prop in props_to_copy:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        new_props[colname] = measurements_to_reconcile_table[colname][highest_weight_index]

    # Set to NaN properties that can't be sensibly combined in any way
    for prop in props_to_nan:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        new_props[colname] = np.NaN

    # Figure out what the shape noise is from the shape errors and weights, and use it to calculate weight
    shape_noise_var = (np.ma.masked_array(measurements_to_reconcile_table[sem_tf.shape_noise], m)**2).mean()
    new_props[sem_tf.weight] = 1. / (0.5 * (new_props[sem_tf.g1_err] ** 2 +
                                            new_props[sem_tf.g2_err] ** 2) + shape_noise_var)

    # Check for any missing properties, and warn and set them to NaN, while we update the output row
    missing_props = []
    for prop in measurements_to_reconcile_table.colnames:
        if prop in new_props:
            output_row[prop] = new_props[prop]
        elif not prop in props_to_ignore:
            missing_props.append(prop)
            try:
                output_row[prop] = np.NaN
            except ValueError:
                # If it's an int, can't set to NaN, so set to -99 as a failure indicator
                output_row[prop] = -99
    if len(missing_props) > 0:
        warn_missing_props(missing_props)

    return
