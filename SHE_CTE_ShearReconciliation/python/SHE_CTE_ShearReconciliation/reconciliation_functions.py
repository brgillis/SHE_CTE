""" @file reconciliation_functions.py

    Created 04 August 2020

    Functions to handle different ways of reconciling different shear estimates.
"""

__updated__ = "2020-08-04"

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

import numpy as np
from SHE_PPT.utility import run_only_once

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
        - measurements_to_reconcile_table is sorted by weight
    
        Return
        ------
        None
    """

    measurements_to_reconcile_table.sort(sem_tf.weight)

    best_row = measurements_to_reconcile_table[-1]

    # Update the output row
    for colname in output_row.colnames:
        output_row[colname] = best_row[colname]

    return

# Data used for weight reconciliation
props_with_independent_errors = ("ra", "dec", "re", "flux", "bulge_frac", "snr", "sersic")
props_to_sum = ("nexp", "weight")
props_to_bitwise_or = ("fit_flags", "val_flags", )
props_to_copy = ("ID", "fit_class")
props_to_nan = ("chi2", "dof", "g1g2_covar", "g1g2_uncal_covar")

@run_only_once
def warn_missing_props(missing_props):
    warning_string = "The following properties have no defined method to combine them in the reconcile_weight function:"
    for prop in missing_props:
        warning_string += " " + prop + ","
    logger.warning(warning_string[:-1])

def reconcile_weight(measurements_to_reconcile_table,
                     output_row,
                     sem_tf,
                     *args, **kwargs):
    """ Reconciliation method which combines measurements based on their weights.
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
        - measurements_to_reconcile_table is sorted by weight
    
        Return
        ------
        None
    """

    # Determine the weight for each row
    if nexp in vars(sem_tf):
        weights = measurements_to_reconcile_table[sem_tf.weight] * measurements_to_reconcile_table[sem_tf.nexp]
    else:
        weights = measurements_to_reconcile_table[sem_tf.weight]
    
    # For any instances of NaN, set the weight to zero
    weights = np.where(np.logical_or(np.isnan(weights),np.isinf(weights)), 0, weights)
    tot_weight = np.sum(weights)
    highest_weight_index = np.argmax(weights)

    new_props = {}

    # TODO - implement special handling for shear, which has non-independent errors

    # Combine each value we can combine simply (where all errors are independent)
    for prop in props_with_independent_errors:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        new_props[colname] = np.sum(measurements_to_reconcile_table[colname]*weights)/tot_weight

        # If this property has an error, calculate that too
        prop_err = colname + "_err"
        if not prop_err in vars(sem_tf):
            continue
        colname_err = getattr(sem_tf, prop_err)
        new_props[colname_err] = np.sqrt(np.sum(np.pow(measurements_to_reconcile_table[colname]*weights,2))/tot_weight)

    # Combine properties we sum up
    for prop in props_to_bitwise_or:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        new_props[colname] = np.sum(measurements_to_reconcile_table[colname])

    # Combine properties bitwise or
    for prop in props_to_sum:
        if not prop in vars(sem_tf):
            continue
        colname = getattr(sem_tf, prop)
        new_props[colname] = np.vectorize(np.bitwise_or)(measurements_to_reconcile_table[colname])

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

    # Check for any missing properties, and warn and set them to NaN, while we update the output row
    for prop in measurements_to_reconcile_table.colnames:
        missing_props = []
        if prop in new_props:
            output_row[prop] = new_props[prop]
        else:
            missing_props.append(prop)
            output_row[prop] = np.NaN
        warn_missing_props(missing_props)

    return


