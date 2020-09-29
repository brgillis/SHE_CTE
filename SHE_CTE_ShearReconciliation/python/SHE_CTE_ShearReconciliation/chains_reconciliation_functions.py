""" @file chains_reconciliation_functions.py

    Created 03 September 2020

    Functions to handle different ways of reconciling different chains.
"""

__updated__ = "2020-09-29"

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
from SHE_PPT.table_formats.she_lensmc_chains import tf as lmcc_tf, len_chain
from SHE_PPT.utility import run_only_once
from astropy.io.ascii.core import masked

import numpy as np

logger = getLogger(__name__)


def reconcile_chains_best(chains_to_reconcile_table,
                          output_row,
                          *args, **kwargs):
    """ Reconciliation method which selects the best (highest-weight) chains.

        Parameters
        ----------
        chains_to_reconcile_table : astropy.table.Table
            The table containing rows of different chains of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled chains
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled chains

        Return
        ------
        None
    """

    # Find the row with the highest weight, ignoring NaN and Inf weights
    best_index = None
    best_weight = -1e99

    for i in range(len(chains_to_reconcile_table)):
        weight = chains_to_reconcile_table[i][lmcc_tf.weight]
        if not (np.isnan(weight) or np.isinf(weight)) and weight > best_weight:
            best_index = i
            best_weight = weight

    # If no good rows are available, just use the first
    if best_index is None:
        best_index = 0

    best_chains_row = chains_to_reconcile_table[best_index]

    # Update the output row
    for colname in output_row.colnames:
        output_row[colname] = best_chains_row[colname]

    return


def reconcile_chains_invvar(chains_to_reconcile_table,
                            output_row,
                            *args, **kwargs):
    """ Reconciliation method which combines chains using the inverse of the supplied
        variance as the weight.


        Parameters
        ----------
        chains_to_reconcile_table : astropy.table.Table
            The table containing rows of different chains of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled chains
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled chains

        Return
        ------
        None
    """

    weights = 1. / (0.5 * chains_to_reconcile_table[lmcc_tf.e_var])

    return reconcile_chains_weight(chains_to_reconcile_table=chains_to_reconcile_table,
                                   output_row=output_row,
                                   weights=weights,
                                   *args, **kwargs)


def reconcile_chains_shape_weight(chains_to_reconcile_table,
                                  output_row,
                                  *args, **kwargs):
    """ Reconciliation method which combines chains based on supplied shape
        measurement weights.


        Parameters
        ----------
        chains_to_reconcile_table : astropy.table.Table
            The table containing rows of different chains of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled chains
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled chains

        Return
        ------
        None
    """

    weights = chains_to_reconcile_table[lmcc_tf.shape_weight]

    return reconcile_chains_weight(chains_to_reconcile_table=chains_to_reconcile_table,
                                   output_row=output_row,
                                   weights=weights,
                                   *args, **kwargs)


# Data used for weight reconciliation
props_with_independent_errors = ("g1", "g2", "ra", "dec", "re", "flux",
                                 "bulge_frac", "snr", "lr", "shape_noise", "e_var")
props_to_sum = ("nexp",)
props_to_bitwise_or = ("fit_flags", "val_flags",)
props_to_copy = ("ID", "fit_class")
props_to_nan = ("chi2", "dof", "acc")


@run_only_once
def warn_missing_props(missing_props):
    warning_string = "The following properties have no defined method to combine them in the reconcile_weight function:"
    for prop in missing_props:
        warning_string += " " + prop + ","
    logger.warning(warning_string[:-1])


def reconcile_chains_weight(chains_to_reconcile_table,
                            output_row,
                            weights,
                            *args, **kwargs):
    """ Reconciliation method which combines chains based on given weights.

        Parameters
        ----------
        chains_to_reconcile_table : astropy.table.Table
            The table containing rows of different chains of the same object
        output_row : row of astropy.table.Row
            The row to which to output the reconciled chains
        weights : astropy.table.Column (floating-point data type)
            The weights to use for each row of the table
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the reconciled chains

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
    chain_m = np.repeat(m, len_chain)

    masked_weights = np.ma.masked_array(weights, m)
    masked_chain_weights = np.repeat(masked_weights, len_chain)

    len_weights = len(weights)

    tot_weight = masked_weights.sum()
    highest_weight_index = np.argmax(weights)

    if tot_weight <= 0:
        logger.warning("Total weight for object " + str(chains_to_reconcile_table[lmcc_tf.ID]) + " with:"
                       "weight = " + str(chains_to_reconcile_table[lmcc_tf.weight]) + " is not positive. " +
                       "Will output copy of 'best' row.")
        return reconcile_best(chains_to_reconcile_tablem,
                              output_row,
                              *args, **kwargs)

    new_props = {}

    # Combine each value we can combine simply (where all errors are independent)
    for prop in props_with_independent_errors:
        colname = getattr(lmcc_tf, prop)
        if not colname in chains_to_reconcile_table.colnames:
            continue
        column = chains_to_reconcile_table[colname]
        if len(column.ravel()) == len_weights:
            masked_column = np.ma.masked_array(chains_to_reconcile_table[colname], m, fill=0)
            new_props[colname] = (masked_column * masked_weights).sum(axis=0) / tot_weight
        elif len(column.ravel()) == len_weights * len_chain:
            masked_column = np.ma.masked_array(chains_to_reconcile_table[colname].data.ravel(), chain_m, fill=0)
            new_props[colname] = (masked_column * masked_chain_weights).reshape((len_weights,
                                                                                 len_chain)).sum(axis=0) / tot_weight
        else:
            raise ValueError("Chains length in table (" + str(len(column.ravel()) / len_weights) + ") not equal to expected length (" +
                             str(len_chain) + ").")

    # Combine properties we sum up
    for prop in props_to_sum:
        colname = getattr(lmcc_tf, prop)
        if not colname in chains_to_reconcile_table.colnames:
            continue
        column = chains_to_reconcile_table[colname]
        masked_column = np.ma.masked_array(chains_to_reconcile_table[colname], m, fill=0)
        new_props[colname] = masked_column.sum(axis=0)

    # Combine properties bitwise or
    for prop in props_to_bitwise_or:
        colname = getattr(lmcc_tf, prop)
        if not colname in chains_to_reconcile_table.colnames:
            continue
        masked_column = np.ma.masked_array(chains_to_reconcile_table[colname], m)
        new_props[colname] = np.bitwise_or.reduce(masked_column[~m], axis=0)

    # Copy the highest-weight property when we just copy
    for prop in props_to_copy:
        colname = getattr(lmcc_tf, prop)
        if not colname in chains_to_reconcile_table.colnames:
            continue
        new_props[colname] = chains_to_reconcile_table[colname][highest_weight_index]

    # Set to NaN properties that can't be sensibly combined in any way
    for prop in props_to_nan:
        colname = getattr(lmcc_tf, prop)
        if not colname in chains_to_reconcile_table.colnames:
            continue
        new_props[colname] = np.NaN

    # Use the shape noise to calculate new weight

    shape_weight = masked_weights.sum(axis=0)
    # new_props[lmcc_tf.shape_weight] = shape_weight
    new_props[lmcc_tf.weight] = 1. / (1. / shape_weight + new_props[lmcc_tf.shape_noise]**2)

    # Check for any missing properties, and warn and set them to NaN, while we update the output row
    for prop in chains_to_reconcile_table.colnames:
        missing_props = []
        if prop in new_props:
            output_row[prop] = new_props[prop]
        else:
            missing_props.append(prop)
            output_row[prop] = np.NaN
        if len(missing_props) > 0:
            warn_missing_props(missing_props)

    return


def reconcile_chains_keep(chains_to_reconcile_table,
                          output_row,
                          *args, **kwargs):
    """ Reconciliation method which keeps all chains, adjusting weights appropriately.

        Parameters
        ----------
        chains_to_reconcile_table : astropy.table.Table
            The table containing rows of different chains of the same object
        output_row : row of astropy.table.Row
            The row to which to output the first chain (with adjusted weights)
        *args, **kwards
            Needed in case a different reconciliation method has an expanded interface

        Side-effects
        ------------
        - output_row is updated with the first chain (with adjusted weights)

        Return
        ------
        List of rows, one for each observation beyond the first, matching the format of output_row
    """

    weights = 1 / (0.5 * chains_to_reconcile_table[lmcc_tf.e_var])

    # For any instances of NaN, set the weight to zero
    weights = np.where(np.logical_or(np.isnan(weights), np.isinf(weights)), 0, weights)
    # Also for any negative weights (we do this in a separate step to avoid a warning)
    weights = np.where(weights < 0, 0, weights)

    # Make a mask of objects with zero weight
    m = weights <= 0
    masked_weights = np.ma.masked_array(weights, m)
    tot_weight = masked_weights.sum()

    # Calculate the normalized weights
    masked_normed_weights = masked_weights / masked_weights.sum()

    # Use the shape noise to calculate new total weight

    shape_noise_var = (
        masked_weights * np.ma.masked_array(chains_to_reconcile_table[lmcc_tf.shape_noise], m)**2).sum() / tot_weight

    total_shape_weight = masked_weights.sum(axis=0)
    total_shear_weight = 1. / (1. / total_shape_weight + shape_noise_var)

    first_row = True
    extra_rows = []

    for i in range(len(chains_to_reconcile_table)):

        row = chains_to_reconcile_table[i]
        row[lmcc_tf.weight] = total_shear_weight * masked_normed_weights[i]

        if not first_row:
            extra_rows.append(row)
            continue

        # For the first row, update the existing row in the catalog
        for colname in row.colnames:
            output_row[colname] = row[colname]

    return extra_rows
