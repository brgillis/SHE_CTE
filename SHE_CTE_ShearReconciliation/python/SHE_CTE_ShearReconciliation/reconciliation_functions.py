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
