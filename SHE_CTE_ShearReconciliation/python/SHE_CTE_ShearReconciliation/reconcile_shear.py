""" @file reconcile_shear.py

    Created 03 August 2020

    Primary execution loop for reconciling shear estimates into a per-tile catalog.
"""

__updated__ = "2020-09-03"

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

import os

from SHE_PPT import products
from SHE_PPT.file_io import (write_xml_product, get_allowed_filename, get_data_filename,
                             read_listfile, find_file, read_xml_product)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import ReconciliationConfigKeys, read_reconciliation_config
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf
from SHE_PPT.table_formats.she_bfd_moments import initialise_bfd_moments_table, tf as bfdm_tf
from SHE_PPT.table_formats.she_ksb_measurements import initialise_ksb_measurements_table, tf as ksbm_tf
from SHE_PPT.table_formats.she_lensmc_measurements import initialise_lensmc_measurements_table, tf as lmcm_tf
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_formats.she_momentsml_measurements import initialise_momentsml_measurements_table, tf as mmlm_tf
from SHE_PPT.table_formats.she_regauss_measurements import initialise_regauss_measurements_table, tf as regm_tf
from SHE_PPT.table_utility import is_in_format
from astropy.table import Table
from galsim import shear

import SHE_CTE
from SHE_CTE_ShearReconciliation.reconciliation_functions import reconcile_best, reconcile_invvar, reconcile_shape_weight
import numpy as np

logger = getLogger(__name__)

reconciliation_methods = {"Best": reconcile_best,
                          "InvVar": reconcile_invvar,
                          "ShapeWeight": reconcile_shape_weight}

default_reconciliation_method = "InvVar"

shear_estimation_method_table_formats = {"KSB": ksbm_tf,
                                         "REGAUSS": regm_tf,
                                         "MomentsML": mmlm_tf,
                                         "LensMC": lmcm_tf,
                                         "BFD": bfdm_tf}

shear_estimation_method_table_initialisers = {"KSB": initialise_ksb_measurements_table,
                                              "REGAUSS": initialise_regauss_measurements_table,
                                              "MomentsML": initialise_momentsml_measurements_table,
                                              "LensMC": initialise_lensmc_measurements_table,
                                              "BFD": initialise_bfd_moments_table}

assert default_reconciliation_method in reconciliation_methods


def store_object_info(new_row, existing_row, ids_to_reconcile, sem_tf, table_initialiser):
    """ Stores an object's info (from its row) in the proper table within the ids_to_reconcile dict, so it
        can be handled later.
    """

    id = new_row[sem_tf.ID]
    assert id == existing_row[sem_tf.ID]

    # Is this the first conflict with this ID?
    if not id in ids_to_reconcile:
        # First conflict, so add the id with a table using the existing row
        t = table_initialiser()
        t.add_row(existing_row)
        ids_to_reconcile[id] = t
    else:
        t = ids_to_reconcile[id]

    # Add the new row to the table for this ID's conflicts
    t.add_row(new_row)

    return


def reconcile_tables(shear_estimates_tables,
                     shear_estimation_method,
                     object_ids_in_tile,
                     reconciliation_function,
                     workdir=None):

    sem_tf = shear_estimation_method_table_formats[shear_estimation_method]
    table_initialiser = shear_estimation_method_table_initialisers[shear_estimation_method]

    # Create a catalog for reconciled measurements. In case of multiple measurements of the same object,
    # this will store the first temporarily and later be updated with the reconciled result
    # TODO - set tile ID in header
    reconciled_catalog = table_initialiser()

    # Set up the object ID to be used as an index for this catalog
    reconciled_catalog.add_index(sem_tf.ID)

    # Create a dict of objects needing reconciliation. Keys are object IDs, and values are tables containing
    # one row for each separate measurement
    ids_to_reconcile = {}

    # Create a set of IDs we've added to the table (faster to access than column indices)
    ids_in_reconciled_catalog = set()

    # Loop through each table
    for estimates_table in shear_estimates_tables:

        if estimates_table is None or estimates_table == "None" or estimates_table == "data/None":
            continue

        if isinstance(estimates_table, str):
            if workdir is None:
                raise ValueError("If a filename is passed to reconcile_tables (\"" + estimates_table + "\"), " +
                                 "the workdir must also be supplied.")
            # It's a filename, so load it in
            qualified_estimates_table_filename = os.path.join(workdir, estimates_table)
            estimates_table = Table.read(qualified_estimates_table_filename)
        else:
            qualified_estimates_table_filename = str(estimates_table)

        # Ensure it's in the right format
        if not is_in_format(estimates_table, sem_tf, verbose=True):
            raise ValueError("Table " + qualified_estimates_table_filename + " is not in expected table format (" +
                             sem_tf.m.table_format + "). See log for details of error.")

        # Loop over the rows of the table
        for row in estimates_table:
            id = row[sem_tf.ID]

            # Skip if this ID isn't in the MER catalog for the Tile
            if id not in object_ids_in_tile:
                continue

            # Check if this ID is already in the reconciled catalog
            if id in ids_in_reconciled_catalog:
                store_object_info(new_row=row,
                                  existing_row=reconciled_catalog.loc[id],
                                  ids_to_reconcile=ids_to_reconcile,
                                  sem_tf=sem_tf,
                                  table_initialiser=table_initialiser)

            else:
                # Otherwise, add it to the reconciled catalog
                reconciled_catalog.add_row(row)
                ids_in_reconciled_catalog.add(id)

        # End looping through rows of this table
    # End looping through tables

    # Now, we need to perform the reconciliation of each id
    for id in ids_to_reconcile:
        reconciliation_function(measurements_to_reconcile_table=ids_to_reconcile[id],
                                output_row=reconciled_catalog.loc[id],
                                sem_tf=sem_tf)
    return reconciled_catalog


def reconcile_shear_from_args(args):
    """ Primary function for performing shear reconciliation
    """

    # Determine the reconciliation method to use

    method = None

    if args.method is not None:
        method = str(args.method)
        logger.info("Using reconciliation method: '" + str(method) + "', passed from command-line.")
    elif args.she_reconciliation_config is not None and args.she_reconciliation_config is not "None":
        # Load in the pipeline configuration and see if the method is supplied there

        pipeline_config = read_reconciliation_config(args.she_reconciliation_config,
                                                     workdir=args.workdir)

        if ReconciliationConfigKeys.REC_METHOD.value in pipeline_config:
            method = str(pipeline_config[ReconciliationConfigKeys.REC_METHOD.value])
            logger.info("Using reconciliation method: '" + str(method) + "', from pipeline configuration file.")

    if method is None:
        # If we get here, it isn't yet determined, so use the default
        method = default_reconciliation_method
        logger.info("Using default reconciliation method: '" + str(method) + "'.")

    if not method in reconciliation_methods:
        allowed_method_str = ""
        for allowed_method in reconciliation_methods:
            allowed_method_str += "\n" + allowed_method
        raise ValueError("Reconciliation method " + method + " is not recognized. Allowed methods are:" +
                         allowed_method_str)
    reconciliation_function = reconciliation_methods[method]

    # Load in the final catalog from MER to get the IDs of objects in this tile
    mer_final_catalog_product = read_xml_product(args.mer_final_catalog, workdir=args.workdir)
    mer_final_catalog_table = Table.read(os.path.join(args.workdir,
                                                      mer_final_catalog_product.get_filename()))
    object_ids_in_tile = frozenset(mer_final_catalog_table[mfc_tf.ID])

    # Load in the filenames of the tables for each method

    validated_shear_estimates_table_filenames = {}
    for shear_estimation_method in shear_estimation_method_table_formats:
        validated_shear_estimates_table_filenames[shear_estimation_method] = []

    # Start by loading in the data products in the listfile
    she_validated_measurements_filename_list = read_listfile(os.path.join(args.workdir,
                                                                          args.she_validated_measurements_listfile))

    # Loop over each product (representing a different observation)
    for she_validated_measurements_filename in she_validated_measurements_filename_list:
        she_validated_measurement_product = read_xml_product(xml_filename=she_validated_measurements_filename,
                                                             workdir=args.workdir)

        # Get the table filename for each method, and if it isn't None, add it to that method's list
        for shear_estimation_method in shear_estimation_method_table_formats:
            estimates_table_filename = she_validated_measurement_product.get_method_filename(shear_estimation_method)
            if estimates_table_filename is not None:
                validated_shear_estimates_table_filenames[shear_estimation_method].append(estimates_table_filename)

    # Create a data product for the output
    reconciled_catalog_product = products.she_reconciled_measurements.create_dpd_she_reconciled_measurements(
        spatial_footprint=mer_final_catalog_product)

    # Get the tile_id and set it for the new product
    tile_id = mer_final_catalog_product.Data.TileIndex
    reconciled_catalog_product.Data.TileIndex = tile_id

    # Loop over each method, and reconcile tables for that method
    for shear_estimation_method in shear_estimation_method_table_formats:

        reconciled_catalog = reconcile_tables(shear_estimates_tables=validated_shear_estimates_table_filenames[shear_estimation_method],
                                              shear_estimation_method=shear_estimation_method,
                                              object_ids_in_tile=object_ids_in_tile,
                                              reconciliation_function=reconciliation_function,
                                              workdir=args.workdir)

        # If we don't have any data, don't write a table at all
        if len(reconciled_catalog) == 0:

            reconciled_catalog_product.set_method_filename(method=shear_estimation_method,
                                                           filename=None)

        else:

            # The output table is now finalized, so output it and store the filename in the output data product
            reconciled_catalog_filename = get_allowed_filename(type_name="REC-SHM-" + shear_estimation_method.upper(),
                                                               instance_id=str(tile_id),
                                                               extension=".fits",
                                                               version=SHE_CTE.__version__)
            reconciled_catalog.write(os.path.join(args.workdir, reconciled_catalog_filename))

            reconciled_catalog_product.set_method_filename(method=shear_estimation_method,
                                                           filename=reconciled_catalog_filename)

    # End looping over methods

    # Output the finalized data product to the desired filename
    write_xml_product(reconciled_catalog_product, args.she_reconciled_measurements, workdir=args.workdir)

    logger.debug("# Exiting reconcile_shear_from_args() successfully.")

    return
