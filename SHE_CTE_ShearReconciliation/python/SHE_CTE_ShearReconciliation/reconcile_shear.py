""" @file reconcile_shear.py

    Created 03 August 2020

    Primary execution loop for reconciling shear estimates into a per-tile catalog.
"""

__updated__ = "2021-08-18"

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
from SHE_PPT.constants.shear_estimation_methods import ShearEstimationMethods, D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS
from SHE_PPT.file_io import (write_xml_product, get_allowed_filename,
                             read_listfile, read_xml_product)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import ReconciliationConfigKeys
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf
from SHE_PPT.table_formats.she_lensmc_chains import tf as lmcc_tf
from SHE_PPT.table_utility import is_in_format
from SHE_PPT.utility import is_any_type_of_none
from astropy.table import Table

import SHE_CTE
from SHE_CTE_ShearReconciliation.chains_reconciliation_functions import (reconcile_chains_best,
                                                                         reconcile_chains_invvar,
                                                                         reconcile_chains_shape_weight,
                                                                         reconcile_chains_keep)
from SHE_CTE_ShearReconciliation.reconciliation_functions import (reconcile_best,
                                                                  reconcile_invvar,
                                                                  reconcile_shape_weight)


logger = getLogger(__name__)

reconciliation_methods = {"Best": reconcile_best,
                          "InvVar": reconcile_invvar,
                          "ShapeWeight": reconcile_shape_weight}

default_reconciliation_method = "InvVar"

chains_reconciliation_methods = {"Best": reconcile_chains_best,
                                 "InvVar": reconcile_chains_invvar,
                                 "ShapeWeight": reconcile_chains_shape_weight,
                                 "Keep": reconcile_chains_keep}

default_chains_reconciliation_method = "Keep"

assert default_reconciliation_method in reconciliation_methods


def store_object_info(new_row, existing_row, ids_to_reconcile, sem_tf, optional_columns):
    """ Stores an object's info (from its row) in the proper table within the ids_to_reconcile dict, so it
        can be handled later.
    """

    object_id = new_row[sem_tf.ID]
    assert object_id == existing_row[sem_tf.ID]
    assert len(new_row) == len(existing_row)

    # Is this the first conflict with this ID?
    if not object_id in ids_to_reconcile:
        # First conflict, so add the id with a table using the existing row, ensuring the right column order
        t = sem_tf.init_table(optional_columns=optional_columns)[optional_columns]
        t.add_row(existing_row)
        ids_to_reconcile[object_id] = t
    else:
        t = ids_to_reconcile[object_id]

    # Add the new row to the table for this ID's conflicts
    t.add_row(new_row)


def reconcile_tables(shear_estimates_tables,
                     shear_estimation_method,
                     object_ids_in_tile,
                     reconciliation_function,
                     workdir=None):

    sem_tf = D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS[shear_estimation_method]

    # We'll properly create the reconciled catalog once we know what optional columns to include in it
    reconciled_catalog = None

    # Create a dict of objects needing reconciliation. Keys are object IDs, and values are tables containing
    # one row for each separate measurement
    ids_to_reconcile = {}

    # Create a set of IDs we've added to the table (faster to access than column indices)
    ids_in_reconciled_catalog = set()

    # Loop through each table
    for estimates_table in shear_estimates_tables:

        if is_any_type_of_none(estimates_table):
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

        # We can now create the reconciled catalog
        if reconciled_catalog is None:

            # Create a catalog for reconciled measurements. In case of multiple measurements of the same object,
            # this will store the first temporarily and later be updated with the reconciled result.
            # Note that we reorder the columns to match input here through the indexing at the end
            # TODO - set tile ID in header

            optional_columns = estimates_table.colnames
            reconciled_catalog = sem_tf.init_table(optional_columns=optional_columns)[optional_columns]

            # Set up the object ID to be used as an index for this catalog
            reconciled_catalog.add_index(sem_tf.ID)

        # Loop over the rows of the table
        for row in estimates_table:
            object_id = row[sem_tf.ID]

            # Skip if this ID isn't in the MER catalog for the Tile
            if object_id not in object_ids_in_tile:
                continue

            # Check if this ID is already in the reconciled catalog
            if object_id in ids_in_reconciled_catalog:
                store_object_info(new_row=row,
                                  existing_row=reconciled_catalog.loc[object_id],
                                  ids_to_reconcile=ids_to_reconcile,
                                  sem_tf=sem_tf,
                                  optional_columns=optional_columns,)

            else:
                # Otherwise, add it to the reconciled catalog
                reconciled_catalog.add_row(row)
                ids_in_reconciled_catalog.add(object_id)

        # End looping through rows of this table
    # End looping through tables

    # Now, we need to perform the reconciliation of each id
    for object_id in ids_to_reconcile:
        reconciliation_function(measurements_to_reconcile_table=ids_to_reconcile[object_id],
                                output_row=reconciled_catalog.loc[object_id],
                                sem_tf=sem_tf)
    return reconciled_catalog


def reconcile_chains(chains_tables,
                     object_ids_in_tile,
                     chains_reconciliation_function,
                     workdir=None):

    # We'll properly create the reconciled chains catalog once we know what optional columns to include in it
    reconciled_chains = None

    # Create a dict of objects needing reconciliation. Keys are object IDs, and values are tables containing
    # one row for each separate measurement
    ids_to_reconcile = {}

    # Create a set of IDs we've added to the table (faster to access than column indices)
    ids_in_reconciled_chains = set()

    # Loop through each table
    for chains_table in chains_tables:

        if is_any_type_of_none(chains_table):
            continue

        if isinstance(chains_table, str):
            if workdir is None:
                raise ValueError("If a filename is passed to reconcile_chains (\"" + chains_table + "\"), " +
                                 "the workdir must also be supplied.")
            # It's a filename, so load it in
            qualified_chains_table_filename = os.path.join(workdir, chains_table)
            chains_table = Table.read(qualified_chains_table_filename)
        else:
            # If it's not a filename, store the string version of the table for logging purposes
            qualified_chains_table_filename = str(chains_table)

        # Ensure it's in the right format
        if not is_in_format(chains_table, lmcc_tf, verbose=True):
            raise ValueError("Table " + qualified_chains_table_filename + " is not in expected table format (" +
                             lmcc_tf.m.table_format + "). See log for details of error.")

        # We can now create the reconciled catalog
        if reconciled_chains is None:

            # Create a catalog for reconciled measurements. In case of multiple measurements of the same object,
            # this will store the first temporarily and later be updated with the reconciled result.
            # Here we need to reorder the columns (done by the indexing at the end) to match the input catalog
            # TODO - set tile ID in header

            optional_columns = chains_table.colnames
            reconciled_chains = lmcc_tf.init_table(optional_columns=optional_columns)[optional_columns]

            # Set up the object ID to be used as an index for this catalog
            reconciled_chains.add_index(lmcc_tf.ID)

        # Loop over the rows of the table
        for row in chains_table:
            object_id = row[lmcc_tf.ID]

            # Skip if this ID isn't in the MER catalog for the Tile
            if object_id not in object_ids_in_tile:
                continue

            # Check if this ID is already in the reconciled catalog
            if object_id in ids_in_reconciled_chains:
                store_object_info(new_row=row,
                                  existing_row=reconciled_chains.loc[object_id],
                                  ids_to_reconcile=ids_to_reconcile,
                                  sem_tf=lmcc_tf,
                                  optional_columns=optional_columns,)

            else:
                # Otherwise, add it to the reconciled catalog
                reconciled_chains.add_row(row)
                ids_in_reconciled_chains.add(object_id)

        # End looping through rows of this table
    # End looping through tables

    # Now, we need to perform the reconciliation of each id
    for object_id in ids_to_reconcile:
        extra_rows = chains_reconciliation_function(chains_to_reconcile_table=ids_to_reconcile[object_id],
                                                    output_row=reconciled_chains.loc[object_id])
        if extra_rows is not None:
            for extra_row in extra_rows:
                reconciled_chains.add_row(extra_row)

    return reconciled_chains


def reconcile_shear_from_args(args):
    """ Primary function for performing shear reconciliation
    """

    # Read in the pipeline config if supplied
    pipeline_config = args.she_reconciliation_config

    # Get the reconciliation method to use from the pipeline_config

    method = pipeline_config[ReconciliationConfigKeys.REC_METHOD]
    logger.info("Using reconciliation method: '" + str(method) + "', from pipeline configuration file.")

    if not method in reconciliation_methods:
        allowed_method_str = ""
        for allowed_method in reconciliation_methods:
            allowed_method_str += "\n" + allowed_method
        raise ValueError("Reconciliation method " + method + " is not recognized. Allowed methods are:" +
                         allowed_method_str)
    reconciliation_function = reconciliation_methods[method]

    # Get the chains reconciliation method to use from the pipeline config

    chains_method = pipeline_config[ReconciliationConfigKeys.CHAINS_REC_METHOD]
    logger.info("Using chains reconciliation method: '" +
                str(chains_method) + "', from pipeline configuration file.")

    # Check we're using a valid chains reconciliation method
    if not chains_method in chains_reconciliation_methods:
        allowed_chains_method_str = ""
        for allowed_chains_method in chains_reconciliation_methods:
            allowed_chains_method_str += "\n" + allowed_chains_method
        raise ValueError("Chains reconciliation method " + chains_method + " is not recognized. Allowed methods are:" +
                         allowed_chains_method_str)
    chains_reconciliation_function = chains_reconciliation_methods[chains_method]

    # Load in the final catalog from MER to get the IDs of objects in this tile
    mer_final_catalog_product = read_xml_product(args.mer_final_catalog, workdir=args.workdir)
    mer_final_catalog_table = Table.read(os.path.join(args.workdir,
                                                      mer_final_catalog_product.get_filename()))
    object_ids_in_tile = frozenset(mer_final_catalog_table[mfc_tf.ID])

    # Load in the filenames of the tables for each method

    validated_shear_estimates_table_filenames = {}
    for shear_estimation_method in ShearEstimationMethods:
        validated_shear_estimates_table_filenames[shear_estimation_method] = []

    # Start by loading in the data products in the listfile
    she_validated_measurements_filename_list = read_listfile(os.path.join(args.workdir,
                                                                          args.she_validated_measurements_listfile))

    # Loop over each product (representing a different observation)
    for she_validated_measurements_filename in she_validated_measurements_filename_list:
        she_validated_measurement_product = read_xml_product(xml_filename=she_validated_measurements_filename,
                                                             workdir=args.workdir)

        # Get the table filename for each method, and if it isn't None, add it to that method's list
        for shear_estimation_method in ShearEstimationMethods:
            estimates_table_filename = she_validated_measurement_product.get_method_filename(shear_estimation_method)
            if not is_any_type_of_none(estimates_table_filename):
                validated_shear_estimates_table_filenames[shear_estimation_method].append(estimates_table_filename)

    # Create a data product for the output
    reconciled_catalog_product = products.she_reconciled_measurements.create_dpd_she_reconciled_measurements(
        spatial_footprint=mer_final_catalog_product)

    # Get the filenames of the chains tables to be reconciled
    lensmc_chains_table_filenames = []
    she_lensmc_chains_filename_list = read_listfile(os.path.join(args.workdir,
                                                                 args.she_lensmc_chains_listfile))
    for she_lensmc_chains_product_filename in she_lensmc_chains_filename_list:
        she_lensmc_chains_product = read_xml_product(she_lensmc_chains_product_filename, workdir=args.workdir)
        she_lensmc_chains_table_filename = she_lensmc_chains_product.get_filename()
        if not is_any_type_of_none(she_lensmc_chains_table_filename):
            lensmc_chains_table_filenames.append(she_lensmc_chains_table_filename)

    # Get the tile_id and set it for the new product
    tile_id = mer_final_catalog_product.Data.TileIndex
    reconciled_catalog_product.Data.TileIndex = tile_id

    # Loop over each method, and reconcile tables for that method
    for shear_estimation_method in ShearEstimationMethods:

        reconciled_catalog = reconcile_tables(shear_estimates_tables=validated_shear_estimates_table_filenames[shear_estimation_method],
                                              shear_estimation_method=shear_estimation_method,
                                              object_ids_in_tile=object_ids_in_tile,
                                              reconciliation_function=reconciliation_function,
                                              workdir=args.workdir)

        # If we don't have any data, don't write a table at all
        if reconciled_catalog is None or len(reconciled_catalog) == 0:

            logger.info("No data to output for method " + shear_estimation_method + ".")

            reconciled_catalog_product.set_method_filename(method=shear_estimation_method,
                                                           filename=None)

        else:

            # The output table is now finalized, so output it and store the filename in the output data product
            reconciled_catalog_filename = get_allowed_filename(type_name="REC-SHM-" + shear_estimation_method.name,
                                                               instance_id=str(tile_id),
                                                               extension=".fits",
                                                               version=SHE_CTE.__version__)

            qualified_reconciled_catalog_filename = os.path.join(args.workdir, reconciled_catalog_filename)

            logger.info(f"Outputting reconciled catalog for method {shear_estimation_method.value} to " +
                        qualified_reconciled_catalog_filename)

            reconciled_catalog.write(qualified_reconciled_catalog_filename)

            reconciled_catalog_product.set_method_filename(method=shear_estimation_method,
                                                           filename=reconciled_catalog_filename)

    # End looping over methods

    # Output the finalized data product to the desired filename
    logger.info("Outputting reconciled catalog data product to " +
                os.path.join(args.workdir, args.she_reconciled_measurements))
    write_xml_product(reconciled_catalog_product, args.she_reconciled_measurements, workdir=args.workdir)

    # Now reconcile the chains
    reconciled_chains_catalog = reconcile_chains(chains_tables=lensmc_chains_table_filenames,
                                                 object_ids_in_tile=object_ids_in_tile,
                                                 chains_reconciliation_function=chains_reconciliation_function,
                                                 workdir=args.workdir)

    # Create the output product
    reconciled_chains_product = products.she_reconciled_lensmc_chains.create_dpd_she_reconciled_lensmc_chains(
        spatial_footprint=mer_final_catalog_product)
    reconciled_chains_product.Data.TileIndex = tile_id

    # If we don't have any data, create an empty table
    if reconciled_chains_catalog is None:

        logger.warning("No reconciled chains data to output.")

        reconciled_chains_catalog = lmcc_tf.init_table()

    # The output table is now finalized, so output it and store the filename in the output data product
    reconciled_chains_catalog_filename = get_allowed_filename(type_name="REC-LMC-CHAINS",
                                                              instance_id=str(tile_id),
                                                              extension=".fits",
                                                              version=SHE_CTE.__version__)

    qualified_reconciled_chains_catalog_filename = os.path.join(args.workdir, reconciled_chains_catalog_filename)

    logger.info("Outputting reconciled chains catalog to " +
                qualified_reconciled_chains_catalog_filename)

    reconciled_chains_catalog.write(qualified_reconciled_chains_catalog_filename)

    reconciled_chains_product.set_filename(reconciled_chains_catalog_filename)

    # Output the finalized chains data product to the desired filename
    logger.info("Outputting reconciled catalog data product to " +
                os.path.join(args.workdir, args.she_reconciled_lensmc_chains))
    write_xml_product(reconciled_chains_product, args.she_reconciled_lensmc_chains, workdir=args.workdir)

    logger.debug("# Exiting reconcile_shear_from_args() successfully.")
