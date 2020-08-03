""" @file reconcile_shear.py

    Created 03 August 2020

    Primary execution loop for reconciling shear estimates into a per-tile catalog.
"""

__updated__ = "2020-08-03"

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

import numpy as np

logger = getLogger(__name__)

reconciliation_methods = {"Best": reconcile_best,
                          "InvVar": reconcile_inv_var}

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


def reconcile_shear_from_args(args):
    """ Primary function for performing shear reconciliation
    """

    # Determine the reconciliation method to use

    method = None

    if args.method is None:
        method = str(args.method)
        logger.info("Using reconciliation method: '" + str(method) + "', passed from command-line.")
    elif args.pipeline_config is not None:
        # Load in the pipeline configuration and see if the method is supplied there

        pipeline_config = read_reconciliation_config(args.pipeline_config,
                                                     workdir=args.workdir)

        if ReconciliationConfigKeys.REC_METHOD.value in pipeline_config:
            method = str(pipeline_config[ReconciliationConfigKeys.REC_METHOD.value])
            logger.info("Using reconciliation method: '" + str(method) + "', from pipeline configuration file.")

    if method is None:
        # If we get here, it isn't yet determined, so use the default
        method = default_reconciliation_method
        logger.info("Using default reconciliation method: '" + str(method) + "'.")

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

    # Loop over each method, and reconcile tables for that method
    for shear_estimation_method in shear_estimation_method_table_formats:

        # Create a catalog for reconciled measurements. In case of multiple measurements of the same object,
        # this will store the first temporarily and later be updated with the reconciled result
        # TODO - set tile ID in header
        reconciled_catalog = shear_estimation_method_table_initialisers[shear_estimation_method]()

        # Create a dict of objects needing reconciliation. Keys are object IDs, and values are tables containing
        # one row for each separate measurement
        ids_to_reconcile = {}

        # Create a set of IDs we've added to the table (faster to access than column indices)
        ids_in_reconciled_catalog = {}

        sem_tf = shear_estimation_method_table_formats[shear_estimation_method]

        # Loop through each table
        for estimates_table_filename in validated_shear_estimates_table_filenames[shear_estimation_method]:
            estimates_table = Table.read(args.workdir, estimates_table_filename)

            if not is_in_format(estimates_table, sem_tf):
                pass

    return

