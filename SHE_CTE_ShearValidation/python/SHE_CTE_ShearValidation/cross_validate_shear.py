""" @file cross_validate_shear.py

    Created 12 Oct 2017

    Function for performing shear validation.
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

from os.path import join

from SHE_PPT import products
from SHE_PPT.file_io import (read_xml_product, write_xml_product,
                             get_allowed_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.she_bfd_moments import tf as bfdm_tf
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_utility import is_in_format
from astropy.table import Table

import SHE_CTE


sm_tfs = {"KSB": sm_tf,
          "REGAUSS": sm_tf,
          "MomentsML": sm_tf,
          "LensMC": sm_tf,
          "BFD": bfdm_tf}


def cross_validate_shear_estimates(primary_shear_estimates_table,
                                   other_shear_estimates_tables=None,
                                   shear_validation_statistics_table=None):
    """
        Stub for validating shear estimates against validation statistics. Presently
        sets everything to "pass".
    """

    # Skip if we don't have the primary table
    if primary_shear_estimates_table is None:
        return

    if other_shear_estimates_tables is None:
        other_shear_estimates_tables = []

    logger = getLogger(__name__)
    logger.debug("Entering validate_shear_estimates")

    # Compare primary table to all others

    for other_shear_estimates_table in other_shear_estimates_tables:
        # TODO Do something to compare with primary
        pass

    # TODO analyse comparisons somehow

    # For now, just say it passed
    primary_shear_estimates_table.meta[sm_tf.m.valid] = 1

    logger.debug("Exiting validate_shear_estimates")

    return


def cross_validate_shear(args, dry_run=False):
    """
        Main function for shear validation.
    """

    logger = getLogger(__name__)

    if dry_run:
        dry_label = " dry"
    else:
        dry_label = ""

    # Load in the files in turn to make sure there aren't any issues with them.

    # Shear estimates product

    logger.info("Reading" + dry_label + " shear estimates product...")

    shear_estimates_prod = read_xml_product(join(args.workdir, args.shear_estimates_product))

    if not isinstance(shear_estimates_prod, products.she_measurements.dpdSheMeasurements):
        raise ValueError("Shear estimates product from " + join(args.workdir, args.shear_estimates_product)
                         + " is invalid type.")

    primary_shear_estimates_table = None
    other_shear_estimates_tables = {}

    for method in sm_tfs:

        filename = shear_estimates_prod.get_method_filename(method)

        if filename is not None and filename != "None" and filename != "data/None":
            shear_estimates_table = Table.read(join(args.workdir, filename), format='fits')
            if not is_in_format(shear_estimates_table, sm_tfs[method], strict=False):
                logger.warning("Shear estimates table from " +
                               join(args.workdir, filename) + " is in invalid format.")
                continue
        else:
            shear_estimates_table = None
            if method == args.primary_method:
                logger.warning("Primary shear estimates file is not available.")

        if method == args.primary_method:
            primary_shear_estimates_table = shear_estimates_table
        else:
            other_shear_estimates_tables[method] = shear_estimates_table

    # Shear validation statistics - Disabled until this exists

    if False:

        logger.info("Reading" + dry_label + " shear validation statistics...")

        shear_validation_stats_prod = read_xml_product(join(args.workdir, args.shear_validation_statistics_table))
        if not isinstance(shear_validation_stats_prod, products.she_expected_shear_validation_statistics.DpdSheExpectedShearValidationStatistics):
            raise ValueError("Shear validation statistics product from " + join(args.workdir, args.shear_validation_stats_product)
                             + " is invalid type.")

        shear_validation_stats_filename = shear_validation_stats_prod.get_filename()

        shear_validation_statistics_table = Table.read(join(args.workdir, shear_validation_stats_filename))

        if not is_in_format(shear_validation_statistics_table, sm_tf):
            raise ValueError("Shear validation statistics table from " +
                             join(args.workdir, filename) + " is in invalid format.")
    else:
        shear_validation_statistics_table = None

    # Perform the validation
    cross_validate_shear_estimates(primary_shear_estimates_table=primary_shear_estimates_table,
                                   other_shear_estimates_tables=other_shear_estimates_tables,
                                   shear_validation_statistics_table=shear_validation_statistics_table)

    # Set up output product

    logger.info("Generating" + dry_label + " validated shear estimates...")

    method_table_filenames = {}

    # Write out the primary table
    if primary_shear_estimates_table is None:
        validated_primary_shear_estimates_filename = "None"
    else:
        validated_primary_shear_estimates_filename = get_allowed_filename("VAL-SHM-" + args.primary_method.upper(), "0",
                                                                          ".fits", version=SHE_CTE.__version__)
        primary_shear_estimates_table.write(
            join(args.workdir, validated_primary_shear_estimates_filename))
        logger.info("Wrote primary shear estimates table to " + validated_primary_shear_estimates_filename)
    method_table_filenames[args.primary_method] = validated_primary_shear_estimates_filename

    # Write out the other tables
    for method in other_shear_estimates_tables:
        method_table = other_shear_estimates_tables[method]
        if method_table is None:
            method_table_filename = "None"
        else:
            method_table_filename = get_allowed_filename("VAL-SHM-" + method.upper(), "0", ".fits",
                                                         version=SHE_CTE.__version__)
            method_table.write(join(args.workdir, method_table_filename))
            logger.info("Wrote shear estimates table to " + method_table_filename)
        method_table_filenames[method] = method_table_filename

    validated_shear_estimates_prod = products.she_validated_measurements.create_dpd_she_validated_measurements(
        spatial_footprint=shear_estimates_prod)

    for method in method_table_filenames:
        validated_shear_estimates_prod.set_method_filename(method, method_table_filenames[method])

    write_xml_product(validated_shear_estimates_prod, args.cross_validated_shear_estimates_product,
                      workdir=args.workdir)

    logger.info("Wrote cross-validated shear measurements product to " + args.cross_validated_shear_estimates_product)

    logger.info("Finished" + dry_label + " shear validation.")

    return
