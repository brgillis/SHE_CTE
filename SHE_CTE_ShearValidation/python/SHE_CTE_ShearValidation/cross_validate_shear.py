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

import math
from os.path import join

import SHE_CTE
from SHE_CTE_ShearValidation import magic_values as mv
from SHE_PPT import products
from SHE_PPT.file_io import (read_xml_product, write_xml_product,
                             get_allowed_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.bfd_moments import tf as bfdtf
from SHE_PPT.table_formats.shear_estimates import tf as setf
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from astropy.io import fits
from astropy.table import Table
import numpy as np


setfs = {"KSB": setf,
         "REGAUSS": setf,
         "MomentsML": setf,
         "LensMC": setf,
         "BFD": bfdtf}

products.calibration_parameters.init()
products.shear_estimates.init()
products.shear_validation_stats.init()
products.validated_shear_estimates.init()


def cross_validate_shear_estimates(primary_shear_estimates_table,
                                   other_shear_estimates_tables=None,
                                   shear_validation_statistics_table=None):
    """
        Stub for validating shear estimates against validation statistics. Presently
        sets everything to "pass".
    """

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
    primary_shear_estimates_table.meta[setf.m.validated] = 1

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

    if not isinstance(shear_estimates_prod, products.shear_estimates.dpdShearMeasurement):
        raise ValueError("Shear estimates product from " + join(args.workdir, args.shear_estimates_product)
                         + " is invalid type.")

    primary_shear_estimates_table = None
    other_shear_estimates_tables = {}

    for method in setfs:

        filename = shear_estimates_prod.get_method_filename(method)

        if filename is not None and filename != "None":
            shear_estimates_table = Table.read(join(args.workdir, filename), format='fits')
            if not is_in_format(shear_estimates_table, setfs[method]):
                raise ValueError("Shear estimates table from " +
                                 join(args.workdir, filename) + " is in invalid format.")
        else:
            shear_estimates_table = None

        if method == args.primary_method:
            primary_shear_estimates_table = shear_estimates_table
        else:
            other_shear_estimates_tables[method] = shear_estimates_table

    # Shear validation statistics - Disabled until this exists

    if False:

        logger.info("Reading" + dry_label + " shear validation statistics...")

        shear_validation_stats_prod = read_xml_product(join(args.workdir, args.shear_validation_statistics_table))
        if not isinstance(shear_validation_stats_prod, products.shear_validation_stats.DpdSheShearValidationStatsProduct):
            raise ValueError("Shear validation statistics product from " + join(args.workdir, args.shear_validation_stats_product)
                             + " is invalid type.")

        shear_validation_stats_filename = shear_validation_stats_prod.get_filename()

        shear_validation_statistics_table = Table.read(join(args.workdir, shear_validation_stats_filename))

        if not is_in_format(shear_validation_statistics_table, setf):
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

    validated_shear_estimates_filename = get_allowed_filename("VAL_SHM", "0", ".fits",
                                                              version=SHE_CTE.__version__)

    validated_shear_estimates_prod = products.validated_shear_estimates.create_validated_shear_estimates_product(
        validated_shear_estimates_filename)

    write_xml_product(validated_shear_estimates_prod,
                      join(args.workdir, args.cross_validated_shear_estimates_product))

    primary_shear_estimates_table.write(validated_shear_estimates_filename)

    logger.info("Finished" + dry_label + " shear validation.")

    return
