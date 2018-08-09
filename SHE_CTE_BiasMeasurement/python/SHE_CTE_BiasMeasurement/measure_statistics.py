""" @file measure_statistics.py

    Created 21 June 2018
    

    Executable for measuring necessary statistics on a set of shear
    measurements.
"""

__updated__ = "2018-08-09"

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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA

import os

from SHE_PPT import products
from SHE_PPT.file_io import read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics
from SHE_CTE_PipelineUtility.archive import archive_product
from astropy.table import Table


def measure_statistics_from_args(args):
    """Workhorse function for measuring necessary statistics to calculate shear estimation bias.
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureStatistics measure_statistics_from_args()')
    logger.debug('#')

    # Open the input files

    # Get the details table

    details_table_product = read_xml_product(os.path.join(args.workdir, args.details_table))
    details_table = Table.read(os.path.join(args.workdir, details_table_product.get_data_filename()))

    # Get the shear estimates product

    shear_estimates_table_product = read_xml_product(os.path.join(args.workdir, args.shear_estimates))

    # Initialise the output product

    shear_bias_statistics_product = products.shear_bias_stats.create_shear_bias_statistics_product()

    # Read in shear estimates, and calculate statistics for each method

    for method in mv.estimation_methods:

        try:

            estimates_table_filename = shear_estimates_table_product.get_method_filename(
                method)

            if estimates_table_filename is None or estimates_table_filename == "":
                continue

            estimates_table = Table.read(
                os.path.join(args.workdir, estimates_table_filename))

            # Calculate statistics
            if not method == "BFD":
                shear_bias_statistics = calculate_shear_bias_statistics(estimates_table, details_table)
            else:
                continue  # FIXME

            # Save these in the data product
            shear_bias_statistics_product.set_method_statistics(method, *shear_bias_statistics)

        except Exception as e:

            logger.warn("Failsafe exception block triggered with exception: " + str(e))
            shear_bias_statistics_product.set_method_statistics(method, None, None)

    # Write out the statistics product
    qualified_statistics_filename = os.path.join(args.workdir, args.shear_bias_statistics)
    write_xml_product(shear_bias_statistics_product, qualified_statistics_filename)

    # Try to archive the product
    if args.archive_dir is not None:
        try:
            archive_product(product_filename=args.shear_bias_statistics,
                            archive_dir=args.archive_dir,
                            workdir=args.workdir)
        except Exception as e:
            logger.warn("Failsafe exception block triggered when trying to save statistics product in archive. " +
                        "Exception was: " + str(e))

    logger.debug('# Exiting SHE_CTE_MeasureStatistics measure_statistics_from_args()')

    return
