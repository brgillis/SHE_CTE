""" @file measure_statistics.py

    Created 21 June 2018
    

    Executable for measuring necessary statistics on a set of shear
    measurements.
"""

__updated__ = "2018-07-27"

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

from os.path import join

from SHE_PPT import products
from SHE_PPT.file_io import read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics
from astropy.table import Table


products.details.init()
products.shear_estimates.init()
products.shear_bias_stats.init()


def measure_statistics_from_args(args):
    """Workhorse function for measuring necessary statistics to calculate shear estimation bias.
    """

    logger = getLogger(mv.logger_name)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MeasureStatistics measure_statistics_from_args()')
    logger.debug('#')

    # Open the input files

    # Get the details table

    details_table_product = read_xml_product(join(args.workdir, args.details_table))
    details_table = Table.read(join(args.workdir, details_table_product.get_data_filename()))

    # Get the shear estimates product

    shear_estimates_table_product = read_xml_product(join(args.workdir, args.shear_estimates))

    # Initialise the output product

    shear_bias_statistics_product = products.shear_bias_stats.create_shear_bias_statistics_product()

    # Read in shear estimates, and calculate statistics for each method

    for method in mv.estimation_methods:

        estimates_table_filename = shear_estimates_table_product.get_method_filename(
            method)

        if estimates_table_filename is None or estimates_table_filename == "":
            continue

        estimates_table = Table.read(
            join(args.workdir, estimates_table_filename))

        # Calculate statistics
        if not method == "BFD":
            shear_bias_statistics = calculate_shear_bias_statistics(estimates_table, details_table)
        else:
            continue

        # Save these in the data product
        shear_bias_statistics_product.set_method_statistics(method, *shear_bias_statistics)

    # Write out the statistics product
    write_xml_product(shear_bias_statistics_product, join(args.workdir, args.shear_bias_statistics))

    logger.debug('# Exiting SHE_CTE_MeasureStatistics measure_statistics_from_args()')

    return
