""" @file measure_statistics.py

    Created 21 June 2018
    
    Executable for measuring necessary statistics on a set of shearmeasurements.
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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA

import os

from SHE_PPT import products
from SHE_PPT.constants.shear_estimation_methods import ShearEstimationMethods
from SHE_PPT.file_io import read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import archive_product, CalibrationConfigKeys
from astropy.table import Table

from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics


def measure_statistics_from_args(args):
    """Workhorse function for measuring necessary statistics to calculate shear estimation bias.
    """

    logger = getLogger(__name__)

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

    shear_bias_statistics_product = products.she_bias_statistics.create_she_bias_statistics_product()

    # Read in shear estimates, and calculate statistics for each method

    for method in ShearEstimationMethods:

        try:

            estimates_table_filename = shear_estimates_table_product.get_method_filename(
                method)

            if estimates_table_filename is None or estimates_table_filename == "":
                continue

            estimates_table = Table.read(
                os.path.join(args.workdir, estimates_table_filename))

            # Calculate statistics
            she_bias_statistics = calculate_shear_bias_statistics(estimates_table, details_table)

            # Save these in the data product
            shear_bias_statistics_product.set_method_bias_statistics(method, she_bias_statistics, workdir=args.workdir)

        except Exception as e:

            logger.warning("Failsafe exception block triggered with exception: " + str(e))
            shear_bias_statistics_product.set_method_bias_statistics(method, None, workdir=args.workdir)

    # Write out the statistics product
    write_xml_product(shear_bias_statistics_product, args.she_bias_statistics, workdir=args.workdir)

    # Try to archive the product

    # First get the pipeline config so we can figure out where to archive it
    pipeline_config = args.pipeline_config

    archive_dir = pipeline_config[CalibrationConfigKeys.MS_ARCHIVE_DIR]
    webdav_dir = pipeline_config[CalibrationConfigKeys.MS_WEBDAV_ARCHIVE]
    webdav_archive = pipeline_config[CalibrationConfigKeys.MS_WEBDAV_DIR]

    # If we're archiving with webdav, determine its mount dir and the full archive directory
    if webdav_archive and archive_dir is not None:
        full_archive_dir = os.path.join(webdav_dir, archive_dir)
    else:
        full_archive_dir = archive_dir

    if archive_dir is not None:
        try:
            archive_product(product_filename=args.she_bias_statistics,
                            archive_dir=full_archive_dir,
                            workdir=args.workdir)
        except Exception as e:
            logger.warning("Failsafe exception block triggered when trying to save statistics product in archive. " +
                           "Exception was: " + str(e))

    logger.debug('# Exiting SHE_CTE_MeasureStatistics measure_statistics_from_args()')
