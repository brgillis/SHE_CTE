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
from shutil import copyfile

from SHE_PPT import products
from SHE_PPT.file_io import read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_BiasMeasurement.statistics_calculation import calculate_shear_bias_statistics
from astropy.table import Table


def archive_product(product_filename, archive_dir, workdir):
    """ Copies an already-written data product to an archive directory.

        Parameters
        ----------
        product_filename : string
            The (unqualified) name of the product to copy
        archive_dir : string
            The root of the archive directory (note, the most-specific part of the workdir path (normally "workspace")
            will be added after this to keep separate runs from conflicting).
        workdir : string
            The working directory for this task

    """

    logger = getLogger(mv.logger_name)

    # Start by figuring out the subdirectory to store it in, based off of the workdir we're using
    subdir = os.path.split(workdir)[1]
    full_archive_dir = os.path.join(archive_dir, subdir)

    # The filename will likely also contain a subdir, so figure that out
    product_subpath = os.path.split(product_filename)[0]

    # Make the directory to store it in
    full_archive_subdir = os.path.join(full_archive_dir, product_subpath)
    full_archive_datadir = os.path.join(full_archive_dir, "data")
    if not os.path.exists(full_archive_subdir):
        os.makedirs(full_archive_subdir)
    if not os.path.exists(full_archive_datadir):
        os.makedirs(full_archive_datadir)

    # Copy the file to the archive
    qualified_filename = os.path.join(workdir, product_filename)
    qualified_archive_product_filename = os.path.join(full_archive_dir, product_filename)
    copyfile(qualified_filename, qualified_archive_product_filename)

    # Copy any files it points to to the archive as well
    try:
        p = read_xml_product(qualified_filename)

        # Remove all files this points to
        if hasattr(p, "get_all_filenames"):
            data_filenames = p.get_all_filenames()
            for data_filename in data_filenames:
                if data_filename is not None and data_filename != "default_filename.fits" and data_filename != "":
                    qualified_data_filename = os.path.join(workdir, data_filename)
                    qualified_archive_data_filename = os.path.join(full_archive_dir, data_filename)
                    copyfile(qualified_data_filename, qualified_archive_data_filename)

        else:
            logger.warn("Product " + qualified_filename + " has no 'get_all_filenames' method.")

    except Exception as e:
        logger.warn("Failsafe exception block triggered when trying to save statistics product in archive. " +
                    "Exception was: " + str(e))

    return


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
