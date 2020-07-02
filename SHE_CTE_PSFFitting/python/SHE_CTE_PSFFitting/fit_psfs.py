""" @file fit_psf.py

    Created 12 Oct 2017

    Function for performing a dry run of mock psf fitting.
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

from copy import deepcopy
import os
from os.path import join

from SHE_PPT import magic_values as ppt_mv
from SHE_PPT import products
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, find_file_in_path)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import get_conditional_product
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from astropy.io import fits
from astropy.table import Table

import SHE_CTE
from SHE_CTE_PSFFitting import magic_values as mv
from SHE_PPT.table_formats.she_psf_tm_state import initialise_psf_field_tm_state_table
from SHE_PPT.table_formats.she_simulated_catalog import tf as simc_tf
import numpy as np

products.psf_field_params.init()

test_mode = True


def fit_psfs(args, dry_run=False):
    """
        Mock run of PSF Fitting.
    """

    logger = getLogger(__name__)

    # Load in the files in turn to make sure there aren't any issues with them.

    # Data images

    if dry_run:
        dry_label = "_dry"
    else:
        dry_label = ""

    logger.info("Reading mock" + dry_label + " data images and she_simulated_catalog files...")

    frame_stack = SHEFrameStack.read(exposure_listfile_filename=args.data_images,
                                     seg_listfile_filename=args.segmentation_images,
                                     detections_listfile_filename=args.she_simulated_catalog_listfile,
                                     workdir=args.workdir,
                                     clean_detections=True,
                                     apply_sc3_fix=True,
                                     memmap=True,
                                     mode='denywrite')

    # AocsTimeSeries products

    if args.aocs_time_series_products is not None:

        logger.info("Reading mock" + dry_label + " aocs time series products...")

        aocs_time_series_product_filenames = read_listfile(join(args.workdir, args.aocs_time_series_products))
        aocs_time_series_products = []

        for i, filename in enumerate(aocs_time_series_product_filenames):
            aocs_time_series_products.append(read_pickled_product(join(args.workdir, filename)))
            if not isinstance(aocs_time_series_products[i], products.le1_aocs_time_series.DpdSheAocsTimeSeriesProduct):
                raise ValueError("AocsTimeSeries product from " + filename + " is invalid type.")

    else:

        aocs_time_series_products = None

    # PSFCalibration product

    psf_calibration_product = get_conditional_product(args.psf_calibration_product, args.workdir)

    if psf_calibration_product is not None and not isinstance(psf_calibration_product,
                                                              products.she_psf_calibration_parameters.dpdShePsfCalibrationParameters):
        raise ValueError("PSFCalibration product from " + filename + " is invalid type.")

    # Set up mock output in the correct format

    logger.info("Outputting mock" + dry_label + " PSF Field Params...")

    field_param_product_filenames = []

    for i in range(len(frame_stack.exposures)):

        # Get a filename for the product
        field_param_product_filename = get_allowed_filename("PSF-FieldParam", str(i), extension=".xml",
                                                            version=SHE_CTE.__version__)
        field_param_product_filenames.append(field_param_product_filename)

        # Get a filename for the table and write it out
        field_param_table_filename = get_allowed_filename("PSF-TelescopeModel", str(i), extension=".fits",
                                                          version=SHE_CTE.__version__)
        field_param_table = initialise_psf_field_tm_state_table()
        field_param_table.add_row({})  # Add a row of zeros

        qualified_field_param_table_filename = join(args.workdir, field_param_table_filename)

        field_param_table.write(join(args.workdir, field_param_table_filename))

        logger.info("Wrote field params table to " + qualified_field_param_table_filename)

        # Create and write the data product
        field_param_product = products.she_psf_field_parameters.create_dpd_she_psf_field_parameters()
        field_param_product.set_zernike_mode_filename(field_param_table_filename)

        qualified_field_param_product_filename = join(args.workdir, field_param_product_filename)
        write_pickled_product(field_param_product, qualified_field_param_product_filename)

        logger.info("Wrote field params product to " + qualified_field_param_product_filename)

    # Write out a listfile of the products
    qualified_listfile_filename = join(args.workdir, args.psf_field_params)
    write_listfile(qualified_listfile_filename, field_param_product_filenames)

    logger.info("Wrote field params listfile to " + qualified_listfile_filename)

    return
