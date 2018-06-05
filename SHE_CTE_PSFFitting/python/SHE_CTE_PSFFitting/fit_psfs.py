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

from astropy.io import fits
from astropy.table import Table

from SHE_CTE_PSFFitting import magic_values as mv
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT import products
from SHE_PPT.file_io import (read_listfile, write_listfile,
                             read_pickled_product, write_pickled_product,
                             get_allowed_filename, find_file_in_path)
from SHE_PPT.logging import getLogger
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.table_utility import is_in_format, table_to_hdu
import numpy as np


test_mode = True

def fit_psfs(args, dry_run = False):
    """
        Mock run of PSF Fitting.
    """

    logger = getLogger(mv.logger_name)

    # Load in the files in turn to make sure there aren't any issues with them.

    # Data images

    if dry_run:
        dry_label = "_dry"
    else:
        dry_label = ""

    logger.info("Reading mock" + dry_label + " data images and detections tables...")

    frame_stack = SHEFrameStack.read(exposure_listfile_filename = args.data_images,
                                     detections_listfile_filename = args.detections_tables,
                                     workdir = args.workdir,
                                     clean_detections = True,
                                     apply_sc3_fix = True)

    # AocsTimeSeries products

    if args.aocs_time_series_products is not None:

        logger.info("Reading mock" + dry_label + " aocs time series products...")

        aocs_time_series_product_filenames = read_listfile(join(args.workdir, args.aocs_time_series_products))
        aocs_time_series_products = []

        for i, filename in enumerate(aocs_time_series_product_filenames):
            aocs_time_series_products.append(read_pickled_product(join(args.workdir, filename)))
            if not isinstance(aocs_time_series_products[i], products.aocs_time_series.DpdSheAocsTimeSeriesProduct):
                raise ValueError("AocsTimeSeries product from " + filename + " is invalid type.")

    else:

        aocs_time_series_products = None

    # PSFCalibration products

    if args.psf_calibration_product is not None:

        logger.info("Reading mock" + dry_label + " PSF calibration products...")

        psf_calibration_product = read_pickled_product(join(args.workdir, args.psf_calibration_product))

        if not isinstance(psf_calibration_product, products.psf_calibration.DpdShePSFCalibrationProduct):
            raise ValueError("PSFCalibration product from " + filename + " is invalid type.")

    else:

        psf_calibration_product = None

    # Set up mock output in the correct format

    logger.info("Outputting mock" + dry_label + " PSF Field Params...")

    field_param_filenames = []

    for i in range(len(frame_stack.exposures)):

        field_param_filename = get_allowed_filename("PSF_FieldParam", str(i))
        field_param_filenames.append(field_param_filename)

        field_param_product = products.psf_field_params.create_dpd_she_psf_field_params()

        write_pickled_product(field_param_product, field_param_filename)

    write_listfile(args.psf_field_params, join(args.workdir, field_param_filenames))

    return


