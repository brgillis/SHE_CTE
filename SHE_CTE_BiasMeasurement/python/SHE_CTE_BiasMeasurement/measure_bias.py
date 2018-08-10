""" @file measure_bias.py

    Created 7 Apr 2017

    Primary execution loop for measuring bias in shear estimates.
"""

__updated__ = "2018-08-10"

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
from SHE_PPT.file_io import read_listfile, read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.math import combine_linregress_statistics, BiasMeasurements

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_PPT.pipeline_utility import archive_product, read_config
import numpy as np


bootstrap_threshold = 50

archive_dir_key = "SHE_CTE_MeasureBias_archive_dir"
webdav_dir_key = "SHE_CTE_MeasureBias_webdav_dir"
webdav_archive_key = "SHE_CTE_MeasureBias_webdav_archive"


class MethodStatisticsList(object):
    """Class to contain lists of g1 and g2 bias statistics.
    """

    def __init__(self):
        self.g1_statistics_list = []
        self.g2_statistics_list = []


def measure_bias_from_args(args):
    """
    @brief
        Perform bias measurement, given arguments from the command-line.

    @param args

    @return None
    """

    logger = getLogger(mv.logger_name)
    logger.debug("Entering measure_bias_from_args.")

    # Get a list of input files
    shear_statistics_files = read_listfile(os.path.join(args.workdir, args.shear_bias_statistics))

    # Load in statistics from each file
    method_shear_statistics_lists = {}
    for method in mv.estimation_methods:
        method_shear_statistics_lists[method] = MethodStatisticsList()

    for shear_statistics_file in shear_statistics_files:

        # This is a merge point, so we get the file as a tuple of length 1 in the listfile
        if isinstance(shear_statistics_file, tuple) or isinstance(shear_statistics_file, list):
            if len(shear_statistics_file) == 1:
                shear_statistics_file = shear_statistics_file[0]
            else:
                raise ValueError("Unexpected format of shear statistics listfile.")

        shear_statistics_prod = read_xml_product(os.path.join(args.workdir, shear_statistics_file))

        for method in mv.estimation_methods:
            method_shear_statistics = shear_statistics_prod.get_method_statistics(method)
            if method_shear_statistics[0] is not None:
                method_shear_statistics_lists[method].g1_statistics_list.append(method_shear_statistics[0])
            if method_shear_statistics[1] is not None:
                method_shear_statistics_lists[method].g2_statistics_list.append(method_shear_statistics[1])

    # Calculate the bias and compile into a data product
    bias_measurement_prod = products.shear_bias_measurements.create_shear_bias_measurements_product()

    for method in mv.estimation_methods:

        if len(method_shear_statistics_lists[method].g1_statistics_list) >= bootstrap_threshold:
            # We have enough data to calculate bootstrap errors
            g1_bias_measurements = calculate_bootstrap_bias_measurements(
                method_shear_statistics_lists[method].g1_statistics_list, seed=args.bootstrap_seed)

        elif len(method_shear_statistics_lists[method].g1_statistics_list) > 0:
            # Not enough for bootstrap errors - calculate simply
            g1_bias_measurements = BiasMeasurements(
                combine_linregress_statistics(method_shear_statistics_lists[method].g1_statistics_list))
        else:
            g1_bias_measurements = None

        if len(method_shear_statistics_lists[method].g2_statistics_list) >= bootstrap_threshold:
            # We have enough data to calculate bootstrap errors
            g2_bias_measurements = calculate_bootstrap_bias_measurements(
                method_shear_statistics_lists[method].g2_statistics_list, seed=args.bootstrap_seed)

        elif len(method_shear_statistics_lists[method].g2_statistics_list) > 0:
            # Not enough for bootstrap errors - calculate simply
            g2_bias_measurements = BiasMeasurements(
                combine_linregress_statistics(method_shear_statistics_lists[method].g2_statistics_list))
        else:
            g2_bias_measurements = None

        bias_measurement_prod.set_method_bias_measurements(method, g1_bias_measurements, g2_bias_measurements)

    write_xml_product(bias_measurement_prod, os.path.join(args.workdir, args.shear_bias_measurements))

    # Try to archive the product

    # First get the pipeline config so we can figure out where to archive it
    try:
        pipeline_config = read_config(args.pipeline_config, workdir=args.workdir)
    except Exception as e:
        logger.warn("Failsafe exception block triggered when trying to read pipeline config. " +
                    "Exception was: " + str(e))
        pipeline_config = {}

    if archive_dir_key in pipeline_config:
        archive_dir = pipeline_config[archive_dir_key]
    else:
        archive_dir = args.archive_dir
        
    if webdav_dir_key in pipeline_config:
        webdav_dir = pipeline_config[webdav_dir_key]
    else:
        webdav_dir = args.webdav_dir

    if webdav_archive_key in pipeline_config:
        webdav_archive = pipeline_config[webdav_archive_key].lower()=="true"
    else:
        webdav_archive = args.webdav_archive

    # If we're archiving with webdav, determine its mount dir and the full archive directory
    if webdav_archive and archive_dir is not None:
        full_archive_dir = os.path.join(webdav_dir, archive_dir)
    else:
        full_archive_dir = archive_dir

    if archive_dir is not None:
        try:
            archive_product(product_filename=args.shear_bias_measurements,
                            archive_dir=full_archive_dir,
                            workdir=args.workdir)
        except Exception as e:
            logger.warn("Failsafe exception block triggered when trying to save bias product in archive. " +
                        "Exception was: " + str(e))

    logger.debug("Exiting measure_bias_from_args.")

    return


def calculate_bootstrap_bias_measurements(statistics_list, n_bootstrap=1000, seed=0):
    """Calculates a BiasMeasurements object using bootstrap errors from a list of statistics objects.
    """

    # Seed the random number generator
    np.random.seed(seed)

    # Get a base object for the m and c calculations
    bias_measurements = BiasMeasurements(combine_linregress_statistics(statistics_list))

    # Bootstrap to get errors on m and c
    n_sample = len(statistics_list)

    m_bs = np.empty(n_bootstrap)
    c_bs = np.empty(n_bootstrap)
    statistics_array = np.array(statistics_list)
    for i in range(n_bootstrap):
        u = np.random.random_integers(0, n_sample - 1, n_sample)
        bias_measurements_bs = BiasMeasurements(combine_linregress_statistics(statistics_array[u]))
        m_bs[i] = bias_measurements_bs.m
        c_bs[i] = bias_measurements_bs.c

    # Update the bias measurements in the output object
    bias_measurements.m_err = np.std(m_bs)
    bias_measurements.c_err = np.std(c_bs)

    return bias_measurements
