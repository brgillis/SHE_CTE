""" @file measure_bias.py

    Created 7 Apr 2017

    Primary execution loop for measuring bias in shear estimates.
"""

__updated__ = "2019-07-08"

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

from _pickle import UnpicklingError
import os

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_BiasMeasurement.find_files import recursive_find_files
from SHE_PPT import products
from SHE_PPT.file_io import read_listfile, read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.math import combine_linregress_statistics, BiasMeasurements, combine_bfd_sum_statistics
from SHE_PPT.pipeline_utility import archive_product, read_config, ConfigKeys
import multiprocessing as mp
import numpy as np


bootstrap_threshold = 2

logger = getLogger(__name__)


class MethodStatistics(object):
    """Class to contain lists of g1 and g2 bias statistics.
    """

    def __init__(self):
        self.g1_statistics = None
        self.g2_statistics = None
        self.bfd_statistics = None


class MethodStatisticsList(object):
    """Class to contain lists of g1 and g2 bias statistics.
    """

    def __init__(self):
        self.g1_statistics_list = []
        self.g2_statistics_list = []
        self.bfd_statistics_list = []


def read_statistics(shear_statistics_file, workdir):

    # This is a merge point, so we get the file as a tuple of length 1 in the listfile
    if isinstance(shear_statistics_file, tuple) or isinstance(shear_statistics_file, list):
        if len(shear_statistics_file) == 1:
            shear_statistics_file = shear_statistics_file[0]
        else:
            raise ValueError("Unexpected format of shear statistics listfile.")

    try:
        shear_statistics_prod = read_xml_product(os.path.join(workdir, shear_statistics_file))
    except (EOFError, UnpicklingError) as e:
        logger.warn("File " + os.path.join(workdir, shear_statistics_file) + " seems to be corrupted: " + str(e))
        return

    all_method_statistics = {}

    for method in mv.estimation_methods:

        method_statistics = MethodStatistics()
        method_shear_statistics = shear_statistics_prod.get_method_statistics(method)

        if not method == "BFD":  # get info for method if not BFD

            method_statistics.g1_statistics = method_shear_statistics[0]
            method_statistics.g2_statistics = method_shear_statistics[1]

        elif method_shear_statistics is not None:  # get info for BFD method
            method_statistics.bfd_statistics = method_shear_statistics

        all_method_statistics[method] = method_statistics

    return all_method_statistics


def measure_bias_from_args(args):
    """
    @brief
        Perform bias measurement, given arguments from the command-line.

    @param args

    @return None
    """
    logger.debug("Entering measure_bias_from_args.")

    # Read in the pipeline config
    try:
        pipeline_config = read_config(args.pipeline_config, workdir=args.workdir)
        if pipeline_config is None:
            pipeline_config = {}
    except Exception as e:
        logger.warn("Failsafe exception block triggered when trying to read pipeline config. " +
                    "Exception was: " + str(e))
        pipeline_config = {}

    if args.number_threads is not None:
        number_threads = args.number_threads
    elif ConfigKeys.MB_NUM_THREADS.value in pipeline_config:
        number_threads = pipeline_config[ConfigKeys.MB_NUM_THREADS.value]
        if number_threads.lower() == "none":
            number_threads = 1
    else:
        number_threads = 1
    elif number_threads < 1:
        number_threads = 1

    # If number_threads is 0 or lower, assume it means this many fewer than the cpu count
    if number_threads <= 0:
        number_threads = max(1, mp.cpu_count() + number_threads)

    if args.archive_dir is not None:
        archive_dir = args.archive_dir
    elif ConfigKeys.MB_ARCHIVE_DIR.value in pipeline_config:
        archive_dir = pipeline_config[ConfigKeys.MB_ARCHIVE_DIR.value]
        if archive_dir.lower() == "none":
            archive_dir = None
    else:
        archive_dir = None

    if args.webdav_dir is not None:
        webdav_dir = args.webdav_dir
    elif ConfigKeys.MB_WEBDAV_DIR.value in pipeline_config:
        webdav_dir = pipeline_config[ConfigKeys.MB_WEBDAV_DIR.value]
        if webdav_dir.lower() == "none":
            webdav_dir = None
    else:
        webdav_dir = None

    if args.webdav_archive is not None:
        webdav_archive = args.webdav_archive
    elif ConfigKeys.MB_WEBDAV_ARCHIVE.value in pipeline_config:
        webdav_archive = pipeline_config[ConfigKeys.MB_WEBDAV_ARCHIVE.value]
        if webdav_archive.lower() == "none":
            webdav_archive = None
    else:
        webdav_archive = None

    # Get a list of input files

    if args.shear_bias_statistics is None or args.shear_bias_statistics == "None":
        # Working in recovery mode, so search within the workdir to find the files
        shear_statistics_files = recursive_find_files(args.workdir, number_threads=number_threads)
    else:
        shear_statistics_files = read_listfile(os.path.join(args.workdir, args.shear_bias_statistics))

    # Load in statistics from each file
    l_method_shear_statistics = []

    # Read simply if number_threads==1
    if number_threads == 1:

        l_method_shear_statistics = [read_statistics(
            shear_statistics_file, workdir=args.workdir) for shear_statistics_file in shear_statistics_files]

    # Otherwise use multiprocessing
    else:

        pool = mp.Pool(processes=number_threads)
        l_method_shear_statistics = [pool.apply(read_statistics, args=(
            shear_statistics_file, args.workdir)) for shear_statistics_file in shear_statistics_files]

    # Combine the statistics into a single list for each method
    method_shear_statistics_lists = {}
    for method in mv.estimation_methods:

        method_shear_statistics_list = MethodStatisticsList()

        method_shear_statistics_list.g1_statistics_list = [l_method_shear_statistics[i]
                                                           [method].g1_statistics for i in range(len(l_method_shear_statistics))]
        method_shear_statistics_list.g2_statistics_list = [l_method_shear_statistics[i]
                                                           [method].g2_statistics for i in range(len(l_method_shear_statistics))]
        method_shear_statistics_list.bfd_statistics_list = [l_method_shear_statistics[i]
                                                            [method].bfd_statistics for i in range(len(l_method_shear_statistics))]

        method_shear_statistics_lists[method] = method_shear_statistics_list

    # Calculate the bias and compile into a data product
    bias_measurement_prod = products.shear_bias_measurements.create_shear_bias_measurements_product()

    for method in mv.estimation_methods:

        if not method == "BFD":
            # do bias measurement for all methods but BFD
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

        else:
            # do bias measurement for BFD
            if len(method_shear_statistics_lists[method].bfd_statistics_list) > 0:
                try:
                    g1_bias_measurements = BiasMeasurements(combine_bfd_sum_statistics(
                        method_shear_statistics_lists[method].bfd_statistics_list, do_g1=True))
                    g2_bias_measurements = BiasMeasurements(combine_bfd_sum_statistics(
                        method_shear_statistics_lists[method].bfd_statistics_list, do_g1=False))
                except np.linalg.linalg.LinAlgError as e:
                    logger.warn("Unable to calculate bias for BFD. Exception was: " + str(e))
                    g1_bias_measurements = None
                    g2_bias_measurements = None
            else:
                g1_bias_measurements = None
                g2_bias_measurements = None

        bias_measurement_prod.set_method_bias_measurements(method, g1_bias_measurements, g2_bias_measurements)

    write_xml_product(bias_measurement_prod, args.shear_bias_measurements, workdir=args.workdir)

    # Try to archive the product

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
