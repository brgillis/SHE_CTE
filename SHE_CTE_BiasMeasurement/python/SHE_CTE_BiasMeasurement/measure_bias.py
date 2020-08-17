""" @file measure_bias.py

    Created 7 Apr 2017

    Primary execution loop for measuring bias in shear estimates.
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

from _pickle import UnpicklingError
import os

from SHE_PPT import products
from SHE_PPT.file_io import read_listfile, read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.math import (combine_linregress_statistics, BiasMeasurements, combine_bfd_sum_statistics,
                          LinregressStatistics, BFDSumStatistics)
from SHE_PPT.pipeline_utility import archive_product, read_config, ConfigKeys

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_BiasMeasurement.find_files import recursive_find_files
from SHE_CTE_BiasMeasurement.print_bias import print_bias_from_product
from SHE_PPT.products.shear_bias_statistics import create_dpd_shear_bias_statistics_from_stats
import multiprocessing as mp
import numpy as np

__updated__ = "2020-08-17"

bootstrap_threshold = 2
default_number_threads = 8

logger = getLogger(__name__)


class MethodStatistics(object):
    """Class to contain lists of g1 and g2 bias statistics.
    """

    def __init__(self):
        self.g1_statistics = None
        self.g2_statistics = None
        self.bfd_statistics = None


class MethodMeasurements(object):
    """Class to contain lists of g1 and g2 bias measurements.
    """

    def __init__(self):
        self.g1_measurements = None
        self.g2_measurements = None


class MethodStatisticsList(object):
    """Class to contain lists of g1 and g2 bias statistics.
    """

    def __init__(self):
        self.g1_statistics_list = []
        self.g2_statistics_list = []
        self.bfd_statistics_list = []


class MethodMeasurementsList(object):
    """Class to contain lists of g1 and g2 bias statistics.
    """

    def __init__(self):
        self.g1_measurements_list = []
        self.g2_measurements_list = []


def read_statistics(shear_statistics_file, workdir, recovery_mode, use_bias_only=False):

    # In recovery mode, adjust the work directory to match where the data will be for each file
    if recovery_mode:
        extra_workdir_portion = os.path.split(shear_statistics_file)[0]

        # Check for if the product was accidentally placed in the data directory
        if len(extra_workdir_portion) >= 4 and extra_workdir_portion[-4:] == "data":
            extra_workdir_portion = extra_workdir_portion[:-5]

        workdir = os.path.join(workdir, extra_workdir_portion)

        shear_statistics_file = shear_statistics_file.replace(extra_workdir_portion, '', 1)
        if shear_statistics_file[0] == '/':
            shear_statistics_file = shear_statistics_file[1:]

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
        return None, None

    all_method_statistics = {}
    all_method_measurements = {}

    for method in mv.estimation_methods:

        if not use_bias_only:

            method_statistics = MethodStatistics()
            try:
                method_shear_statistics = shear_statistics_prod.get_method_bias_statistics(method, workdir=workdir)
            except Exception as e:
                logger.warning("Shear bias statistics in product " + os.path.join(workdir, shear_statistics_file) +
                               " for method " + method + " appear to be corrupted. Exception was: " + str(e))
                return None, None

            if not method == "BFD":  # get info for method if not BFD

                method_statistics.g1_statistics = method_shear_statistics[0]
                method_statistics.g2_statistics = method_shear_statistics[1]

            elif method_shear_statistics is not None:  # get info for BFD method
                method_statistics.bfd_statistics = method_shear_statistics

            all_method_statistics[method] = method_statistics

        else:

            all_method_statistics[method] = []

        method_measurements = MethodMeasurements()
        method_bias_measurements = shear_statistics_prod.get_method_bias_measurements(method, workdir=workdir)

        if method_bias_measurements is None:
            method_measurements.g1_measurements = None
            method_measurements.g2_measurements = None
        else:
            method_measurements.g1_measurements = method_bias_measurements[0]
            method_measurements.g2_measurements = method_bias_measurements[1]

        all_method_measurements[method] = method_measurements

    return all_method_statistics, all_method_measurements


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
            number_threads = default_number_threads
    else:
        number_threads = default_number_threads

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
        logger.info("No list of statistics files supplied, so searching for files which match the filename '" +
                    "shear_bias_statistics.xml" + "'.")
        shear_statistics_files = recursive_find_files(base_dir=args.workdir,
                                                      bias_statistics_filename=args.recovery_bias_statistics_filename,
                                                      bias_measurements_filename=args.recovery_bias_measurements_filename,
                                                      number_threads=number_threads)
        recovery_mode = True
    else:
        shear_statistics_files = read_listfile(os.path.join(args.workdir, args.shear_bias_statistics))
        recovery_mode = False

    # Load in statistics from each file

    # Read simply if number_threads==1
    if number_threads == 1:

        lt_method_shear_statistics = [read_statistics(
            shear_statistics_file, workdir=args.workdir,
            recovery_mode=recovery_mode, use_bias_only=args.use_bias_only) for shear_statistics_file in shear_statistics_files]

    # Otherwise use multiprocessing
    else:

        pool = mp.Pool(processes=number_threads)
        lt_method_shear_statistics = [pool.apply(read_statistics, args=(
            shear_statistics_file, args.workdir, recovery_mode, args.use_bias_only)) for shear_statistics_file in shear_statistics_files]

    l_method_shear_statistics, l_method_bias_measurements = zip(*lt_method_shear_statistics)

    # Combine the statistics into a single list for each method
    method_shear_statistics_lists = {}
    method_bias_measurements_lists = {}

    def remove_values_from_list(the_list, val):
        return [value for value in the_list if value != val]

    l_method_shear_statistics = remove_values_from_list(l_method_shear_statistics, None)
    l_method_bias_measurements = remove_values_from_list(l_method_bias_measurements, None)

    have_some_data = False
    missing_shear_statistics = False

    if not args.use_bias_only:

        for method in mv.estimation_methods:

            method_shear_statistics_list = MethodStatisticsList()

            l_g1_statistics_list = remove_values_from_list([l_method_shear_statistics[i]
                                                            [method].g1_statistics for i in range(len(l_method_shear_statistics))], None)
            l_g2_statistics_list = remove_values_from_list([l_method_shear_statistics[i]
                                                            [method].g2_statistics for i in range(len(l_method_shear_statistics))], None)
            l_bfd_statistics_list = remove_values_from_list([l_method_shear_statistics[i]
                                                             [method].bfd_statistics for i in range(len(l_method_shear_statistics))], None)

            method_shear_statistics_list.g1_statistics_list = []
            method_shear_statistics_list.g2_statistics_list = []
            method_shear_statistics_list.bfd_statistics_list = []

            # Compress the lists to be 1D
            for uncompressed_list, final_list in ((l_g1_statistics_list, method_shear_statistics_list.g1_statistics_list),
                                                  (l_g2_statistics_list, method_shear_statistics_list.g2_statistics_list),
                                                  (l_bfd_statistics_list, method_shear_statistics_list.bfd_statistics_list),):
                for item in uncompressed_list:
                    if isinstance(item, list):
                        final_list += item
                        if len(item) == 0:
                            missing_shear_statistics = True
                    elif isinstance(item, LinregressStatistics) or isinstance(item, BFDSumStatistics):
                        final_list.append(item)
                    else:
                        raise ValueError("Unexpected type of bias statistics: " + str(type(item)))

            method_shear_statistics_lists[method] = method_shear_statistics_list

            if (len(method_shear_statistics_list.g1_statistics_list) > 0 or
                len(method_shear_statistics_list.g2_statistics_list) > 0 or
                    len(method_shear_statistics_list.bfd_statistics_list) > 0):
                have_some_data = True

    if missing_shear_statistics or args.use_bias_only:

        for method in mv.estimation_methods:

            # Check if we have some bias measurements

            def remove_missing_bms_from_list(the_list):
                new_list = []
                for value in the_list:
                    try:
                        if (value.m != '' and value.m_err != '' and value.m_err > 0
                            and value.c != '' and value.c_err != '' and value.c_err > 0):
                            new_list.append(value)
                    except (TypeError, ValueError) as _:
                        pass
                return new_list

            method_bias_measurements_list = MethodMeasurementsList()

            method_bias_measurements_list.g1_measurements_list = remove_missing_bms_from_list(remove_values_from_list([l_method_bias_measurements[i]
                                                            [method].g1_measurements for i in range(len(l_method_bias_measurements))], None))
            method_bias_measurements_list.g2_measurements_list = remove_missing_bms_from_list(remove_values_from_list([l_method_bias_measurements[i]
                                                            [method].g2_measurements for i in range(len(l_method_bias_measurements))], None))

            method_bias_measurements_lists[method] = method_bias_measurements_list

            if (len(method_bias_measurements_list.g1_measurements_list) > 0 or
                len(method_bias_measurements_list.g2_measurements_list) > 0):
                have_some_data = True

    if not have_some_data:
        raise RuntimeError("No shear bias statistics or measurements are available; aborting.")
    elif missing_shear_statistics and not args.use_bias_only:
        logger.warn("Some shear statistics are missing; relying on shear bias measurements only.")

    # Calculate the bias and compile into a data product
    if args.store_measurements_only or missing_shear_statistics or args.use_bias_only:
        bias_measurement_prod = create_dpd_shear_bias_statistics_from_stats(BFD_bias_statistics=[],
                                                                            KSB_bias_statistics=([],
                                                                                                 []),
                                                                            LensMC_bias_statistics=([],
                                                                                                    []),
                                                                            MomentsML_bias_statistics=([],
                                                                                                       []),
                                                                            REGAUSS_bias_statistics=([],
                                                                                                     []),
                                                                            workdir=args.workdir)
    else:
        bias_measurement_prod = create_dpd_shear_bias_statistics_from_stats(BFD_bias_statistics=method_shear_statistics_lists["BFD"].bfd_statistics_list,
                                                                            KSB_bias_statistics=(method_shear_statistics_lists["KSB"].g1_statistics_list,
                                                                                                 method_shear_statistics_lists["KSB"].g2_statistics_list),
                                                                            LensMC_bias_statistics=(method_shear_statistics_lists["LensMC"].g1_statistics_list,
                                                                                                    method_shear_statistics_lists["LensMC"].g2_statistics_list),
                                                                            MomentsML_bias_statistics=(method_shear_statistics_lists["MomentsML"].g1_statistics_list,
                                                                                                       method_shear_statistics_lists["MomentsML"].g2_statistics_list),
                                                                            REGAUSS_bias_statistics=(method_shear_statistics_lists["REGAUSS"].g1_statistics_list,
                                                                                                     method_shear_statistics_lists["REGAUSS"].g2_statistics_list),
                                                                            workdir=args.workdir)

    for method in mv.estimation_methods:

        if missing_shear_statistics or args.use_bias_only:

            # Calculate from bias measurements

            if len(method_bias_measurements_lists[method].g1_measurements_list) >= bootstrap_threshold:

                logger.info("Calculating bootstrap bias measurements for method " + method + " from list of bias measurements.")

                # We have enough data to calculate bootstrap errors
                g1_bias_measurements = calculate_bootstrap_bias_measurements_from_measurements(
                    method_bias_measurements_lists[method].g1_measurements_list, seed=args.bootstrap_seed)

            elif len(method_bias_measurements_lists[method].g1_measurements_list) >= 0:

                logger.info("Calculating non-bootstrap bias measurements for method " + method + " from list of bias measurements.")

                # Not enough for bootstrap errors - calculate simply
                g1_bias_measurements = combine_bias_measurements(
                    method_bias_measurements_lists[method].g1_measurements_list)

            else:
                g1_bias_measurements = None

            if len(method_bias_measurements_lists[method].g2_measurements_list) >= bootstrap_threshold:
                # We have enough data to calculate bootstrap errors
                g2_bias_measurements = calculate_bootstrap_bias_measurements_from_measurements(
                    method_bias_measurements_lists[method].g2_measurements_list, seed=args.bootstrap_seed)
            elif len(method_bias_measurements_lists[method].g2_measurements_list) >= 0:
                # Not enough for bootstrap errors - calculate simply
                g2_bias_measurements = combine_bias_measurements(
                    method_bias_measurements_lists[method].g2_measurements_list)
            else:
                g2_bias_measurements = None

        elif not method == "BFD":
            # do bias measurement for all methods but BFD
            if len(method_shear_statistics_lists[method].g1_statistics_list) >= bootstrap_threshold:

                logger.info("Calculating bootstrap bias measurements for method " + method + " from list of bias statistics.")

                # We have enough data to calculate bootstrap errors
                g1_bias_measurements = calculate_bootstrap_bias_measurements_from_statistics(
                    method_shear_statistics_lists[method].g1_statistics_list, seed=args.bootstrap_seed)

            elif len(method_shear_statistics_lists[method].g1_statistics_list) > 0:

                logger.info("Calculating non-bootstrap bias measurements for method " + method + " from list of bias statistics.")

                # Not enough for bootstrap errors - calculate simply
                g1_bias_measurements = BiasMeasurements(
                    combine_linregress_statistics(method_shear_statistics_lists[method].g1_statistics_list))
            else:
                g1_bias_measurements = None

            if len(method_shear_statistics_lists[method].g2_statistics_list) >= bootstrap_threshold:
                # We have enough data to calculate bootstrap errors
                g2_bias_measurements = calculate_bootstrap_bias_measurements_from_statistics(
                    method_shear_statistics_lists[method].g2_statistics_list, seed=args.bootstrap_seed)

            elif len(method_shear_statistics_lists[method].g2_statistics_list) > 0:
                # Not enough for bootstrap errors - calculate simply
                g2_bias_measurements = BiasMeasurements(
                    combine_linregress_statistics(method_shear_statistics_lists[method].g2_statistics_list))
            else:
                g2_bias_measurements = None

        else:
            # do bias measurement for BFD
            try:
                if len(method_shear_statistics_lists[method].bfd_statistics_list) >= bootstrap_threshold:

                    logger.info("Calculating bootstrap bias measurements for method " + method + " from list of bias statistics.")

                    # We have enough data to calculate bootstrap errors
                    g1_bias_measurements, g2_bias_measurements = calculate_bfd_bootstrap_bias_measurements_from_statistics(
                        method_shear_statistics_lists[method].bfd_statistics_list, seed=args.bootstrap_seed)

                elif len(method_shear_statistics_lists[method].bfd_statistics_list) > 0:

                    logger.info("Calculating bootstrap bias measurements for method " + method + " from list of bias statistics.")

                    g1_bias_measurements = BiasMeasurements(combine_bfd_sum_statistics(
                        method_shear_statistics_lists[method].bfd_statistics_list, do_g1=True))
                    g2_bias_measurements = BiasMeasurements(combine_bfd_sum_statistics(
                        method_shear_statistics_lists[method].bfd_statistics_list, do_g1=False))

                else:

                    g1_bias_measurements = None
                    g2_bias_measurements = None

            except np.linalg.linalg.LinAlgError as e:
                logger.warn("Unable to calculate bias for BFD. Exception was: " + str(e))
                g1_bias_measurements = None
                g2_bias_measurements = None

        bias_measurement_prod.set_method_bias_measurements(
            method, (g1_bias_measurements, g2_bias_measurements), workdir=args.workdir)

    # Print the bias measurements
    print_bias_from_product(bias_measurement_prod, workdir=args.workdir)

    logger.info("Writing combined bias measurments to " + os.path.join(args.workdir, args.shear_bias_measurements))
    write_xml_product(bias_measurement_prod, args.shear_bias_measurements, workdir=args.workdir)

    # Try to archive the product

    # If we're archiving with webdav, determine its mount dir and the full archive directory
    if webdav_archive and archive_dir is not None:
        full_archive_dir = os.path.join(webdav_dir, archive_dir)
    else:
        full_archive_dir = archive_dir

    if archive_dir is not None:
        try:
            logger.info("Archiving combined bias measurments to " + os.path.join(full_archive_dir, args.workdir, args.shear_bias_measurements))
            archive_product(product_filename=args.shear_bias_measurements,
                            archive_dir=full_archive_dir,
                            workdir=args.workdir)
        except Exception as e:
            logger.warn("Failsafe exception block triggered when trying to save bias product in archive. " +
                        "Exception was: " + str(e))

    logger.debug("Exiting measure_bias_from_args.")

    return


def combine_bias_statistics(l_bias_statistics):
    """ Simple passthrough function to combine stats into a BiasMeasurements object.
    """

    return BiasMeasurements(combine_linregress_statistics(l_bias_statistics))


def combine_bias_measurements(l_bias_measurements):
    """Function to combine a list of bias measurements into a single value, trusting the errors and weights are all correct
    and there aren't any correlations between errors and bias values.
    
    Parameters
    ----------
    l_bias_measurements : iterable<measure_bias.MethodMeasurements>
        The list of bias measurements to combine
    
    Return
    ------
    bias_measurements : <BiasMeasurements>
    
    """

    num_bias_measurements = len(l_bias_measurements)

    lm = np.array([l_bias_measurements[i].m for i in range(num_bias_measurements)])
    lm_err = np.array([l_bias_measurements[i].m_err for i in range(num_bias_measurements)])
    lc = np.array([l_bias_measurements[i].c for i in range(num_bias_measurements)])
    lc_err = np.array([l_bias_measurements[i].c_err for i in range(num_bias_measurements)])
    lmc_covar = np.array([l_bias_measurements[i].mc_covar for i in range(num_bias_measurements)])

    # Calculate weight arrays
    lwm = lm_err ** (-2)
    lwc = lc_err ** (-2)
    try:
        lwmc = np.where(np.logical_or(np.isinf(lmc_covar), np.isnan(lmc_covar)), 0, lmc_covar ** (-1))
        bad_covar = False
    except TypeError as e:
        # If we have bad covars, get covar weighting from regular vars
        bad_covar = True

    # Determine a mask for both m and c
    m_mask = np.logical_or(np.logical_or(np.isnan(lwm), np.isinf(lwm)), lwm <= 0)
    c_mask = np.logical_or(np.logical_or(np.isnan(lwc), np.isinf(lwc)), lwc <= 0)
    mc_mask = np.logical_or(m_mask, c_mask)

    # Calculate total weights
    wm = lwm[~mc_mask].sum()
    wc = lwc[~mc_mask].sum()
    if not bad_covar:
        wmc = lwmc[~mc_mask].sum()

    # Check for zero weight, otherwise calculate bias values

    bias_measurements = BiasMeasurements()

    if wm <= 0:
        bias_measurements.m = 0
        bias_measurements.m_err = np.nan
    else:
        bias_measurements.m = (lm * lwm)[~mc_mask].sum() / wm
        bias_measurements.m_err = wm ** (-0.5)

    if wc <= 0:
        bias_measurements.c = 0
        bias_measurements.c_err = np.nan
    else:
        bias_measurements.c = (lc * lwc)[~mc_mask].sum() / wc
        bias_measurements.c_err = wc ** (-0.5)

    if bad_covar or wmc <= 0:
        bias_measurements.mc_covar = np.nan
    else:
        bias_measurements.mc_covar = wmc ** (-1)

    return bias_measurements


def calculate_bootstrap_bias_measurements_from_statistics(statistics_list, n_bootstrap=1000, seed=0):
    """Calculates a BiasMeasurements object using bootstrap errors from a list of statistics objects.
    """
    return calculate_bootstrap_bias_measurements(objects_list=statistics_list, combine_func=combine_bias_statistics,
                                                 n_bootstrap=n_bootstrap, seed=seed)


def calculate_bootstrap_bias_measurements_from_measurements(measurements_list, n_bootstrap=1000, seed=0):
    """Calculates a BiasMeasurements object using bootstrap errors from a list of measurements objects.
    """
    return calculate_bootstrap_bias_measurements(objects_list=measurements_list, combine_func=combine_bias_measurements,
                                                 n_bootstrap=n_bootstrap, seed=seed)


def calculate_bootstrap_bias_measurements(objects_list, combine_func, n_bootstrap=1000, seed=0):
    """Calculates a BiasMeasurements object using bootstrap errors from a list of statistics objects.
    """

    # Seed the random number generator
    np.random.seed(seed)

    # Get a base object for the m and c calculations
    bias_measurements = combine_func(objects_list)

    # Bootstrap to get errors on m and c
    n_sample = len(objects_list)

    m_bs = np.empty(n_bootstrap)
    c_bs = np.empty(n_bootstrap)
    objects_array = np.array(objects_list)
    for i in range(n_bootstrap):
        u = np.random.random_integers(0, n_sample - 1, n_sample)
        bias_measurements_bs = combine_func(objects_array[u])
        m_bs[i] = bias_measurements_bs.m
        c_bs[i] = bias_measurements_bs.c

    # Update the bias measurements in the output object
    bias_measurements.m_err = np.std(m_bs)
    bias_measurements.c_err = np.std(c_bs)

    return bias_measurements


def calculate_bfd_bootstrap_bias_measurements_from_statistics(statistics_list, n_bootstrap=1000, seed=0):
    """Calculates a BiasMeasurements object using bootstrap errors from a list of BFD statistics objects.
    """

    # Seed the random number generator
    np.random.seed(seed)

    # Get a base object for the m and c calculations
    g1_bias_measurements = BiasMeasurements(combine_bfd_sum_statistics(statistics_list, do_g1=True))
    g2_bias_measurements = BiasMeasurements(combine_bfd_sum_statistics(statistics_list, do_g1=False))

    # Bootstrap to get errors on m and c
    n_sample = len(statistics_list)

    m1_bs = np.empty(n_bootstrap)
    c1_bs = np.empty(n_bootstrap)
    m2_bs = np.empty(n_bootstrap)
    c2_bs = np.empty(n_bootstrap)
    statistics_array = np.array(statistics_list)
    for i in range(n_bootstrap):
        u = np.random.random_integers(0, n_sample - 1, n_sample)

        bias_measurements_b1s = BiasMeasurements(combine_bfd_sum_statistics(statistics_array[u], do_g1=True))
        bias_measurements_b2s = BiasMeasurements(combine_bfd_sum_statistics(statistics_array[u], do_g1=False))

        m1_bs[i] = bias_measurements_b1s.m
        c1_bs[i] = bias_measurements_b1s.c

        m2_bs[i] = bias_measurements_b2s.m
        c2_bs[i] = bias_measurements_b2s.c

    # Update the bias measurements in the output objects
    g1_bias_measurements.m_err = np.std(m1_bs)
    g1_bias_measurements.c_err = np.std(c1_bs)
    g2_bias_measurements.m_err = np.std(m2_bs)
    g2_bias_measurements.c_err = np.std(c2_bs)

    return g1_bias_measurements, g2_bias_measurements
