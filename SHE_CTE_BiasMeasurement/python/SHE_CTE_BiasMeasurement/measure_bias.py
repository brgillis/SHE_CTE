""" @file measure_bias.py

    Created 7 Apr 2017

    Primary execution loop for measuring bias in shear estimates.
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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

from os.path import join

from SHE_PPT import products
from SHE_PPT.file_io import read_listfile, read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.math import combine_linregress_statistics, BiasMeasurements

from SHE_CTE_BiasMeasurement import magic_values as mv

import numpy as np


products.shear_bias_stats.init()
products.shear_bias_measurements.init()


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
    shear_statistics_files = read_listfile(join(args.workdir, args.shear_bias_statistics))

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

        shear_statistics_prod = read_xml_product(join(args.workdir, shear_statistics_file))

        for method in mv.estimation_methods:
            method_shear_statistics = shear_statistics_prod.get_method_statistics(method)
            if method_shear_statistics[0] is not None:
                method_shear_statistics_lists[method].g1_statistics_list.append(method_shear_statistics[0])
            if method_shear_statistics[1] is not None:
                method_shear_statistics_lists[method].g2_statistics_list.append(method_shear_statistics[1])

    # Calculate the bias and compile into a data product
    bias_measurement_prod = products.shear_bias_measurements.create_shear_bias_measurements_product()

    for method in mv.estimation_methods:

        if len(method_shear_statistics_lists[method].g1_statistics_list) > 50:
            # We have enough data to calculate bootstrap errors
            g1_bias_measurements = calculate_bootstrap_bias_measurements(
                method_shear_statistics_lists[method].g1_statistics_list, seed=args.bootstrap_seed)

        elif len(method_shear_statistics_lists[method].g1_statistics_list) > 0:
            # Not enough for bootstrap errors - calculate simply
            g1_bias_measurements = BiasMeasurements(
                combine_linregress_statistics(method_shear_statistics_lists[method].g1_statistics_list))
        else:
            g1_bias_measurements = None

        if len(method_shear_statistics_lists[method].g2_statistics_list) > 50:
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

    write_xml_product(bias_measurement_prod, join(args.workdir, args.shear_bias_measurements))

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
    for i in range(n_bootstrap):
        u = np.random.random_integers(0, n_sample - 1, n_sample)
        bias_measurements_bs = BiasMeasurements(combine_linregress_statistics(statistics_list[u]))
        m_bs[i] = bias_measurements_bs.m
        c_bs[i] = bias_measurements_bs.c

    # Update the bias measurements in the output object
    bias_measurements.m_err = np.std(m_bs)
    bias_measurements.c_err = np.std(c_bs)

    return bias_measurements
