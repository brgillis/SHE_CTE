""" @file PrintBias.py

    Created 16 July 2018

    Main program for printing out bias of shear estimates
"""

__updated__ = "2019-01-22"

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

import argparse
import os

from SHE_PPT import products
from SHE_PPT.file_io import read_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string

from SHE_CTE_BiasMeasurement import magic_values as mv
import numpy as np


products.shear_bias_measurements.init()

methods = ["BFD",
           "KSB",
           "LensMC",
           "MomentsML",
           "REGAUSS"]


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_PrintBias defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')

    # Input data
    parser.add_argument('--shear_bias_measurements', type=str,
                        help='Desired name of the output shear bias statistics data product')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('# Exiting SHE_CTE_PrintBias defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to measure bias.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_PrintBias mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE 0.6.9 SHE_CTE_PrintBias",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    p = read_xml_product(os.path.join(args.workdir, args.shear_bias_measurements), allow_pickled=True)

    for method in methods:

        g1_bias, g2_bias = p.get_method_bias_measurements(method)

        if g1_bias is None or g2_bias is None:
            logger.warn("No bias measurements available for method " + method)
            continue

        if np.isnan(g1_bias.m) or np.isnan(g2_bias.m) or np.isinf(g1_bias.m) or np.isinf(g2_bias.m):
            logger.warn("Bad bias measurements for method " + method)
            continue

        logger.info('#')
        logger.info('# Bias measurements for method "' + method + '":')
        logger.info('#')

        for bias, label in ((g1_bias, "G1"),
                            (g2_bias, "G2")):
            logger.info(label + " component:")
            logger.info("m = " + str(bias.m) + " +/- " + str(bias.m_err))
            logger.info("c = " + str(bias.c) + " +/- " + str(bias.c_err))

    logger.debug('# Exiting SHE_CTE_PrintBias mainMethod()')

    return


def main():
    """
    @brief
        Alternate entry point for non-Elements execution.
    """

    parser = defineSpecificProgramOptions()

    args = parser.parse_args()

    mainMethod(args)

    return


if __name__ == "__main__":
    main()
