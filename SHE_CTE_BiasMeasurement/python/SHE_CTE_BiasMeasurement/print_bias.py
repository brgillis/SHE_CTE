""" @file PrintBias.py

    Created 16 July 2018

    Main program for printing out bias of shear estimates
"""

__updated__ = "2019-12-11"

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

import SHE_CTE
from SHE_CTE_BiasMeasurement import magic_values as mv
import numpy as np


methods = ["BFD",
           "KSB",
           "LensMC",
           "MomentsML",
           "REGAUSS"]

logger = getLogger(__name__)

def print_bias_from_product_filename(product_filename, workdir):

    p = read_xml_product(os.path.join(workdir, product_filename), allow_pickled=True)
    
    print_bias_from_product(p, workdir)
    
    return

def print_bias_from_product(p, workdir):

    for method in methods:

        g1_bias, g2_bias = p.get_method_bias_measurements(method, workdir=workdir)

        if (g1_bias is None or g2_bias is None or g1_bias.m is None or g2_bias.m is None or
                g1_bias.m == "" or g2_bias.m == ""):
            logger.warn('#')
            logger.warn("No bias measurements available for method " + method)
            logger.warn('#')
            continue

        if (np.isnan(g1_bias.m) or np.isnan(g2_bias.m) or np.isinf(g1_bias.m) or np.isinf(g2_bias.m) or
                g1_bias.m == "NaN" or g2_bias.m == "NaN" or g1_bias.m == "Inf" or g2_bias.m == "Inf"):
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

    return
