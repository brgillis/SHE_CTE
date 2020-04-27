""" @file bfd_integrate.py

    Created 03 Dec 2019

    Performs BFD Integration step
"""

__updated__ = "2020-04-27"


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

import copy
import os
import pdb
from SHE_PPT import magic_values as ppt_mv
from SHE_PPT import mdb
from SHE_PPT import products
from SHE_PPT.file_io import (read_xml_product, write_xml_product, get_allowed_filename, get_data_filename,
                             read_listfile, find_file)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import ConfigKeys, read_config, get_conditional_product
from SHE_PPT.table_formats.bfd_moments import initialise_bfd_moments_table, tf as setf_bfd
from SHE_PPT.table_utility import is_in_format, table_to_hdu
from SHE_PPT.utility import hash_any
from astropy.io import fits
from astropy.table import Table

import SHE_CTE
from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_CTE_ShearEstimation.bfd_functions import bfd_measure_moments, bfd_perform_integration, bfd_load_training_data
import numpy as np

estimation_methods = ["KSB",
                      "REGAUSS",
                      "MomentsML",
                      "LensMC",
                      "BFD"]

def perform_bfd_integration(args, dry_run=False):
    """
    @brief
        Perform integration step for BFD method to produce probabilities

    @param kwargs <dict>

    @return None
    """

    logger = getLogger(__name__)

    logger.debug("Entering perform_bfd_integration")

    # Load in the files in turn

    if dry_run:
        dry_label = " mock dry"
    else:
        dry_label = ""

    # Load in the MDB
    if args.mdb is None:
        logger.warn("No MDB file provided as input. Default values will be used where necessary.")
        mdb.init(find_file("WEB/SHE_PPT/sample_mdb.xml"))
    elif args.mdb[-5:] == ".json":
        mdb_files = read_listfile(os.path.join(args.workdir, args.mdb))
        mdb.init(mdb_files=mdb_files, path=args.workdir)
    elif args.mdb[-4:] == ".xml":
        mdb.init(mdb_files=args.mdb, path=args.workdir)
    else:
        logger.warn("Unrecognized format for MDB file: " + os.path.splitext(args.mdb)[-1] +
                    ". Expected '.xml' or '.json'. Will attempt to proceed with default values.")

    if not dry_run:
        # Set up BFD training Data
        
        bfd_training_data_filename=args.bfd_training_data
        bfd_training_data=bfd_load_training_data(bfd_training_data_filename,workdir=args.workdir)

        # Read in Shear Estimates Product
        shear_estimates_prod_table = read_xml_product(os.path.join(args.workdir,args.shear_estimates_product))


        method_shear_estimates = {}

        # Check for methods in the pipeline options
        pipeline_config = read_config(args.pipeline_config, workdir=args.workdir)

        # Use methods specified in the command-line first
        if args.methods is not None and len(args.methods) > 0:
            methods = args.methods
        elif ConfigKeys.ES_METHODS.value in pipeline_config:
            methods = pipeline_config[ConfigKeys.ES_METHODS.value].split()
        else:
            # Default to using all methods
            methods = estimation_methods

        for method in methods:
            print(method)

            if method != "BFD":
                logger.info("Skipping method " + method + " because does not need BFD Integration")
                continue

            logger.info("Performing Integration for method " + method + "...")
            # get BFD filename for shear estimates

            shear_estimates_filename = shear_estimates_prod_table.get_method_filename(method)

            # Perform the Integration

            try:
                pmem = os.popen('ps -p ' + str(os.getpid()) + ' -o pmem').readlines()[-1].split()[0]
                logger.debug("Memory used before deletion: " + pmem + "%")
                pmem = os.popen('ps -p ' + str(os.getpid()) + ' -o pmem').readlines()[-1].split()[0]
                logger.debug("Memory used after deletion: " + pmem + "%")
                bfd_perform_integration(target_file=os.path.join(args.workdir, shear_estimates_filename),
                                        template_file=os.path.join(args.workdir,bfd_training_data))
            except Exception as e:
                logger.warn("Failsafe exception block triggered with exception: " + str(e))

    else:  # Dry run

        logger.info("Dry Run No BFD Integration Run")

    write_xml_product(shear_estimates_prod_table, args.shear_estimates_product_update, workdir=args.workdir)

    logger.info("Finished BFD Integration")

    return
