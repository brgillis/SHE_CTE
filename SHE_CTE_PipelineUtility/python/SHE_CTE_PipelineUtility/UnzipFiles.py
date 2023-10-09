#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#


"""
:file: python/SHE_CTE_PipelineUtility/UnzipFiles.py

:date: 31/08/23
:author: Gordon Gibb

"""

import argparse
import pathlib
import shutil
import gzip

from ElementsKernel import Logging
from ST_DM_DmUtils import DmUtils

logger = Logging.getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--workdir", type=str, default=".", help="Working directory")

    parser.add_argument("--logdir", type=str, default=".", help="Logging directory")

    parser.add_argument("--input_product", required=True, type=str, help="Product to unzip")

    parser.add_argument("--output_product", required=True, type=str, help="Name of product this is written to")

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger.info("#")
    logger.info("# Entering SHE_CTE_UnzipFiles mainMethod()")
    logger.info("#")

    workdir = pathlib.Path(args.workdir)

    convert_product(workdir, args.input_product, args.output_product)

    logger.info("#")
    logger.info("# Exiting SHE_CTE_UnzipFiles mainMethod()")
    logger.info("#")


def convert_product(workdir, input_product, output_product):
    """Takes an input product and (g)unzips all its gzipped files,
    outputting the product with corrected filenames in it"""

    logger.info("Reading input product %s", workdir / input_product)

    datadir = workdir / "data"

    dpd = DmUtils.read_product_metadata(workdir / input_product)

    # NOTE: This may not work when moving to xsData.
    # The below code was influenced by DmUtils.get_fits_file_names - once
    # the transition to xsData is complete, it may be worth looking at that
    # function to see how to implement the below code using xsData.

    dom = dpd.toDOM()

    for el in dom.getElementsByTagName("FileName"):
        fits_filename = el.firstChild.nodeValue
        fits_path = datadir / fits_filename

        if fits_path.suffix == ".gz":
            unzipped_fits_path = fits_path.with_suffix("")
            if unzipped_fits_path.exists():
                logger.info("Unzipped file %s already exists. Skipping.", unzipped_fits_path)
            else:
                try:
                    gunzip_file(fits_path, unzipped_fits_path)
                    logger.info("Gunzipped %s to %s", fits_path, unzipped_fits_path)
                except FileNotFoundError:
                    logger.warning("Cannot find file %s to be gunzipped - skipping", fits_path)

            # set the new name regardless (in case the file has already been unzipped somewhere else)
            el.firstChild.nodeValue = unzipped_fits_path.name

    output_filename = workdir / output_product

    logger.info("Writing product to %s", output_filename)
    with open(output_filename, "w") as f:
        f.write(dom.toprettyxml())


def gunzip_file(zipped_file, unzipped_file):
    """gunzips a file"""
    with gzip.open(zipped_file, "rb") as f_in:
        with open(unzipped_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
