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
:file: python/SHE_CTE_PipelineUtility/ReorderListfile.py

:date: 19/02/24
:author: Gordon Gibb

This executable takes two input listfiles and matches their entries together.

"""

import argparse
import pathlib
import json

import ElementsKernel.Logging as log

from ST_DM_DmUtils import DmUtils


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    def directory(s):
        p = pathlib.Path(s)
        if not p.is_dir():
            raise ValueError(f"Directory {s} does not exit")
        return p

    parser = argparse.ArgumentParser()

    # Args required by the pipeline

    parser.add_argument("--workdir", type=directory, default=".", help="The workdir for the pipeline")

    parser.add_argument("--logdir", type=directory, default=".", help="The workdir for the pipeline")

    # Input file arguments

    parser.add_argument(
        "--reference_listfile", type=str, required=True, help="Listfile used as a reference when sorting"
    )

    parser.add_argument(
        "--input_listfile", type=str, required=True, help="Listfile to be sorted against the reference listfile"
    )

    # Output file arguments

    parser.add_argument("--output_listfile", type=str, required=True, help="Name of the sorted listfiler")

    # Options

    parser.add_argument(
        "--reference_path",
        type=str,
        required=True,
        help="Path inside the reference listfile to match against, e.g. Data.ObservationId",
    )

    parser.add_argument(
        "--input_path",
        type=str,
        required=True,
        help="Path inside the reference listfile to match against, e.g. Data.ObservationId",
    )

    parser.add_argument(
        "--allow_duplicates",
        action="store_true",
        help="If set, the reference product may have duplicate entries, and the output product may contain duplicate entries",
    )

    return parser


def rgetattr(input_obj, attr_string):
    """Recursive form of getattr. rgetattr(obj,"a.b.c") will return obj.a.b.c"""
    attrs = attr_string.split(".")

    obj = input_obj

    for attr in attrs:
        obj = getattr(obj, attr)

    return obj


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger("ReorderListfile")

    logger.info("#")
    logger.info("# Entering ReorderListfile mainMethod()")
    logger.info("#")

    workdir = args.workdir

    with open(workdir / args.reference_listfile) as f:
        reference_files = json.load(f)
    reference_products = [DmUtils.read_product_metadata(workdir / p) for p in reference_files]
    logger.info("Read %d products from reference listfile", len(reference_files))

    with open(workdir / args.input_listfile) as f:
        input_files = json.load(f)
    input_products = [DmUtils.read_product_metadata(workdir / p) for p in input_files]
    logger.info("Read %d projects from input listfile", len(input_files))

    reference_vals = [rgetattr(p, args.reference_path) for p in reference_products]
    logger.info("Reference values are %s", str(reference_vals))

    input_vals = [rgetattr(p, args.input_path) for p in input_products]
    logger.info("Input values are %s", str(input_vals))

    if not args.allow_duplicates:
        if len(set(reference_vals)) != len(reference_vals):
            raise ValueError(f"There are duplicate reference values: {reference_vals}")

        if len(input_products) != len(reference_products):
            raise ValueError(
                f"Unexpected number of input products. Expected {len(reference_products)} but got {len(input_products)}"
            )

    if len(set(input_vals)) != len(input_vals):
        raise ValueError(f"Some of the input values are duplicates - unable to match. {input_vals}")

    input_dict = {value: prod for value, prod in zip(input_vals, input_files)}

    try:
        output_files = [input_dict[reference_val] for reference_val in reference_vals]
    except KeyError as e:
        raise KeyError("Missing key in reference product...") from e

    with open(workdir / args.output_listfile, "w") as f:
        json.dump(output_files, f)

    logger.info("#")
    logger.info("# Exiting ReorderListfile mainMethod()")
    logger.info("#")
