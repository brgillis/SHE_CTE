""" @file ObjectIdSplit.py

    Created 14 Mar 2019

    Split point executable for splitting up processing of objects into batches.
"""

__updated__ = "2021-08-18"

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

from SHE_PPT.argument_parser import dir_path
from SHE_PPT.logging import getLogger

from .object_id_split import object_id_split_from_args

logger = getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ObjectIdSplit defineSpecificProgramOptions()')
    logger.debug('#')

    def integer_list(s):
        """
        Takes a comma separated string (e.g. "1,2,3,4") and returns a list of integers from
        that string, (e.g. [1,2,3,4])
        """
        if len(s) == 0:
            return []
        return [int(i) for i in s.split(",")]

    parser = argparse.ArgumentParser()

    # Pipeline args

    parser.add_argument("--workdir", type=dir_path, default=".", help="Workdir")

    parser.add_argument("--logdir", type=dir_path, default=".", help="Logging directory")

    # Input arguments

    parser.add_argument('--mer_final_catalog_tables', required=True, type=str)

    parser.add_argument(
        '--data_images', type=str, required=True,
        help='.json listfile containing filenames of data image products.'
    )

    # Output arguments

    parser.add_argument(
        '--object_ids', type=str, required=True,
        help='Desired filename for output .json listfile of object ID list products.'
    )

    parser.add_argument(
        '--batch_mer_catalogs', type=str, required=True,
        help='Desired filename for output .json listfile of MER catalogues for each batch of objects.'
    )

    # option arguments

    parser.add_argument(
        "--batch_size", type=int,
        help="Maximum number of objects per batch. The actual returned batchsize will be at least half of this"
    )

    parser.add_argument(
        "--max_batches", type=int,
        help="Maximum number of batches. If this is zero, the maximum number is unlimited"
    )

    parser.add_argument(
        "--grouping_radius", type=float,
        help="The grouping radius in arcseconds (objects closer than this are put into groups)"
    )

    parser.add_argument(
        "--include_vis_non_detections", action="store_true",
        help="If set, objects with VIS_DET=0 will be included"
    )

    parser.add_argument(
        "--ids_to_use", type=integer_list,
        help="Comma separated list of object_ids to use. If unset, all objects are considered"
    )

    logger.debug('# Exiting SHE_CTE_ObjectIdSplit defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    object_id_split_from_args(args)


def main():
    """
    @brief
        Alternate entry point for non-Elements execution.
    """

    parser = defineSpecificProgramOptions()

    args = parser.parse_args()

    mainMethod(args)


if __name__ == "__main__":
    main()
