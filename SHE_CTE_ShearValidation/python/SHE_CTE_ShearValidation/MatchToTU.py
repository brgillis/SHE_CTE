""" @file MatchToTu.py

    Created 10 May 2019

    Executable for matching the output of the Analysis pipeline to SIM's True Universe catalogs.
"""

__updated__ = "2021-02-12"

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

from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string

import SHE_CTE
from SHE_CTE_ShearValidation.match_to_tu import match_to_tu_from_args


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program, using all possible configurations.

    @return
        An  ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MatchToTU defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Input filenames
    parser.add_argument('--she_measurements_product', type=str,
                        help='Filename for shear estimates data product (XML data product)')
    parser.add_argument('--tu_galaxy_catalog_list', type=str, default=None,
                        help='Filename for True Universe Galaxy Catalog listfile (.json).')
    parser.add_argument('--tu_star_catalog_list', type=str, default=None,
                        help='Filename for True Universe Star Catalog listfile (.json).')
    parser.add_argument('--tu_galaxy_catalog', type=str, default=None,
                        help='Filename for True Universe Galaxy Catalog data product (XML data product)')
    parser.add_argument('--tu_star_catalog', type=str, default=None,
                        help='Filename for True Universe Star Catalog data product (XML data product)')

    # Output filenames
    parser.add_argument('--matched_catalog', type=str,
                        help='Desired filename for output matched catalog data product (XML data product).')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    # Optional arguments (can't be used with pipeline runner)
    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')
    parser.add_argument('--sim_path', type=str, default="/mnt/cephfs/share/SC8/SIM",
                        help="Path to where the SIM data is stored")
    parser.add_argument('--match_threshold', type=float, default=0.3,
                        help="Maximum distance allowed for a match in units of arcsec.")

    logger.debug('Exiting SHE_CTE_MatchToTU defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to generate galaxy images.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_MatchToTU mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_MatchToTU",
                                    store_true=["profile"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    if args.profile:
        import cProfile
        cProfile.runctx("match_to_tu_from_args(args,dry_run=dry_run)", {},
                        {"match_to_tu_from_args": match_to_tu_from_args,
                         "args": args, },
                        filename="match_to_tu_from_args.prof")
    else:
        match_to_tu_from_args(args)

    logger.debug('Exiting SHE_CTE_MatchToTU mainMethod()')

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
