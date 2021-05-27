""" @file ObjectIdSplit.py

    Created 14 Mar 2019

    Split point executable for splitting up processing of objects into batches.
"""

__updated__ = "2021-05-27"

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
from SHE_CTE_PipelineUtility.object_id_split import object_id_split_from_args


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

    parser = argparse.ArgumentParser()

    # Input arguments
    parser.add_argument('--mer_final_catalog_tables', type=str)

    parser.add_argument('--data_images', type=str, default=None,
                        help='.json listfile containing filenames of data image products.')

    parser.add_argument("--pipeline_config", default=None, type=str,
                        help="Pipeline-wide configuration file.")

    # Output arguments
    parser.add_argument('--object_ids', type=str, default="object_ids.json",
                        help='Desired filename for output .json listfile of object ID list products.')
    parser.add_argument('--batch_mer_catalogs', type=str, default="batch_mer_catalogs.json",
                        help='Desired filename for output .json listfile of MER catalogues for each batch of objects.')

    # Required pipeline arguments
    parser.add_argument('--workdir', type=str,)
    parser.add_argument('--logdir', type=str,)
    parser.add_argument('--debug', action='store_true',
                        help="Set to enable debugging protocols")
    parser.add_argument('--profile', action='store_true')

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

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ObjectIdSplit mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd=f"E-Run SHE_CTE {SHE_CTE.__version__} SHE_CTE_ObjectIdSplit",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    try:

        if args.profile:
            import cProfile
            cProfile.runctx("object_id_split_from_args(args,sub_batch=False)", {},
                            {"object_id_split_from_args": object_id_split_from_args,
                             "args": args, },
                            filename="object_id_split.prof")
        else:
            object_id_split_from_args(args,
                                      sub_batch=False)

    except Exception as e:
        # logger.warning(f"Failsafe exception block triggered with exception: {e}")
        raise

    logger.debug('# Exiting SHE_CTE_ObjectIdSplit mainMethod()')

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
