""" @file SubObjectIdSplit.py

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

from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string

import SHE_CTE
from .ObjectIdSplit import defineSpecificProgramOptions
from .object_id_split import object_id_split_from_args


logger = getLogger(__name__)


# defineSpecificProgramOptions is imported from ObjectIdSplit in order to keep a constant interface


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_SubObjectIdSplit mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd=f"E-Run SHE_CTE {SHE_CTE.__version__} SHE_CTE_SubObjectIdSplit",
                                    store_true=["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    try:

        if args.profile:
            import cProfile
            cProfile.runctx("object_id_split_from_args(args,sub_batch=True)", {},
                            {"object_id_split_from_args": object_id_split_from_args,
                             "args": args, },
                            filename="sub_object_id_split.prof")
        else:
            object_id_split_from_args(args,
                                      sub_batch=True)

    except Exception as e:
        # logger.warning(f"Failsafe exception block triggered with exception: {e}")
        raise

    logger.debug('# Exiting SHE_CTE_SubObjectIdSplit mainMethod()')


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
