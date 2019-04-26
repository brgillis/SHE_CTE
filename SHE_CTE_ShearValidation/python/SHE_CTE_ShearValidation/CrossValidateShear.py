""" @file CrossValidateShear.py

    Created 12 Oct 2017

    Executable for cross-validating shear measurements and creating a combined catalogue.
"""

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

from ElementsKernel.Logging import getLogger
import SHE_CTE
from SHE_CTE.magic_values import force_dry_run
from SHE_CTE_ShearValidation import magic_values as mv
from SHE_CTE_ShearValidation.cross_validate_shear import cross_validate_shear
from SHE_PPT.utility import get_arguments_string


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program, using all possible configurations.

    @return
        An  ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_CrossValidateShear defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Option for profiling
    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')
    parser.add_argument('--dry_run', action='store_true',
                        help='Dry run (no data processed).')

    # Input filenames
    parser.add_argument('--shear_estimates_product', type=str,
                        help='Filename for shear estimates data product (XML data product)')

#     parser.add_argument('--shear_validation_statistics_table',type=str, # Disabled until it exists
#                         help='Filename for table of shear validation statistics.')

    # Output filenames
    parser.add_argument('--cross_validated_shear_estimates_product', type=str,
                        help='Desired filename for output shear estimates table (multi-HDU fits file).')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    # Optional arguments (can't be used with pipeline runner)
    parser.add_argument('--primary_method', type=str, default="LensMC",
                        help="Shear measurement method to consider primary, and compare against others.")

    logger.debug('Exiting SHE_CTE_CrossValidateShear defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_CrossValidateShear mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_CrossValidateShear",
                                    store_true=["profile", "debug", "dry_run"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    dry_run = args.dry_run or force_dry_run

    if args.profile:
        import cProfile
        cProfile.runctx("cross_validate_shear(args,dry_run=dry_run)", {},
                        {"cross_validate_shear": cross_validate_shear,
                         "args": args,
                         "dry_run": dry_run, },
                        filename="cross_validate_shear.prof")
    else:
        cross_validate_shear(args)

    logger.debug('Exiting SHE_CTE_CrossValidateShear mainMethod()')

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
