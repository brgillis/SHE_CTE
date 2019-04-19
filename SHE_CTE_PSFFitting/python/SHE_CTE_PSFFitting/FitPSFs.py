""" @file ModelPSFs.py

    Created 12 Oct 2017

    Executable for fitting PSFs and generating images of them.
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
from SHE_CTE.magic_values import force_dry_run
from SHE_CTE_PSFFitting import magic_values as mv
from SHE_CTE_PSFFitting.fit_psfs import fit_psfs
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
    logger.debug('# Entering SHE_CTE_FitPSFs defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Option for profiling
    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')
    parser.add_argument('--dry_run', action='store_true',
                        help='Dry run (no data processed).')

    # Input filenames
    parser.add_argument('--data_images', type=str,
                        help='Filename for list of data images (.json listfile)')
    parser.add_argument('--detections_tables', type=str,
                        help='Filename for list of detections tables (.json listfile)')
    parser.add_argument('--aocs_time_series_products', type=str, default=None,  # Allowing to use without for SC4
                        help='Filename for list of AOCS data series data products (.json listfile)')
    parser.add_argument('--psf_calibration_product', type=str, default=None,  # Allowing to use without for SC4
                        help='Filename for PSF calibration data product (.xml data product)')
    parser.add_argument('--segmentation_images', type=str, default=None,
                        help='Filename for the list of segmentation maps (.json listfile)')
    parser.add_argument('--mdb', type=str, default=None,
                        help='Filename for the Mission Database (MDB) file to use as input')

    # Output filenames
    parser.add_argument('--psf_field_params', type=str,
                        help='Desired filename for output PSF field parameters listfile.')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('Exiting SHE_CTE_FitPSFs defineSpecificProgramOptions()')

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
    logger.debug('# Entering SHE_CTE_FitPSFs mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE 0.7 SHE_CTE_FitPSFs",
                                    store_true=["profile", "dry_run"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    dry_run = args.dry_run or force_dry_run

    if args.profile:
        import cProfile
        cProfile.runctx("fit_psfs(args,dry_run=dry_run)", {},
                        {"fit_psfs": fit_psfs,
                         "args": args,
                         "dry_run": dry_run, },
                        filename="fit_psfs.prof")
    else:
        fit_psfs(args, dry_run=dry_run)

    logger.debug('Exiting SHE_CTE_FitPSFs mainMethod()')

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
