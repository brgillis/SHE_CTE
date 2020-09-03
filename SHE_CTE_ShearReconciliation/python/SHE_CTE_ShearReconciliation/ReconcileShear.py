""" File: python/SHE_CTE_ShearReconciliation/ReconcileShear.py
    
    Created on: 3 August 2020
    Author: Bryan Gillis
"""

__updated__ = "2020-09-03"

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

import argparse

from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string

import SHE_CTE
from SHE_CTE.magic_values import force_dry_run
from SHE_CTE_ShearReconciliation.reconcile_shear import reconcile_shear_from_args

logger = getLogger('SHE_CTE_ReconcileShear')


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details
        See the Elements documentation for more details.
    @return
        An  ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')
    parser.add_argument('--dry_run', action='store_true',
                        help='Dry run (no data processed).')
    parser.add_argument('--debug', action='store_true',
                        help='Enables debug mode - currently no functional difference.')

    # Required input arguments

    parser.add_argument('--she_validated_measurements_listfile', type=str,
                        help='.json listfile containing filenames of DpdSheValidatedMeasurements data products ' +
                             'to be reconciled and combined.')

    parser.add_argument('--she_lensmc_chains_listfile', type=str,
                        help='.json listfile containing filenames of DpdSheLensMcChains data products ' +
                             'to be reconciled and combined.')

    parser.add_argument('--mer_final_catalog', type=str,
                        help='DpdMerFinalCatalog data product for this tile, which is used to determine the objects ' +
                             'to include in the output catalog.')

    parser.add_argument('--she_reconciliation_config', type=str, default=None,
                        help='DpdSheReconciliationConfig data product, which stores configuration options for this ' +
                             'executable, such as the reconciliation method to use.')

    # Optional input arguments (cannot be used in pipeline)

    parser.add_argument('--method', type=str, default=None,
                        help="Which reconciliation method to use. If not specified, will take value from pipeline " +
                             "config if available. If that's not available either, will use default of InvVar " +
                             "(Inverse Variance) combination.")

    # Output arguments

    parser.add_argument('--she_reconciled_measurements', type=str,
                        help='Desired filename to contain the output DpdSheReconciledMeasurements data product.')

    parser.add_argument('--she_reconciled_lensmc_chains', type=str,
                        help='Desired filename to contain the output DpdSheReconciledLensMcChains data product.')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    return parser


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ReconcileShear mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_ReconcileShear",
                                    store_true=["profile", "debug", "dry_run"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    dry_run = args.dry_run or force_dry_run

    if args.profile:
        import cProfile
        cProfile.runctx("reconcile_shear_from_args(args)", {},
                        {"reconcile_shear_from_args": reconcile_shear_from_args,
                         "args": args,
                         "dry_run": dry_run},
                        filename="reconcile_shear.prof")
    else:
        reconcile_shear_from_args(args, dry_run)

    logger.debug('# Exiting SHE_CTE_ReconcileShear mainMethod() successfully.')
