""" File: python/SHE_CTE_ShearReconciliation/ReconcileShear.py

    Created on: 3 August 2020
    Author: Bryan Gillis
"""

__updated__ = "2021-08-18"

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

from typing import Any, Dict, Tuple, Type, Union

from SHE_CTE.executor import CteLogOptions, SheCteExecutor
from SHE_CTE_ShearReconciliation.reconcile_shear import reconcile_shear_from_args
from SHE_PPT.argument_parser import SheArgumentParser
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.executor import ReadConfigArgs
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (ConfigKeys, ReconciliationConfigKeys)

# Set up dicts for pipeline config defaults and types
D_REC_SHEAR_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    ReconciliationConfigKeys.REC_METHOD       : "InvVar",
    ReconciliationConfigKeys.CHAINS_REC_METHOD: "Keep",
    }

D_REC_SHEAR_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    ReconciliationConfigKeys.REC_METHOD       : str,
    ReconciliationConfigKeys.CHAINS_REC_METHOD: str,
    }

D_REC_SHEAR_CONFIG_CLINE_ARGS: Dict[ConfigKeys, str] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    ReconciliationConfigKeys.REC_METHOD       : "method",
    ReconciliationConfigKeys.CHAINS_REC_METHOD: "chains_method",
    }

EXEC_NAME = "SHE_CTE_ReconcileShear"

logger = getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details
        See the Elements documentation for more details.
    @return
        An  ArgumentParser.
    """

    logger.debug('#')
    logger.debug(f'# Entering {EXEC_NAME} defineSpecificProgramOptions()')
    logger.debug('#')

    parser = SheArgumentParser()

    # Required input arguments

    parser.add_input_arg('--she_validated_measurements_listfile', type = str,
                         help = '.json listfile containing filenames of DpdSheValidatedMeasurements data products ' +
                                'to be reconciled and combined.')

    parser.add_input_arg('--she_lensmc_chains_listfile', type = str,
                         help = '.json listfile containing filenames of DpdSheLensMcChains data products ' +
                                'to be reconciled and combined.')

    parser.add_input_arg('--mer_final_catalog', type = str,
                         help = 'DpdMerFinalCatalog data product for this tile, which is used to determine the objects '
                                'to include in the output catalog.')

    # TODO: Remove this, and update port in pipeline script and package def
    parser.add_input_arg('--she_reconciliation_config', type = str,
                         help = 'DpdSheReconciliationConfig data product, which stores configuration options for this '
                                'executable, such as the reconciliation method to use.')

    # Output arguments

    parser.add_output_arg('--she_reconciled_measurements', type = str,
                          help = 'Desired filename to contain the output DpdSheReconciledMeasurements data product.')

    parser.add_output_arg('--she_reconciled_lensmc_chains', type = str,
                          help = 'Desired filename to contain the output DpdSheReconciledLensMcChains data product.')

    # Optional input arguments (cannot be used in pipeline)

    parser.add_option_arg('--method', type = str,
                          help = "Which reconciliation method to use. If not specified, will take value from pipeline "
                                 "config if available. If that's not available either, will use default of InvVar " +
                                 "(Inverse Variance) combination.")

    parser.add_option_arg('--chains_method', type = str,
                          help = "Which chains reconciliation method to use. If not specified, will take value from "
                                 "pipeline config if available. If that's not available either, will use default of "
                                 "keeping all chains.")

    logger.debug(f'# Exiting {EXEC_NAME} defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    executor = SheCteExecutor(run_from_args_function = reconcile_shear_from_args,
                              log_options = CteLogOptions(executable_name = EXEC_NAME),
                              config_args = ReadConfigArgs(d_config_defaults = D_REC_SHEAR_CONFIG_DEFAULTS,
                                                           d_config_types = D_REC_SHEAR_CONFIG_TYPES,
                                                           d_config_cline_args = D_REC_SHEAR_CONFIG_CLINE_ARGS,
                                                           s_config_keys_types = {ReconciliationConfigKeys},
                                                           ))

    executor.run(args, logger = logger, pass_args_as_dict = False)


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
