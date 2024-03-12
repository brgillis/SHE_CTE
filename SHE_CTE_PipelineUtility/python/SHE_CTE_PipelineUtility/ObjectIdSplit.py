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

from typing import Any, Dict, Tuple, Type, Union

from SHE_CTE.executor import CteLogOptions, SheCteExecutor
from SHE_PPT.argument_parser import SheArgumentParser
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.executor import ReadConfigArgs, RunArgs
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import AnalysisConfigKeys, ConfigKeys
from .object_id_split import object_id_split_from_args

# Set up dicts for pipeline config defaults and types
D_OID_SPLIT_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    AnalysisConfigKeys.OID_BATCH_SIZE : 400,
    AnalysisConfigKeys.OID_MAX_BATCHES: 0,
    AnalysisConfigKeys.OID_IDS        : None,
    AnalysisConfigKeys.OID_GROUPING_RADIUS: 1.0
    }

D_OID_SPLIT_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    AnalysisConfigKeys.OID_BATCH_SIZE : int,
    AnalysisConfigKeys.OID_MAX_BATCHES: int,
    AnalysisConfigKeys.OID_IDS        : (list, int),
    AnalysisConfigKeys.OID_GROUPING_RADIUS: float
    }

D_OID_SPLIT_CONFIG_CLINE_ARGS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    AnalysisConfigKeys.OID_BATCH_SIZE : None,
    AnalysisConfigKeys.OID_MAX_BATCHES: None,
    AnalysisConfigKeys.OID_IDS        : None,
    AnalysisConfigKeys.OID_GROUPING_RADIUS: None
    }

EXEC_NAME = "SHE_CTE_ObjectIdSplit"

logger = getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_(Sub)ObjectIdSplit defineSpecificProgramOptions()')
    logger.debug('#')

    parser = SheArgumentParser()

    # Input arguments
    parser.add_input_arg('--mer_final_catalog_tables', type = str)

    parser.add_input_arg('--data_images', type = str, default = None,
                         help = '.json listfile containing filenames of data image products.')

    # Output arguments
    parser.add_output_arg('--object_ids', type = str, default = "object_ids.json",
                          help = 'Desired filename for output .json listfile of object ID list products.')
    parser.add_output_arg('--batch_mer_catalogs', type = str, default = "batch_mer_catalogs.json",
                          help = 'Desired filename for output .json listfile of MER catalogues for each batch of '
                                 'objects.')
    
    parser.add_argument(
        "--skip_vis_non_detections", action="store_true",
        help="If set, objects with VIS_DET=0 will be excluded"
    )

    logger.debug('# Exiting SHE_CTE_(Sub)ObjectIdSplit defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    executor = SheCteExecutor(run_from_args_function = object_id_split_from_args,
                              log_options = CteLogOptions(executable_name = EXEC_NAME),
                              config_args = ReadConfigArgs(d_config_defaults = D_OID_SPLIT_CONFIG_DEFAULTS,
                                                           d_config_types = D_OID_SPLIT_CONFIG_TYPES,
                                                           d_config_cline_args = D_OID_SPLIT_CONFIG_CLINE_ARGS,
                                                           s_config_keys_types = {AnalysisConfigKeys},
                                                           ),
                              run_args = RunArgs(d_run_kwargs = {"sub_batch": False}))

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
