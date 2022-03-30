""" @file SubObjectIdSplit.py

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
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.executor import ReadConfigArgs, RunArgs
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import AnalysisConfigKeys, ConfigKeys
from .ObjectIdSplit import defineSpecificProgramOptions
from .object_id_split import object_id_split_from_args

# Set up dicts for pipeline config defaults and types
D_SOID_SPLIT_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    AnalysisConfigKeys.SOID_BATCH_SIZE : 20,
    AnalysisConfigKeys.SOID_MAX_BATCHES: 0,
    AnalysisConfigKeys.SOID_IDS        : None,
    }

D_SOID_SPLIT_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    AnalysisConfigKeys.SOID_BATCH_SIZE : int,
    AnalysisConfigKeys.SOID_MAX_BATCHES: int,
    AnalysisConfigKeys.SOID_IDS        : (list, int),
    }

D_SOID_SPLIT_CONFIG_CLINE_ARGS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    AnalysisConfigKeys.SOID_BATCH_SIZE : None,
    AnalysisConfigKeys.SOID_MAX_BATCHES: None,
    AnalysisConfigKeys.SOID_IDS        : None,
    }

EXEC_NAME = "SHE_CTE_SubObjectIdSplit"

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

    executor = SheCteExecutor(run_from_args_function = object_id_split_from_args,
                              log_options = CteLogOptions(executable_name = EXEC_NAME),
                              config_args = ReadConfigArgs(d_config_defaults = D_SOID_SPLIT_CONFIG_DEFAULTS,
                                                           d_config_types = D_SOID_SPLIT_CONFIG_TYPES,
                                                           d_config_cline_args = D_SOID_SPLIT_CONFIG_CLINE_ARGS,
                                                           s_config_keys_types = {AnalysisConfigKeys},
                                                           ),
                              run_args = RunArgs(d_run_kwargs = {"sub_batch": True}))

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
