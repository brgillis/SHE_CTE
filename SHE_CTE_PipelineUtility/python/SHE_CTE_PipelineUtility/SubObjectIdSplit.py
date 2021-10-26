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

import os
from typing import Any, Dict, Tuple, Type, Union

import SHE_CTE
from EL_PythonUtils.utilities import get_arguments_string
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import AnalysisConfigKeys, ConfigKeys, GlobalConfigKeys, read_analysis_config
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

    exec_cmd = get_arguments_string(args, cmd = f"E-Run SHE_CTE {SHE_CTE.__version__} SHE_CTE_SubObjectIdSplit",
                                    store_true = ["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    args.pipeline_config = read_analysis_config(config_filename = args.pipeline_config,
                                                workdir = args.workdir,
                                                d_defaults = D_SOID_SPLIT_CONFIG_DEFAULTS,
                                                d_cline_args = D_SOID_SPLIT_CONFIG_CLINE_ARGS,
                                                parsed_args = args,
                                                d_types = D_SOID_SPLIT_CONFIG_TYPES)

    # check if profiling is to be enabled from the pipeline config
    profiling = args.pipeline_config[GlobalConfigKeys.PIP_PROFILE]

    if args.profile or profiling:
        import cProfile

        logger.info("Profiling enabled")
        filename = os.path.join(args.workdir, args.logdir, "sub_object_id_split.prof")
        logger.info("Writing profiling data to %s", filename)

        cProfile.runctx("object_id_split_from_args(args,sub_batch=True)", {},
                        {"object_id_split_from_args": object_id_split_from_args,
                         "args"                     : args, },
                        filename = filename)
    else:
        logger.info("Profiling disabled")
        object_id_split_from_args(args,
                                  sub_batch = True)

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
