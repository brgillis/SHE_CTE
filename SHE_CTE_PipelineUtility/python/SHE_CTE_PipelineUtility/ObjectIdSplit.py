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

import argparse
import os
from typing import Any, Dict, Tuple, Type, Union

import SHE_CTE
from EL_PythonUtils.utilities import get_arguments_string
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import AnalysisConfigKeys, ConfigKeys, GlobalConfigKeys, read_analysis_config
from .object_id_split import object_id_split_from_args

# Set up dicts for pipeline config defaults and types
D_OID_SPLIT_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    AnalysisConfigKeys.OID_BATCH_SIZE : 400,
    AnalysisConfigKeys.OID_MAX_BATCHES: 0,
    AnalysisConfigKeys.OID_IDS        : None,
    }

D_OID_SPLIT_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    AnalysisConfigKeys.OID_BATCH_SIZE : int,
    AnalysisConfigKeys.OID_MAX_BATCHES: int,
    AnalysisConfigKeys.OID_IDS        : (list, int),
    }

D_OID_SPLIT_CONFIG_CLINE_ARGS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    AnalysisConfigKeys.OID_BATCH_SIZE : None,
    AnalysisConfigKeys.OID_MAX_BATCHES: None,
    AnalysisConfigKeys.OID_IDS        : None,
    }

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

    parser = argparse.ArgumentParser()

    # Input arguments
    parser.add_argument('--mer_final_catalog_tables', type = str)

    parser.add_argument('--data_images', type = str, default = None,
                        help = '.json listfile containing filenames of data image products.')

    parser.add_argument("--pipeline_config", default = None, type = str,
                        help = "Pipeline-wide configuration file.")

    # Output arguments
    parser.add_argument('--object_ids', type = str, default = "object_ids.json",
                        help = 'Desired filename for output .json listfile of object ID list products.')
    parser.add_argument('--batch_mer_catalogs', type = str, default = "batch_mer_catalogs.json",
                        help = 'Desired filename for output .json listfile of MER catalogues for each batch of '
                               'objects.')

    # Required pipeline arguments
    parser.add_argument('--workdir', type = str, )
    parser.add_argument('--logdir', type = str, )
    parser.add_argument('--debug', action = 'store_true',
                        help = "Set to enable debugging protocols")
    parser.add_argument('--profile', action = 'store_true')

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

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_ObjectIdSplit mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd = f"E-Run SHE_CTE {SHE_CTE.__version__} SHE_CTE_ObjectIdSplit",
                                    store_true = ["profile", "debug"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    args.pipeline_config = read_analysis_config(config_filename = args.pipeline_config,
                                                workdir = args.workdir,
                                                d_defaults = D_OID_SPLIT_CONFIG_DEFAULTS,
                                                d_cline_args = D_OID_SPLIT_CONFIG_CLINE_ARGS,
                                                parsed_args = args,
                                                d_types = D_OID_SPLIT_CONFIG_TYPES)

    # check if profiling is to be enabled from the pipeline config
    profiling = args.pipeline_config[GlobalConfigKeys.PIP_PROFILE]

    if args.profile or profiling:
        import cProfile

        logger.info("Profiling enabled")
        filename = os.path.join(args.workdir, args.logdir, "object_id_split.prof")
        logger.info("Writing profiling data to %s", filename)

        cProfile.runctx("object_id_split_from_args(args,sub_batch=False)", {},
                        {"object_id_split_from_args": object_id_split_from_args,
                         "args"                     : args, },
                        filename = filename)
    else:
        logger.info("Profiling disabled")
        object_id_split_from_args(args,
                                  sub_batch = False)

    logger.debug('# Exiting SHE_CTE_ObjectIdSplit mainMethod()')


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
