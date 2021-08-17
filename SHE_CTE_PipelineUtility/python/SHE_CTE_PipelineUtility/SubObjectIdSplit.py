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

import os

from SHE_PPT.logging import getLogger
from EL_PythonUtils.utilities import get_arguments_string
from SHE_PPT.pipeline_utility import read_analysis_config, AnalysisConfigKeys

import SHE_CTE
from .ObjectIdSplit import defineSpecificProgramOptions
from .object_id_split import object_id_split_from_args, defaults_dict_sub_batch


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

    pipeline_config = read_analysis_config(config_filename=args.pipeline_config,
                                            workdir=args.workdir,
                                            cline_args=None,
                                            defaults=defaults_dict_sub_batch)
    
    #set args.pipeline_config to the read-in pipeline_config
    args.pipeline_config = pipeline_config
    
    #check if profiling is to be enabled from the pipeline config
    profiling = pipeline_config[AnalysisConfigKeys.PIP_PROFILE.value].lower() in ['true', 't']

    try:

        if args.profile or profiling:
            import cProfile

            logger.info("Profiling enabled")
            filename = os.path.join(args.workdir,args.logdir,"sub_object_id_split.prof")
            logger.info("Writing profiling data to %s",filename)

            cProfile.runctx("object_id_split_from_args(args,sub_batch=True)", {},
                            {"object_id_split_from_args": object_id_split_from_args,
                             "args": args, },
                            filename=filename)
        else:
            logger.info("Profiling disabled")
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
