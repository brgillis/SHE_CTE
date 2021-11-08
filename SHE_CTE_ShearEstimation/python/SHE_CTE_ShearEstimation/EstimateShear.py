""" @file EstimateShears.py

    Created 27 Mar 2017

    Main program for estimating shears on simulation data.
"""

__updated__ = "2021-09-15"

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
from SHE_CTE_ShearEstimation.estimate_shears import estimate_shears_from_args
from SHE_PPT.argument_parser import SheArgumentParser
from SHE_PPT.constants.classes import ShearEstimationMethods
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (AnalysisConfigKeys, CalibrationConfigKeys, ConfigKeys, GlobalConfigKeys,
                                      read_config, )

# Set up dicts for pipeline config defaults and types
D_EST_SHEAR_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    AnalysisConfigKeys.ES_METHODS               : list(ShearEstimationMethods),
    AnalysisConfigKeys.ES_CHAINS_METHOD         : ShearEstimationMethods.LENSMC,
    AnalysisConfigKeys.ES_MEMMAP_IMAGES         : False,
    AnalysisConfigKeys.LENSMC_STAMP_SIZE        : 384,
    AnalysisConfigKeys.LENSMC_X_BUFFER          : 3,
    AnalysisConfigKeys.LENSMC_Y_BUFFER          : 3,
    AnalysisConfigKeys.LENSMC_DILATE_MASK       : True,
    AnalysisConfigKeys.LENSMC_HL_TO_EXP         : 0.15,
    AnalysisConfigKeys.LENSMC_N_BULGE           : 1,
    AnalysisConfigKeys.LENSMC_N_DISC            : 1,
    AnalysisConfigKeys.LENSMC_E_MAX             : 0.99,
    AnalysisConfigKeys.LENSMC_RE_MAX            : 2.,
    AnalysisConfigKeys.LENSMC_DELTA_MAX         : 0.3,
    AnalysisConfigKeys.LENSMC_E_FLAG            : 0.98,
    AnalysisConfigKeys.LENSMC_RE_FLAG           : 1.8,
    AnalysisConfigKeys.LENSMC_DELTA_FLAG        : 0.28,
    AnalysisConfigKeys.LENSMC_DISC_ONLY         : False,
    AnalysisConfigKeys.LENSMC_PSF_OVERSAMPLING  : 5,
    AnalysisConfigKeys.LENSMC_SEED              : None,
    AnalysisConfigKeys.LENSMC_SHAPE_NOISE       : 0.25,
    AnalysisConfigKeys.LENSMC_FAST_MODE         : False,
    AnalysisConfigKeys.LENSMC_ONLY_VIS_DETECTED : True,
    AnalysisConfigKeys.LENSMC_MONITOR           : False,
    }

D_EST_SHEAR_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    AnalysisConfigKeys.ES_METHODS               : (list, ShearEstimationMethods),
    AnalysisConfigKeys.ES_CHAINS_METHOD         : ShearEstimationMethods,
    AnalysisConfigKeys.ES_MEMMAP_IMAGES         : bool,
    AnalysisConfigKeys.LENSMC_STAMP_SIZE        : int,
    AnalysisConfigKeys.LENSMC_X_BUFFER          : int,
    AnalysisConfigKeys.LENSMC_Y_BUFFER          : int,
    AnalysisConfigKeys.LENSMC_DILATE_MASK       : bool,
    AnalysisConfigKeys.LENSMC_HL_TO_EXP         : float,
    AnalysisConfigKeys.LENSMC_N_BULGE           : float,
    AnalysisConfigKeys.LENSMC_N_DISC            : float,
    AnalysisConfigKeys.LENSMC_E_MAX             : float,
    AnalysisConfigKeys.LENSMC_RE_MAX            : float,
    AnalysisConfigKeys.LENSMC_DELTA_MAX         : float,
    AnalysisConfigKeys.LENSMC_E_FLAG            : float,
    AnalysisConfigKeys.LENSMC_RE_FLAG           : float,
    AnalysisConfigKeys.LENSMC_DELTA_FLAG        : float,
    AnalysisConfigKeys.LENSMC_DISC_ONLY         : bool,
    AnalysisConfigKeys.LENSMC_PSF_OVERSAMPLING  : int,
    AnalysisConfigKeys.LENSMC_SEED              : int,
    AnalysisConfigKeys.LENSMC_SHAPE_NOISE       : float,
    AnalysisConfigKeys.LENSMC_FAST_MODE         : bool,
    AnalysisConfigKeys.LENSMC_ONLY_VIS_DETECTED : bool,
    AnalysisConfigKeys.LENSMC_MONITOR           : bool,
    }

D_EST_SHEAR_CONFIG_CLINE_ARGS: Dict[ConfigKeys, str] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    AnalysisConfigKeys.ES_METHODS               : "methods",
    AnalysisConfigKeys.ES_CHAINS_METHOD         : "chains_method",
    AnalysisConfigKeys.ES_MEMMAP_IMAGES         : "memmap_images",
    AnalysisConfigKeys.LENSMC_STAMP_SIZE        : "lensmc_stamp_size",
    AnalysisConfigKeys.LENSMC_X_BUFFER          : "lensmc_x_buffer",
    AnalysisConfigKeys.LENSMC_Y_BUFFER          : "lensmc_y_buffer",
    AnalysisConfigKeys.LENSMC_DILATE_MASK       : "lensmc_dilate_mask",
    AnalysisConfigKeys.LENSMC_HL_TO_EXP         : "lensmc_hl_to_exp",
    AnalysisConfigKeys.LENSMC_N_BULGE           : "lensmc_n_bulge",
    AnalysisConfigKeys.LENSMC_N_DISC            : "lensmc_n_disc",
    AnalysisConfigKeys.LENSMC_E_MAX             : "lensmc_e_max",
    AnalysisConfigKeys.LENSMC_RE_MAX            : "lensmc_re_max",
    AnalysisConfigKeys.LENSMC_DELTA_MAX         : "lensmc_delta_max",
    AnalysisConfigKeys.LENSMC_E_FLAG            : "lensmc_e_flag",
    AnalysisConfigKeys.LENSMC_RE_FLAG           : "lensmc_re_flag",
    AnalysisConfigKeys.LENSMC_DELTA_FLAG        : "lensmc_delta_flag",
    AnalysisConfigKeys.LENSMC_DISC_ONLY         : "lensmc_disc_only",
    AnalysisConfigKeys.LENSMC_PSF_OVERSAMPLING  : "lensmc_psf_oversampling",
    AnalysisConfigKeys.LENSMC_SEED              : "lensmc_seed",
    AnalysisConfigKeys.LENSMC_SHAPE_NOISE       : "lensmc_shape_noise",
    AnalysisConfigKeys.LENSMC_FAST_MODE         : "lensmc_fast_mode",
    AnalysisConfigKeys.LENSMC_ONLY_VIS_DETECTED : "lensmc_only_vis_detected",
    AnalysisConfigKeys.LENSMC_MONITOR           : "lensmc_monitor",
    }


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_EstimateShear defineSpecificProgramOptions()')
    logger.debug('#')

    parser = SheArgumentParser()

    parser.add_argument('--debug', action = 'store_true',
                        help = 'Enables debug mode - only process first 1000 galaxies.')

    # Required input arguments

    parser.add_argument('--data_images', type = str,
                        help = '.json listfile containing filenames of data image products.')

    parser.add_argument('--stacked_image', type = str,
                        help = 'Filename of the stacked image product')

    parser.add_argument('--psf_images_and_tables', type = str,
                        help = '.json listfile containing filenames of psf image products.')

    parser.add_argument('--segmentation_images', type = str,
                        help = '.json listfile containing filenames of segmentation map images.')

    parser.add_argument('--stacked_segmentation_image', type = str,
                        help = 'Filename of the stacked segmentation map image.')

    parser.add_argument('--detections_tables', type = str,
                        help = '.json listfile containing filenames of detections table products.')

    parser.add_argument('--object_ids', type = str,
                        help = 'XML dataproduct that contains within it a list of object ids to loop over')

    parser.add_argument('--ksb_training_data', type = str, default = None,
                        # Use default in case we don't use it for SC4
                        help = 'Data product for KSB training data.')

    parser.add_argument('--lensmc_training_data', type = str, default = None,
                        # Use default in case we don't use it for SC4
                        help = 'Data product for LensMC training data.')

    parser.add_argument('--momentsml_training_data', type = str, default = None,
                        # Use default in case we don't use it for SC4
                        help = 'Data product for MomentsML training data.')

    parser.add_argument('--regauss_training_data', type = str, default = None,
                        # Use default in case we don't use it for SC4
                        help = 'Data product for REGAUSS training data.')

    parser.add_argument('--mdb', type = str, default = None,  # Use default to allow simple running with default values
                        help = 'Mission Database .xml file')

    parser.add_argument('--galaxy_population_priors_table', type = str, default = None,
                        # Use default in case we don't use it for SC4
                        help = '.json listfile containing filenames of detections table products.')

    parser.add_argument('--calibration_parameters_product', type = str, default = None,
                        # Use default in case we don't use it for SC4
                        help = 'Filename of calibration parameters product (XML data product).')

    # Optional input arguments (cannot be used in pipeline)

    parser.add_argument('--methods', type = str, nargs = '*', default = None,
                        help = 'Which shear estimation methods to apply. If not specified, all will be run.')

    parser.add_argument('--chains_method', type = str, default = None,
                        help = 'Which shear estimation method to generate chains with.')

    parser.add_argument('--fast_mode', action = 'store_true',
                        help = 'Enable LensMC fast mode')

    parser.add_argument('--memmap_images', action = 'store_true',
                        help = 'Memory-map images rather than reading them on-demand.')

    # Output arguments

    parser.add_argument('--shear_estimates_product', type = str,
                        help = 'XML data product to contain file links to the shear estimates tables.')

    parser.add_argument('--she_lensmc_chains', type = str,
                        help = 'XML data product to contain LensMC chains data.')

    logger.debug('# Exiting SHE_CTE_EstimateShear defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to estimate shears.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_EstimateShear mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd = "E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_EstimateShear",
                                    store_true = ["profile", "debug", "dry_run", "fast_mode", "memmap_images"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    dry_run = args.dry_run

    # load the pipeline config in
    args.pipeline_config = read_config(args.pipeline_config, workdir = args.workdir,
                                       config_keys = (AnalysisConfigKeys, CalibrationConfigKeys),
                                       d_cline_args = D_EST_SHEAR_CONFIG_CLINE_ARGS,
                                       d_defaults = D_EST_SHEAR_CONFIG_DEFAULTS, d_types = D_EST_SHEAR_CONFIG_TYPES,
                                       parsed_args = args)

    # check if profiling is to be enabled from the pipeline config
    profiling = args.pipeline_config[GlobalConfigKeys.PIP_PROFILE]

    if profiling:
        import cProfile
        logger.info("Profiling enabled")

        filename = os.path.join(args.workdir, args.logdir, "estimate_shears.prof")
        logger.info("Writing profiling data to %s", filename)

        cProfile.runctx("estimate_shears_from_args(args)", {},
                        {"estimate_shears_from_args": estimate_shears_from_args,
                         "args"                     : args,
                         "dry_run"                  : dry_run},
                        filename = filename)
    else:
        logger.info("Profiling disabled")
        estimate_shears_from_args(args, dry_run)

    logger.debug('# Exiting SHE_CTE_EstimateShear mainMethod() successfully.')


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
