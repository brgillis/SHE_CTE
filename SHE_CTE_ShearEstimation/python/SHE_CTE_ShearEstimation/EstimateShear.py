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

from typing import Any, Dict, Tuple, Type, Union

from SHE_CTE.executor import CteLogOptions, SheCteExecutor
from SHE_CTE_ShearEstimation.estimate_shears import estimate_shears_from_args
from SHE_PPT.argument_parser import SheArgumentParser
from SHE_PPT.constants.classes import ShearEstimationMethods
from SHE_PPT.constants.config import D_GLOBAL_CONFIG_CLINE_ARGS, D_GLOBAL_CONFIG_DEFAULTS, D_GLOBAL_CONFIG_TYPES
from SHE_PPT.executor import ReadConfigArgs
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import (AnalysisConfigKeys, CalibrationConfigKeys, ConfigKeys, )

EXEC_NAME = "SHE_CTE_EstimateShear"

# Set up dicts for pipeline config defaults and types
D_EST_SHEAR_CONFIG_DEFAULTS: Dict[ConfigKeys, Any] = {
    **D_GLOBAL_CONFIG_DEFAULTS,
    AnalysisConfigKeys.ES_METHODS                   : list(ShearEstimationMethods),
    AnalysisConfigKeys.ES_CHAINS_METHOD             : ShearEstimationMethods.LENSMC,
    AnalysisConfigKeys.ES_MEMMAP_IMAGES             : False,
    AnalysisConfigKeys.LENSMC_STAMP_SIZE            : 384,
    AnalysisConfigKeys.LENSMC_X_BUFFER              : 3,
    AnalysisConfigKeys.LENSMC_Y_BUFFER              : 3,
    AnalysisConfigKeys.LENSMC_NO_MASK_DILATION      : False,
    AnalysisConfigKeys.LENSMC_HL_TO_EXP             : 0.15,
    AnalysisConfigKeys.LENSMC_N_BULGE               : 1,
    AnalysisConfigKeys.LENSMC_N_DISC                : 1,
    AnalysisConfigKeys.LENSMC_E_MAX                 : 0.99,
    AnalysisConfigKeys.LENSMC_RE_MAX                : 2.,
    AnalysisConfigKeys.LENSMC_DELTA_MAX             : 0.3,
    AnalysisConfigKeys.LENSMC_E_FLAG                : 0.98,
    AnalysisConfigKeys.LENSMC_RE_FLAG               : 1.8,
    AnalysisConfigKeys.LENSMC_DELTA_FLAG            : 0.28,
    AnalysisConfigKeys.LENSMC_DISC_ONLY             : False,
    AnalysisConfigKeys.LENSMC_PSF_OVERSAMPLING      : 5,
    AnalysisConfigKeys.LENSMC_SEED                  : -1,
    AnalysisConfigKeys.LENSMC_SHAPE_NOISE           : 0.25,
    AnalysisConfigKeys.LENSMC_RETURN_CHAINS         : False,
    AnalysisConfigKeys.LENSMC_FAST_MODE             : False,
    AnalysisConfigKeys.LENSMC_INCLUDE_VIS_UNDETECTED: False,
    AnalysisConfigKeys.LENSMC_MONITOR               : False,
    }

D_EST_SHEAR_CONFIG_TYPES: Dict[ConfigKeys, Union[Type, Tuple[Type, Type]]] = {
    **D_GLOBAL_CONFIG_TYPES,
    AnalysisConfigKeys.ES_METHODS                   : (list, ShearEstimationMethods),
    AnalysisConfigKeys.ES_CHAINS_METHOD             : ShearEstimationMethods,
    AnalysisConfigKeys.ES_MEMMAP_IMAGES             : bool,
    AnalysisConfigKeys.LENSMC_STAMP_SIZE            : int,
    AnalysisConfigKeys.LENSMC_X_BUFFER              : int,
    AnalysisConfigKeys.LENSMC_Y_BUFFER              : int,
    AnalysisConfigKeys.LENSMC_NO_MASK_DILATION      : bool,
    AnalysisConfigKeys.LENSMC_HL_TO_EXP             : float,
    AnalysisConfigKeys.LENSMC_N_BULGE               : float,
    AnalysisConfigKeys.LENSMC_N_DISC                : float,
    AnalysisConfigKeys.LENSMC_E_MAX                 : float,
    AnalysisConfigKeys.LENSMC_RE_MAX                : float,
    AnalysisConfigKeys.LENSMC_DELTA_MAX             : float,
    AnalysisConfigKeys.LENSMC_E_FLAG                : float,
    AnalysisConfigKeys.LENSMC_RE_FLAG               : float,
    AnalysisConfigKeys.LENSMC_DELTA_FLAG            : float,
    AnalysisConfigKeys.LENSMC_DISC_ONLY             : bool,
    AnalysisConfigKeys.LENSMC_PSF_OVERSAMPLING      : int,
    AnalysisConfigKeys.LENSMC_SEED                  : int,
    AnalysisConfigKeys.LENSMC_SHAPE_NOISE           : float,
    AnalysisConfigKeys.LENSMC_RETURN_CHAINS         : bool,
    AnalysisConfigKeys.LENSMC_FAST_MODE             : bool,
    AnalysisConfigKeys.LENSMC_INCLUDE_VIS_UNDETECTED: bool,
    AnalysisConfigKeys.LENSMC_MONITOR               : bool,
    }

D_EST_SHEAR_CONFIG_CLINE_ARGS: Dict[ConfigKeys, str] = {
    **D_GLOBAL_CONFIG_CLINE_ARGS,
    AnalysisConfigKeys.ES_METHODS                   : "methods",
    AnalysisConfigKeys.ES_CHAINS_METHOD             : "chains_method",
    AnalysisConfigKeys.ES_MEMMAP_IMAGES             : "memmap_images",
    AnalysisConfigKeys.LENSMC_STAMP_SIZE            : "lensmc_stamp_size",
    AnalysisConfigKeys.LENSMC_X_BUFFER              : "lensmc_x_buffer",
    AnalysisConfigKeys.LENSMC_Y_BUFFER              : "lensmc_y_buffer",
    AnalysisConfigKeys.LENSMC_NO_MASK_DILATION      : "lensmc_no_mask_dilation",
    AnalysisConfigKeys.LENSMC_HL_TO_EXP             : "lensmc_hl_to_exp",
    AnalysisConfigKeys.LENSMC_N_BULGE               : "lensmc_n_bulge",
    AnalysisConfigKeys.LENSMC_N_DISC                : "lensmc_n_disc",
    AnalysisConfigKeys.LENSMC_E_MAX                 : "lensmc_e_max",
    AnalysisConfigKeys.LENSMC_RE_MAX                : "lensmc_re_max",
    AnalysisConfigKeys.LENSMC_DELTA_MAX             : "lensmc_delta_max",
    AnalysisConfigKeys.LENSMC_E_FLAG                : "lensmc_e_flag",
    AnalysisConfigKeys.LENSMC_RE_FLAG               : "lensmc_re_flag",
    AnalysisConfigKeys.LENSMC_DELTA_FLAG            : "lensmc_delta_flag",
    AnalysisConfigKeys.LENSMC_DISC_ONLY             : "lensmc_disc_only",
    AnalysisConfigKeys.LENSMC_PSF_OVERSAMPLING      : "lensmc_psf_oversampling",
    AnalysisConfigKeys.LENSMC_SEED                  : "lensmc_seed",
    AnalysisConfigKeys.LENSMC_SHAPE_NOISE           : "lensmc_shape_noise",
    AnalysisConfigKeys.LENSMC_RETURN_CHAINS         : "lensmc_return_chains",
    AnalysisConfigKeys.LENSMC_FAST_MODE             : "lensmc_fast_mode",
    AnalysisConfigKeys.LENSMC_INCLUDE_VIS_UNDETECTED: "lensmc_include_vis_undetected",
    AnalysisConfigKeys.LENSMC_MONITOR               : "lensmc_monitor",
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
    logger.debug('# Entering SHE_CTE_EstimateShear defineSpecificProgramOptions()')
    logger.debug('#')

    parser = SheArgumentParser()

    # Required input arguments

    parser.add_input_arg('--data_images', type = str,
                         help = '.json listfile containing filenames of data image products.')

    parser.add_input_arg('--stacked_image', type = str,
                         help = 'Filename of the stacked image product')

    parser.add_input_arg('--psf_images_and_tables', type = str,
                         help = '.json listfile containing filenames of psf image products.')

    parser.add_input_arg('--segmentation_images', type = str,
                         help = '.json listfile containing filenames of segmentation map images.')

    parser.add_input_arg('--stacked_segmentation_image', type = str,
                         help = 'Filename of the stacked segmentation map image.')

    parser.add_input_arg('--detections_tables', type = str,
                         help = '.json listfile containing filenames of detections table products.')

    parser.add_input_arg('--object_ids', type = str,
                         help = 'XML dataproduct that contains within it a list of object ids to loop over')

    parser.add_input_arg('--ksb_training_data', type = str,
                         help = 'Data product for KSB training data.')

    parser.add_input_arg('--lensmc_training_data', type = str,
                         help = 'Data product for LensMC training data.')

    parser.add_input_arg('--momentsml_training_data', type = str,
                         help = 'Data product for MomentsML training data.')

    parser.add_input_arg('--regauss_training_data', type = str,
                         help = 'Data product for REGAUSS training data.')

    parser.add_input_arg('--mdb', type = str,
                         help = 'Mission Database .xml file')

    parser.add_input_arg('--galaxy_population_priors_table', type = str,
                         help = '.json listfile containing filenames of detections table products.')

    parser.add_input_arg('--calibration_parameters_product', type = str,
                         help = 'Filename of calibration parameters product (XML data product).')

    # Option arguments (cannot be used in pipeline)

    parser.add_option_arg('--debug', action = 'store_true',
                          help = 'Enables debug mode - only process first 1000 galaxies.')

    parser.add_option_arg('--methods', type = str, nargs = '*',
                          help = 'Which shear estimation methods to apply. If not specified, all will be run.')

    parser.add_option_arg('--chains_method', type = str,
                          help = 'Which shear estimation method to generate chains with.')

    parser.add_option_arg('--memmap_images', action = 'store_true',
                          help = 'Memory-map images rather than reading them on-demand.')

    # LensMC-specific option arguments

    parser.add_option_arg('--lensmc_stamp_size',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_STAMP_SIZE],
                          help = 'LensMC: Requested stamp size in pixels.')

    parser.add_option_arg('--lensmc_x_buffer',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_X_BUFFER],
                          help = 'LensMC: Do not fit object if closer to edge of detector by this number of pixels '
                                 'along the x coordinate.')

    parser.add_option_arg('--lensmc_y_buffer',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_Y_BUFFER],
                          help = 'LensMC: Do not fit object if closer to edge of detector by this number of pixels '
                                 'along the y coordinate.')

    parser.add_option_arg('--lensmc_no_mask_dilation',
                          action = 'store_true',
                          help = 'LensMC: Do not dilate mask by one pixel.')

    parser.add_option_arg('--lensmc_hl_to_exp',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_HL_TO_EXP],
                          help = 'LensMC: Half-light radius of the bulge to exponential scalelength of the disc.')

    parser.add_option_arg('--lensmc_n_bulge',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_N_BULGE],
                          help = 'LensMC: Bulge Sersic index; available: n=(1., 1.5, 2., 2.5, 3., 3.5, 4.).')

    parser.add_option_arg('--lensmc_n_disc',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_N_DISC],
                          help = 'LensMC: Disc Sersic index; available: n=1.')

    parser.add_option_arg('--lensmc_e_max',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_E_MAX],
                          help = 'LensMC: Hard upper bound on ellipticity.')

    parser.add_option_arg('--lensmc_re_max',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_RE_MAX],
                          help = 'LensMC: Hard upper bound on effective radius.')

    parser.add_option_arg('--lensmc_delta_max',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_DELTA_MAX],
                          help = 'LensMC: Hard upper bound on position offset.')

    parser.add_option_arg('--lensmc_e_flag',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_E_FLAG],
                          help = 'LensMC: Flagging threshold for ellipticity.')

    parser.add_option_arg('--lensmc_re_flag',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_RE_FLAG],
                          help = 'LensMC: Flagging threshold for effective radius.')

    parser.add_option_arg('--lensmc_delta_flag',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_DELTA_FLAG],
                          help = 'LensMC: Flagging threshold for position offset.')

    parser.add_option_arg('--lensmc_disc_only',
                          action = 'store_true',
                          help = 'LensMC: Whether to fit only for a disc component.')

    parser.add_option_arg('--lensmc_psf_oversampling',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_PSF_OVERSAMPLING],
                          help = 'LensMC: PSF oversampling factor.')

    parser.add_option_arg('--lensmc_seed',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_SEED],
                          help = 'LensMC: Seed the random sampler.')

    parser.add_option_arg('--lensmc_shape_noise',
                          type = D_EST_SHEAR_CONFIG_TYPES[AnalysisConfigKeys.LENSMC_SHAPE_NOISE],
                          help = 'LensMC: Shape noise standard deviation if not provided by training data.')

    parser.add_option_arg('--lensmc_return_chains',
                          action = 'store_true',
                          help = 'LensMC: Whether to return the chains.')

    parser.add_option_arg('--lensmc_fast_mode',
                          action = 'store_true',
                          help = 'Enable LensMC fast mode. Override any sampling settings and produce MAP estimate '
                                 '(without MCMC/errors/intcal).')

    parser.add_option_arg('--lensmc_include_vis_undetected',
                          action = 'store_true',
                          help = 'LensMC: Measure all objects, including those that have not been detected by VIS.')

    parser.add_option_arg('--lensmc_monitor',
                          action = 'store_true',
                          help = 'LensMC: Monitor data and visually check consistency of input data.')

    # Output arguments

    parser.add_output_arg('--shear_estimates_product', type = str,
                          help = 'XML data product to contain file links to the shear estimates tables.')

    parser.add_output_arg('--she_lensmc_chains', type = str,
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

    executor = SheCteExecutor(run_from_args_function = estimate_shears_from_args,
                              log_options = CteLogOptions(executable_name = EXEC_NAME),
                              config_args = ReadConfigArgs(d_config_defaults = D_EST_SHEAR_CONFIG_DEFAULTS,
                                                           d_config_types = D_EST_SHEAR_CONFIG_TYPES,
                                                           d_config_cline_args = D_EST_SHEAR_CONFIG_CLINE_ARGS,
                                                           s_config_keys_types = {AnalysisConfigKeys,
                                                                                  CalibrationConfigKeys},
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
