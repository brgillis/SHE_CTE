""" @file package_definition.py

    Created 29 Aug 2017

    Package definition for the OU-SHE pipeline.

    ---------------------------------------------------------------------

    Copyright (C) 2012-2020 Euclid Science Ground Segment      
       
    This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General    
    Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)    
    any later version.    
       
    This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied    
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more    
    details.    
       
    You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to    
    the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
"""

from euclidwf.framework.taskdefs import Executable, Input, Output, ComputingResources

she_prepare_configs = Executable(command="E-Run SHE_CTE_PrepareConfigs",
                                 inputs=[Input("config_template"),Input("prep_config")],
                                 outputs=[Output("configs_list", mime_type="json", content_type="listfile")])

she_simulate_images = Executable(command="E-Run SHE_CTE_SimulateImages",
                                 inputs=[Input("config")],
                                 outputs=[Output("data_image_list", mime_type="json", content_type="listfile"),
                                          Output("psf_image_list", mime_type="json", content_type="listfile"),
                                          Output("noise_image_list", mime_type="json", content_type="listfile"),
                                          Output("mask_image_list", mime_type="json", content_type="listfile"),
                                          Output("segmentation_image_list", mime_type="json", content_type="listfile"),
                                          Output("detections_table_list", mime_type="json", content_type="listfile"),
                                          Output("details_table")])

she_estimate_shear = Executable(command="E-Run SHE_CTE_EstimateShear",
                                 inputs=[Input("data_images"),
                                         Input("psf_images_and_tables"),
                                         Input("segmentation_images"),
                                         Input("detections_tables"),
                                         Input("galaxy_population_priors_table"),
                                         Input("calibration_parameters_product")],
                                 outputs=[Output("shear_estimates_product")])

she_validate_shear = Executable(command="E-Run SHE_CTE_ValidateShear",
                                inputs=[Input("shear_estimates_product"),
                                        Input("shear_validation_statistics_table")],
                                outputs=[Output("validated_shear_estimates_table")])

she_measure_statistics = Executable(command="E-Run SHE_CTE_MeasureStatistics",
                                 inputs=[Input("details_table"),
                                         Input("shear_measurements_product")],
                                 outputs=[Output("statistics_product")])

she_measure_bias = Executable(command="E-Run SHE_CTE_MeasureBias",
                                 inputs=[Input("statistics_product_list",mime_type="json",content_type="listfile")],
                                 outputs=[Output("bias_measurements")])
