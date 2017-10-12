""" @file analysis_pipeline.py

    Created 29 Aug 2017

    Pipeline script for the shear-estimation-only pipeline.

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

from SHE_CTE_Pipeline.package_definition import she_estimate_shear, she_validate_shear

@pipeline(outputs=('validated_shear_estimates_table'))
def shear_analysis_pipeline( data_images,
                             psf_calibration_product,
                             psf_calibration_listfile,
                             segmentation_images,
                             detections_tables,
                             astrometry_product,
                             aocs_time_series_product,
                             mission_time_product,
                             galaxy_population_priors_table,
                             calibration_parameters_product,
                             shear_validation_statistics_table ):
    
    psf_images_and_tables = she_fit_psf( data_images = data_images, 
                                         detections_tables = detections_tables,
                                         astrometry_product = astrometry_product,
                                         aocs_time_series_product = aocs_time_series_product,
                                         mission_time_product = mission_time_product,
                                         psf_calibration_product = psf_calibration_product,
                                         psf_calibration_listfile = psf_calibration_listfile, )
    
    shear_estimates_product, shear_estimates_listfile = she_estimate_shear( data_images = data_images,
                                                              psf_images_and_tables = psf_images_and_tables,
                                                              segmentation_images = segmentation_images,
                                                              detections_tables = detections_tables,
                                                              galaxy_population_priors_table = galaxy_population_priors_table,
                                                              calibration_parameters_product = calibration_parameters_product )
    
    validated_shear_estimates_table = she_validate_shear( shear_estimates_product = shear_estimates_product,
                                                          shear_estimates_listfile = shear_estimates_listfile,
                                                          validation_statistics_product = validation_statistics_product )
    
    return validated_shear_estimates_table

if __name__ == '__main__':
    from euclidwf.framework.graph_builder import build_graph
    from euclidwf.utilities import visualizer
    pydron_graph=build_graph(shear_analysis_pipeline)
    visualizer.visualize_graph(pydron_graph)