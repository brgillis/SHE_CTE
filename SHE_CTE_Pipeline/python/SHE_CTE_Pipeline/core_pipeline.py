""" @file core_pipeline.py

    Created 29 Aug 2017

    Pipeline script for the shear-measurement-only pipeline.

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

from SHE_CTE_Pipeline.package_definition import she_estimate_shear

@pipeline(outputs=('shear_measurements_table',))
def shear_measurement_pipeline( data_image, psf_image, noise_image, mask_image, detections_table ):
    
    shear_measurements_table = she_estimate_shear( data_image=data_image,
                                                   psf_image=psf_image,
                                                   noise_image=noise_image,
                                                   mask_image=mask_image,
                                                   detections_table=detections_table )
    
    return shear_measurements_table

if __name__ == '__main__':
    from euclidwf.framework.graph_builder import build_graph
    from euclidwf.utilities import visualizer
    pydron_graph=build_graph(correction_pipeline)
    visualizer.visualize_graph(pydron_graph)