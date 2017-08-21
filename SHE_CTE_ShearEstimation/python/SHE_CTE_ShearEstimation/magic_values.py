""" @file magic_values.py

    Created 27 Mar 2017

    Magic values for estimating shears on simulation data.

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

from SHE_GST_GalaxyImageGeneration import magic_values as sim_mv

logger_name = "SHE_CTE_EstimateShear"

image_tail = ".fits"
output_tail = "_shear_measurements.fits"

default_stamp_size_px = 256

fits_table_ID_label = sim_mv.detections_table_ID_label
fits_table_gal_x_label = "GAL_X"
fits_table_gal_y_label = "GAL_Y"
fits_table_gal_g1_label = "GAL_EST_G1"
fits_table_gal_g2_label = "GAL_EST_G2"
fits_table_gal_gerr_label = "GAL_GERR"
