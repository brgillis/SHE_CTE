""" @file magic_values.py

    Created 7 Apr 2017

    Magic values for measuring bias of shear estimates.

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
from SHE_CTE_ShearEstimation import magic_values as est_mv

logger_name = "SHE_CTE_MeasureBias"

default_output_filename = "shear_biases.fits"
default_output_format = "fits"

fits_table_ID_label = sim_mv.detections_table_ID_label
fits_table_sim_g1_label = "GAL_SIM_G1"
fits_table_sim_g2_label = "GAL_SIM_G2"
fits_table_est_g1_label = est_mv.fits_table_gal_g1_label
fits_table_est_g2_label = est_mv.fits_table_gal_g2_label
fits_table_est_gerr_label = est_mv.fits_table_gal_gerr_label
fits_table_est_g1_err_label = "GAL_G1_ERR"
fits_table_est_g2_err_label = "GAL_G2_ERR"
fits_table_est_e1_err_label = "GAL_ERR_E1"
fits_table_est_e2_err_label = "GAL_ERR_E2"

details_table_g_label = 'shear_magnitude'
details_table_beta_label = 'shear_angle'

var_e = {"m2": 0.30241731263684996,
         "m1": 0.27099491059570091,
         "p0": 0.2422044791810781,
         "p1": 0.21333834512347458,
         "p2": 0.18556590758327751}
