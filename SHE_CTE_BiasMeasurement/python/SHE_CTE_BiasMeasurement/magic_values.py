""" @file magic_values.py

    Created 7 Apr 2017

    Magic values for measuring bias of shear estimates.
"""

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

from SHE_GST_GalaxyImageGeneration import magic_values as sim_mv

from SHE_PPT.shear_estimates_format import *

logger_name = "SHE_CTE_MeasureBias"

fits_table_sim_g1_label = "SIM_G1"
fits_table_sim_g2_label = "SIM_G2"

default_output_filename = "shear_biases.fits"
default_output_format = "fits"

var_e = {"m2": 0.30241731263684996,
         "m1": 0.27099491059570091,
         "p0": 0.2422044791810781,
         "p1": 0.21333834512347458,
         "p2": 0.18556590758327751}
