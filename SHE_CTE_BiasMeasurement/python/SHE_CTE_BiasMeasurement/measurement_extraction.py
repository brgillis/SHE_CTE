""" @file measurement_extraction.py

    Created 10 Apr 2017

    Functions to get shear measurements and actual values from tables.
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

from astropy.table import Table, join, vstack
import numpy as np

from SHE_PPT.details_table_format import tf as datf
from SHE_PPT.shear_estimates_table_format import tf as setf
from SHE_PPT.table_utility import is_in_format

from SHE_CTE_BiasMeasurement import magic_values as mv

def isolate_measurements(measurements_table, var_e = mv.var_e["p0"]):
        
    # Trim the joined table to only have the needed columns
    colnames = measurements_table.colnames
    for colname in colnames:
        if colname not in [mv.shear_estimates_ID_label,
                           mv.shear_estimates_gal_g1_label,
                           mv.shear_estimates_gal_g2_label,
                           mv.shear_estimates_gal_g1_err_label,
                           mv.shear_estimates_gal_g2_err_label,
                           mv.shear_estimates_gal_e1_err_label,
                           mv.shear_estimates_gal_e2_err_label,]:
            measurements_table.remove_column(colname)
        else:
            # Handle specific columns as necessary
            if colname == mv.shear_estimates_gal_e1_err_label:
                # Use to calculate g1_err
                measurements_table[mv.shear_estimates_gal_g1_err_label] = np.sqrt(var_e+measurements_table[colname]**2)
                measurements_table.remove_column(colname)
            elif colname == mv.shear_estimates_gal_e2_err_label:
                # Use to calculate g2_err
                measurements_table[mv.shear_estimates_gal_g2_err_label] = np.sqrt(var_e+measurements_table[colname]**2)
                measurements_table.remove_column(colname)
    

def get_all_shear_measurements(input_files, var_e = mv.var_e["p0"]):
    """
    @brief
        Extract shear measurements and actual values from the input files.
    @param input_files
        <list<tuple>> List of tuples of detections and details files.
    @return <astropy.table.Table>
        Table of shear measurements.
    """
    
    joined_tables = []
    
    for measurements_table_filename, details_table_filename in input_files:
        
        # Load each table
        measurements_table = Table.read(measurements_table_filename)
        if not is_in_format(measurements_table,setf):
            raise Exception("Shear estimates table " + measurements_table_filename + " is in wrong format.")
        details_table = Table.read(details_table_filename)
        if not is_in_format(details_table,datf):
            raise Exception("Details table " + details_table + " is in wrong format.")
        
        # Get just the needed columns for the measurements table
        isolate_measurements(measurements_table, var_e)
        
        # Join the tables
        joined_table = join(measurements_table,
                            details_table,
                            keys=datf.ID)
        
        # Check that the tables were joined properly
        if len(joined_table) != len(measurements_table):
            raise Exception(measurements_table_filename + " could not be joined to " + details_table_filename + ".")
        
        # Trim the joined table to only have the needed columns
        
        colnames = joined_table.colnames
        for colname in colnames:
            if colname not in [datf.ID,
                               setf.gal_g1,
                               setf.gal_g2,
                               setf.gal_g1_err,
                               setf.gal_g2_err,
                               datf.shear_magnitude,
                               datf.shear_angle,]:
                joined_table.remove_column(colname)
                
        joined_table[mv.fits_table_sim_g1_label] = ( joined_table[datf.shear_magnitude] *
                                                     np.cos(np.pi/90.*joined_table[datf.shear_angle]) )
        joined_table[mv.fits_table_sim_g2_label] = ( joined_table[datf.shear_magnitude] *
                                                     np.sin(np.pi/90.*joined_table[datf.shear_angle]) )
        
        joined_table.remove_column(datf.shear_magnitude)
        joined_table.remove_column(datf.shear_angle)
        
        joined_tables.append(joined_table)
        
    # Append to the comparisons table
    comparison_table = vstack(joined_tables)
    
    return comparison_table