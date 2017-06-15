""" @file output_shear_estimates.py

    Created 27 Mar 2017

    Function to output shear estimates.

    ---------------------------------------------------------------------

    Copyright (C) 2017 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from astropy.table import Table
from SHE_CTE_ShearEstimation import magic_values as mv

def output_shear_estimates(stamps, output_file_name, galaxies_image_file_name):
    """
    @brief
        Outputs shear estimates into the desired fits table.
        
    @param stamps <list> List of stamp objects containing shear estimates
    
    @param kwargs <dict> Parsed command-line arguments
    
    @return None
    """
    
    # Get the desired output filename
    if output_file_name is None:
        # Replace the image tail with output tail
        output_file_name = galaxies_image_file_name.replace(mv.image_tail,mv.output_tail)
        
        # If that failed, make sure we don't overwrite the image
        if output_file_name == galaxies_image_file_name:
            output_file_name = galaxies_image_file_name + mv.output_tail
    
    # Initialize a table for output        
    otable = Table(names=[mv.fits_table_ID_label,
                          mv.fits_table_gal_x_label,
                          mv.fits_table_gal_y_label,
                          mv.fits_table_gal_g1_label,
                          mv.fits_table_gal_g2_label,
                          mv.fits_table_gal_gerr_label,],
                   dtype=[int,
                          float,
                          float,
                          float,
                          float,
                          float,])
    
    # Add each stamp's data to it in turn
    for stamp in stamps:
        otable.add_row({mv.fits_table_ID_label       : stamp.ID,
                        mv.fits_table_gal_x_label    : stamp.center.x,
                        mv.fits_table_gal_y_label    : stamp.center.y,
                        mv.fits_table_gal_g1_label   : stamp.shear_estimate.g1,
                        mv.fits_table_gal_g2_label   : stamp.shear_estimate.g2,
                        mv.fits_table_gal_gerr_label : stamp.shear_estimate.gerr,})
            
    # Output the table
    otable.write(output_file_name,format='fits',overwrite=True)
    
    return