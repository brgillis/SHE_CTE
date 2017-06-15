""" @file extract_galaxy_stamps.py

    Created 27 Mar 2017

    Function to extract galaxy or psf stamps from an image.

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

import galsim

from SHE_SIM_galaxy_image_generation import magic_values as sim_mv

from SHE_CTE_ShearEstimation import magic_values as mv

class Stamp(object):
    pass

def extract_stamps(detections_table, image_hdulist,
                   x_colname, y_colname):
    
    image_hdu = image_hdulist[0]
    
    # Get the stamp size from the header
    try:
        stamp_size_px = image_hdu.header["STAMP_PX"]
    except KeyError as _e:
        stamp_size_px = mv.default_stamp_size_px
    
    # Get the image shape
    im_nx = image_hdu.header["NAXIS1"]
    im_ny = image_hdu.header["NAXIS2"]
    
    # Get the image
    image = galsim.fits.read(hdu_list=image_hdulist)
    
    # Get the x and y position columns
    IDs = detections_table[sim_mv.detections_table_ID_label]
    xcs = detections_table[x_colname]
    ycs = detections_table[y_colname]
    
    stamps = []
    for i in range(len(IDs)):
        
        # Get the initial bounds
        xmin = xcs[i] - stamp_size_px // 2
        xmax = xmin + stamp_size_px - 1
        ymin = ycs[i] - stamp_size_px // 2
        ymax = ymin + stamp_size_px - 1
        
        # Ensure the bounds are within the image bounds
        if xmin < 1:
            x_shift = 1-xmin
            xmin += x_shift
            xmax += x_shift
        if xmax > im_nx:
            x_shift = xmax - im_nx
            xmin -= x_shift
            xmax -= x_shift
        if ymin < 1:
            y_shift = 1-ymin
            ymin += y_shift
            ymax += y_shift
        if ymax > im_ny:
            y_shift = ymax - im_ny
            ymin -= y_shift
            ymax -= y_shift
            
        # Set up bounds for the stamp
        stamp_bounds = galsim.BoundsI(xmin=xmin,
                                      xmax=xmax,
                                      ymin=ymin,
                                      ymax=ymax,)
        
        # Get a subimage from these bounds
        subimage = image.subImage(stamp_bounds)
        
        # Set this up, along with its center, as an object to output
        stamp = Stamp()
        stamp.ID = IDs[i]
        stamp.image = subimage
        stamp.center = stamp_bounds.center()
        
        stamps.append(stamp)
            
    return stamps
    