""" @file bfd_measure_moments.py

    Created September 7, 2017

    Provides functions to measure the BFD moments of galaxies

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

from math import sqrt
import numpy as np

from SHE_GST_IceBRGpy.logging import getLogger

from SHE_GST_GalaxyImageGeneration.noise import get_var_ADU_per_pixel

from SHE_PPT.detections_table_format import tf as detf
from SHE_PPT.bfd_moments_table_format import initialize_bfd_moments_table, tf as setf

from SHE_PPT.she_image import SHEImage
from SHE_PPT.magic_values import scale_label, stamp_size_label

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_BFD_CalculateMoments import bfd



def get_bfd_moments(gal_stamp, psf_stamp, sky_var, is_target):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_bfd_moments")
    
    # create wt function - eventually get from method_data
    wt = bfd.KSigmaWeight(n=4.0,sigma=3.5)
    
    kdata = bfd.simpleImageTEST(gal_stamp,center,psf_stamp,\
                                pixel_scale=gal_stamp.header[scale_label],\
                                psf_pixel_scale=psf_stamp.header[scale_label],\
                                convolve_psf_with_pixel=True,\
                                pixel_noise=np.sqrt(sky_var))

    mc = bfd.MomentCalculator(kdata,wt)
    if is_target == True:
        t = mc.make_templates(**params['TEMPLATE'])
        xyshift=None
    else:
        # for targets get moments at MX=MY=0
        xyshift,error,msg=mc.recenter()

    logger.debug("Exiting get_bfd_moments")
    
    return mc,xyshift



def measure_moments( data_stack, method, shape_noise_var ):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering measuring BFD moments")
    
    # Get lists of exposures, PSF images, and detection tables
    data_images = data_stack.get_data_images()
    psf_images = data_stack.get_psf_images()
    detection_tables = data_stack.get_detection_tables_images()
    
    # Calculate the sky variances
    sky_vars = []
    for detection_table in detection_tables:
        sky_vars.append(get_var_ADU_per_pixel(pixel_value_ADU=0.,
                                                sky_level_ADU_per_sq_arcsec=detections_table[detf.m.subtracted_sky_level],
                                                read_noise_count=detections_table[detf.m.read_noise],
                                                pixel_scale=data_image.scale,
                                                gain=detections_table[detf.m.gain]))
    
    num_tables = len(detections_tables)
    
    # Get all unique IDs
    IDs = None
    for table_index in range(num_tables):
        if IDs is None:
            IDs = set(detections_tables[table_index][detf.ID])
        else:
            IDs = set.union(IDs,detections_tables[table_index][detf.ID])
            
    # Set the ID as an index for each table
    for detection_table in detection_tables:
        detection_table.add_index(detf.ID)
    
    bfd_moments_table = initialize_bfd_moments_table(detections_tables[0],
                                                             optional_columns=[setf.e1_err,setf.e2_err])
    
    for ID in IDs:
        
        xys = []
        moments = []
        if is_target = False:
            moments_deriv=[]

        for detections_table in detections_tables:
            
            try:
                
                # Get the row for this ID
                row = detections_table.loc[ID]
            
                # Get galaxy and PSF stamps
                gal_stamp = data_image.extract_stamp(row[detf.gal_x],
                                                     row[detf.gal_y],
                                                     data_image.header[stamp_size_label])
                psf_stamp = psf_image.extract_stamp(row[detf.psf_x],
                                                    row[detf.psf_y],
                                                    psf_image.header[stamp_size_label])
        
                bfd_moments = get_bfd_moments(gal_stamp, psf_stamp, sky_vars[table_index],is_target)
                
                xys.append(bfd_moments[1])
                moments.append(bfd_moments[0].get_moment(0,0))
                if is_target == False:
                    moments_deriv.append(bfd_moments)
                
            except KeyError as e:
                if "No matches found for key" in e:
                    pass # ID isn't present in this table, so just skip it
                else:
                    raise
                
        g1s = np.array(g1s)
        g2s = np.array(g2s)
        gerrs = np.array(gerrs)
        
            
        # Add this row to the estimates table
        bfd_moments_table.add_row({ setf.ID : detections_table[detf.ID][i],
                                        setf.g1 : g1,
                                        setf.g2 : g2,
                                        setf.e1_err : gerr1,
                                        setf.e2_err : gerr2,
                                       })
        
    
    logger.debug("Exiting GS_estimate_shear")
    
    return bfd_moments_table, None # No MCMC chains for this method
    
