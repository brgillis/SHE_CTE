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


from SHE_PPT.she_image import SHEImage
from SHE_PPT.magic_values import scale_label, stamp_size_label

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_BFD_CalculateMoments import bfd



def get_bfd_moments(gal_stamp, psf_stamp, sky_var, is_target):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_bfd_moments")
    
    # create wt function - eventually get from method_data
    wt = bfd.KSigmaWeight(**method_data['WEIGHT'])
    
    kdata = bfd.simpleImageTEST(gal_stamp,center,psf_stamp,\
                                pixel_scale=gal_stamp.header[scale_label],\
                                psf_pixel_scale=psf_stamp.header[scale_label],\
                                convolve_psf_with_pixel=True,\
                                pixel_noise=np.sqrt(sky_var))

    mc = bfd.MomentCalculator(kdata,wt)
#    if is_target == False:
#        t = mc.make_templates(**method_data['TEMPLATE'])
#        xyshift=None
#    else:
#        # for targets get moments at MX=MY=0
#        xyshift,error,msg=mc.recenter()

    logger.debug("Exiting get_bfd_moments")
    
    return mc



def measure_moments( data_stack, method_data ):

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
    
    if method_data['isTarget'] == True:
        bfd_moments_table = initialize_shear_estimation_table(detections_tables[0],
                                                              optional_columns= \
                                                              [setf.x,
                                                               setf.y,
                                                               setf.bfd_moments,
                                                               setf.bfd_cov_even,
                                                               setf.bfd_cov_odd,
                                                               setf.bfd_pqr])
        bfd_moments_table.header['WT_N'] = method_data['WEIGHT']['N']
        bfd_moments_table.header['WT_SIGMA'] = method_data['WEIGHT']['SIGMA']
        nlost=0

    else:
        bfd_moments_table = initialize_shear_estimation_table(detections_tables[0],
                                                         optional_columns= \
                                                         [setf.x,
                                                          setf.y,
                                                          setf.bfd_moments,
                                                          setf.bfd_deriv_moments_dg1,
                                                          setf.bfd_deriv_moments_dg2,
                                                          setf.bfd_deriv_moments_dmu,
                                                          setf.bfd_2ndderiv_moments_dg1dg1,
                                                          setf.bfd_2ndderiv_moments_dg1dg2,
                                                          setf.bfd_2ndderiv_moments_dg2dg2,
                                                          setf.bfd_2ndderiv_moments_dg1dmu,
                                                          setf.bfd_2ndderiv_moments_dg2dmu,
                                                          setf.bfd_2ndderiv_moments_dmudmu,
                                                          setf.bfd_template_weight,
                                                          setf.bfd_jsuppress])

        bfd_moments_table.header['WT_N'] = method_data['WEIGHT']['N']
        bfd_moments_table.header['WT_SIGMA'] = method_data['WEIGHT']['SIGMA']    

        bfd_moments_table.header['TMPL_SNMIN'] = method_data['TEMPLATE']['SNMIN']
        bfd_moments_table.header['TMPL_SIGMA_XY'] = method_data['TEMPLATE']['SIGMA_XY']
        bfd_moments_table.header['TMPL_SIGMA_FLUX'] = method_data['TEMPLATE']['SIGMA_FLUX']
        bfd_moments_table.header['TMPL_SIGMA_STEP'] = method_data['TEMPLATE']['SIGMA_STEP']
        bfd_moments_table.header['TMPL_SIGMA_MAX'] = method_data['TEMPLATE']['SIGMA_MAX']
        bfd_moments_table.header['TMPL_XY_MAX'] = method_data['TEMPLATE']['XY_MAX']

    for ID in IDs:
        
        x = 0
        y = 0
        moments = 0

        if method_data['isTarget'] == False:
            dm_dg1=0
            dm_dg2=0
            dm_dmu=0
            d2m_dg1dg1=0
            d2m_dg1dg2=0
            d2m_dg2dg2=0
            d2m_dg1dmu=0
            d2m_dg2dmu=0
            d2m_dmudmu=0
            jsuppress=0
            weight=0

        for detections_table in detections_tables:
            
            gal_stamps=[]
            psf_stamps=[]
            sky_var_list=[]

            try:
                
                # Get the row for this ID
                row = detections_table.loc[ID]
            
                # Get galaxy and PSF stamps
                gal_stamps.append(data_image.extract_stamp(row[detf.gal_x],
                                                     row[detf.gal_y],
                                                     data_image.header[stamp_size_label]))
                psf_stamps.append(psf_image.extract_stamp(row[detf.psf_x],
                                                    row[detf.psf_y],
                                                    psf_image.header[stamp_size_label]))

                sky_var_list.append(sky_vars[table_index])
            stacked_stamp = make_stacked_stamp(gal_stamps,psf_stamps,sky_var_list)
                                
            bfd_moments = get_bfd_moments(stacked_stamp.galaxy, stacked_stamp.psf, stacked_stamp.sky_var, method_data)

            if method_data['isTarget'] == True:
                xyshift,error,msg=bfd_moments.recenter()
                if error:
                    nlost+=1
                else:
                    x=xyshift[0]
                    y=xyshift[1]
                    galmoment=bfd_moments.get_moment(0,0)
                    galcov = bfd_moments.get_covariance()
                    moments=np.append(galmoment.even,galmoment.odd)
                    cov_even=np.append(galcov[0][0,0:5],galcov[0][1,1:5],galcov[0][2,2:5],galcov[0][3,3:5],galcov[0][4,4:5])
                    cov_odd=np.append(galcov[1][0,0:2],galcov[1][1,1:2])
            else:
                t=mc.make_template(**method_data['TEMPLATES'])
                for tmpl in t:
                        moments.append(np.append(bfd_moments.even,bfd_moments.odd))
                        dm_dg1.append(np.append(bfd))
                        #etc
                
            except KeyError as e:
                if "No matches found for key" in e:
                    pass # ID isn't present in this table, so just skip it
                else:
                    raise
                
        
        
        # Add this row to the estimates table
        bfd_moments_table.add_row({ setf.ID : detections_table[detf.ID][i],
                                    setf.bfd_moments : moments,
                                    setf.x : x,
                                    setf.y : y,
                                    setf.bfd_cov_even : cov_even,
                                    setf.bfd_cov_odd : cov_odd})
        
        
    logger.debug("Exiting GS_estimate_shear")
    
    return bfd_moments_table, None # No MCMC chains for this method
    
