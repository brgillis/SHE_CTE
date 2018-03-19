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
import pdb
from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_PPT.noise import get_var_ADU_per_pixel
from SHE_PPT.shear_utility import get_g_from_e
from SHE_PPT.logging import getLogger
from SHE_PPT.magic_values import scale_label, stamp_size_label
from SHE_PPT.she_image import SHEImage
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.table_formats.psf import tf as pstf
from SHE_PPT.table_formats.bfd_moments import initialise_bfd_moments_table, tf as setf

from SHE_BFD_CalculateMoments import bfd
from SHE_BFD_CalculateMoments.measure_moments import get_BFD_info



def bfd_measure_moments( data_stack, method_data=None ):

    logger = getLogger(mv.logger_name)
    logger.debug("Entering measuring BFD moments")

    # check to see if method_data has been given, if None go get 
    # standard file in SHE_BFD auxdir

    if method_data is None:
        method_data = bfd_load_method_data('/home/user/Work/Projects/SHE_BFD/SHE_BFD_CalculateMoments/auxdir/bfd.config')
    
    # Get lists of exposures, PSF images, and detection tables
    data_images = []
    detection_tables = []
    bulge_psf_images = []
    disk_psf_images = []
    psf_tables = []
    
    num_exposures = len(data_stack.exposures)

    for i in range(num_exposures):
        data_images.append(data_stack.exposures[i].science_image)
        detection_tables.append(data_stack.exposures[i].detections_table)
        bulge_psf_images.append(data_stack.exposures[i].bpsf_image)
        disk_psf_images.append(data_stack.exposures[i].dpsf_image)
        psf_tables.append(data_stack.exposures[i].psf_table)
 
    # Calculate the sky variances
    sky_vars = []
    for table_index in range(num_exposures):
        detections_table = detection_tables[table_index]
        sky_vars.append(get_var_ADU_per_pixel(pixel_value_ADU=0.,
                                                sky_level_ADU_per_sq_arcsec=detections_table.meta[detf.m.subtracted_sky_level],
                                                read_noise_count=detections_table.meta[detf.m.read_noise],
                                                pixel_scale=data_images[table_index].header[scale_label],
                                                gain=detections_table.meta[detf.m.gain]))


    # Get all unique IDs
    IDs = None
    for table_index in range(num_exposures):

        if IDs is None:
            IDs = set(detection_tables[table_index][detf.ID])
        else:
            IDs = set.union(IDs,detection_tables[table_index][detf.ID])

        # Set the ID as an index for each table
        detection_tables[table_index].add_index(detf.ID)
        psf_tables[table_index].add_index(pstf.ID)

    if method_data['isTarget'] == True:
        bfd_moments_table = initialise_bfd_moments_table(detection_tables[0],
                                                              optional_columns= \
                                                              [setf.x,
                                                               setf.y,
                                                               setf.bfd_moments,
                                                               setf.bfd_cov_even,
                                                               setf.bfd_cov_odd,
                                                               setf.bfd_pqr])

        bfd_moments_table.meta['WT_N'] = method_data['WEIGHT']['N']
        bfd_moments_table.meta['WT_SIGMA'] = method_data['WEIGHT']['SIGMA']
        nlost=0

    else:
        bfd_moments_table = initialise_bfd_moments_table(detection_tables[0],
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

        bfd_moments_table.meta['WT_N'] = method_data['WEIGHT']['N']
        bfd_moments_table.meta['WT_SIGMA'] = method_data['WEIGHT']['SIGMA']    
        # TO DO ADD TEMPLATE SPECIFIC STUFF
        #bfd_moments_table.meta['TMPL_SNMIN'] = method_data['TEMPLATE']['SNMIN']
        #bfd_moments_table.meta['TMPL_SIGMA_XY'] = method_data['TEMPLATE']['SIGMA_XY']
        #bfd_moments_table.meta['TMPL_SIGMA_FLUX'] = method_data['TEMPLATE']['SIGMA_FLUX']
        #bfd_moments_table.meta['TMPL_SIGMA_STEP'] = method_data['TEMPLATE']['SIGMA_STEP']
        #bfd_moments_table.meta['TMPL_SIGMA_MAX'] = method_data['TEMPLATE']['SIGMA_MAX']
        #bfd_moments_table.meta['TMPL_XY_MAX'] = method_data['TEMPLATE']['XY_MAX']

    for ID in IDs:
        
        x         = 0
        y         = 0
        moments   = 0
        cov_even  = 0
        cov_odd   = 0
        if method_data['isTarget'] == False:
            dm_dg1     = 0
            dm_dg2     = 0
            dm_dmu     = 0
            d2m_dg1dg1 = 0
            d2m_dg1dg2 = 0
            d2m_dg2dg2 = 0
            d2m_dg1dmu = 0
            d2m_dg2dmu = 0
            d2m_dmudmu = 0
            jsuppress  = 0
            weight     = 0

        for table_index in range(num_exposures):
            
            gal_stamps=[]
            bulge_psf_stamps=[]
            disk_psf_stamps=[]
            sky_var_list=[]

            try:
                
                # Get the row for this ID
                g_row = detection_tables[table_index].loc[ID]
                p_row = psf_tables[table_index].loc[ID]

                # Get galaxy and PSF stamps
                gal_stamps.append(data_images[table_index].extract_stamp(g_row[detf.gal_x],
                                                                         g_row[detf.gal_y],
                                                                         data_images[table_index].header[stamp_size_label],
                                                                         keep_header=True))

                bulge_psf_stamps.append(bulge_psf_images[table_index].extract_stamp(p_row[pstf.psf_x],
                                                                                    p_row[pstf.psf_y],
                                                                                    bulge_psf_images[table_index].header[stamp_size_label],
                                                                                    keep_header=True))

                disk_psf_stamps.append(disk_psf_images[table_index].extract_stamp(p_row[pstf.psf_x],
                                                                                  p_row[pstf.psf_y],
                                                                                  disk_psf_images[table_index].header[stamp_size_label],
                                                                                  keep_header=True))


                sky_var_list.append(sky_vars[table_index])

            except KeyError as e:
                if "No matches found for key" in e:
                    pass # ID isn't present in this table, so just skip it
                else:
                    raise
        
        # will try to stack stamps eventually ,for now use 1st
        #stacked_stamp = make_stacked_stamp(gal_stamps,psf_stamps,sky_var_list)

        bfd_info = get_bfd_info(gal_stamps[0], bulge_psf_stamps[0], sky_var_list[0], method_data)
        if method_data['isTarget']==True:
            if bfd_info.lost:
                nlost+=1
            else:
                cov_even=[]
                for i in xrange(5):
                    cov_even=np.append(cov_even,bfd_info.cov_even[i,i:5])

                cov_odd=np.append(bfd_info.cov_odd[0,0:2],bfd_info.cov_odd[1,1:2])

                # Add this row to the estimates table
                bfd_moments_table.add_row({ setf.ID : ID,
                                            setf.bfd_moments :bfd_info.m,
                                            setf.x : bfd_info.x,
                                            setf.y : bfd_info.y,
                                            setf.bfd_cov_even : cov_even,
                                            setf.bfd_cov_odd : cov_odd})

        else:
            if bfd_info.template_lost==False:
                for i, tmpl in enumerate(bfd_info.templates):
                    bfd_info.get_template(i)
                    bfd_moments_table.add_row({setf.ID:ID,
                                               setf.x:bfd_info.x,
                                               setf.y:bfd_info.y,
                                               setf.bfd_moments:bfd_info.m,
                                               setf.bfd_dm_dg1:bfd_info.dm_dg1,
                                               setf.bfd_dm_dg2:bfd_info.dm_dg2,
                                               setf.bfd_dm_dmu:bfd_info.dm_dmu,
                                               setf.bfd_d2m_dg1_dg1:bfd_info.d2m_dg1dg1,
                                               setf.bfd_d2m_dg1_dg2:bfd_info.d2m_dg1dg2,
                                               setf.bfd_d2m_dg2_dg2:bfd_info.d2m_dg2dg2,
                                               setf.bfd_d2m_dg1_dmu:bfd_info.d2m_dg1dmu,
                                               setf.bfd_d2m_dg2_dmu:bfd_info.d2m_dg2dmu,
                                               setf.bfd_d2m_dmu_dmu:bfd_info.d2m_dmudmu,
                                               setf.bfd_jsuppression:bfd_info.jsuppression,
                                               setf.bfd_weight:bfd_info.weight})
        if method_data['isTarget'] == True:
            bfd_moments_table.meta['NLOST'] = nlost

    logger.debug("Exiting BFD_measure_moments")

    return bfd_moments_table # No MCMC chains for this method

def bfd_load_method_data(method_data_filename):
    
    # set up dictionary structure with defaults
    method_data={'isTarget': True,
                 'WEIGHT':{'N':4,
                           'SIGMA':3.5},
                 'TEMPLATE':{}
             }

    # read in configuration file
    lines = [line.rstrip('\n') for line in open(method_data_filename)]

    # loop through configuration file to add parameters
    for l in lines:
        param, val = l.split(" ")

        if param.find('isTarget') != -1:
            if val=='True':
                method_data['isTarget']==True
            elif val=='False':
                method_data['isTarget']==False
            else:
                raise Exception("isTarget must be True or False")

        if param.find('WT') != -1:
            method_data['WEIGHT'][param.split("_")[1]]=np.float(val)
        
        if param.find('TMPL') != -1:
            raise Exception("Need to setup template info")


    return method_data
