""" @file bfd_measure_moments.py

    Created September 7, 2017

    Provides functions to measure the BFD moments of galaxies
"""

__updated__ = "2019-08-26"

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

from math import sqrt
import os
import subprocess

from SHE_PPT.file_io import read_xml_product, find_file, get_data_filename
from SHE_PPT.logging import getLogger
from SHE_PPT.magic_values import scale_label, stamp_size_label
from SHE_PPT.noise import get_var_ADU_per_pixel
from SHE_PPT.pipeline_utility import get_conditional_product
from SHE_PPT.she_image import SHEImage
from SHE_PPT.table_formats.bfd_moments import initialise_bfd_moments_table, tf as setf
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.table_formats.psf import tf as pstf
from SHE_PPT.table_utility import is_in_format
from astropy.table import Table
import pdb

from ElementsKernel.Auxiliary import getAuxiliaryPath, getAuxiliaryLocations
from SHE_BFD_CalculateMoments import bfd
from SHE_BFD_CalculateMoments.return_moments import get_bfd_info, load_bfd_configuration
from SHE_CTE_ShearEstimation import magic_values as mv
import numpy as np


#bfd = None
#get_bfd_info = None
#load_bfd_configuration = None
stamp_size = 128  # hardcoded for now
x_buffer = -5
y_buffer = -5


def bfd_measure_moments(data_stack, training_data, calibration_data, workdir, debug=False, *args, **kwargs):
    # not using training data or calibration data yet

    logger = getLogger(__name__)
    logger.info("Entering measuring BFD moments")

    # get configuration info
    config_data = load_bfd_configuration()
#    config_data = load_bfd_configuration(mode='training')
    # initialize the table
    if config_data['isTarget'] == True:
        bfd_moments_table = initialise_bfd_moments_table(optional_columns=[setf.bfd_moments,
                                                                           setf.bfd_cov_even,
                                                                           setf.bfd_cov_odd,
                                                                           setf.bfd_pqr])

        bfd_moments_table.meta['WT_N'] = config_data['WEIGHT']['N']
        bfd_moments_table.meta['WT_SIGMA'] = config_data['WEIGHT']['SIGMA']
        nlost = 0
    else:
        bfd_moments_table = initialise_bfd_moments_table(optional_columns=[setf.bfd_moments,
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

        bfd_moments_table.meta['WT_N'] = config_data['WEIGHT']['N']
        bfd_moments_table.meta['WT_SIGMA'] = config_data['WEIGHT']['SIGMA']
        bfd_moments_table.meta['T_SNMIN'] = config_data['TEMPLATE']['TMPL_SN_MIN']
        bfd_moments_table.meta['T_SIGXY'] = config_data['TEMPLATE']['TMPL_SIGMA_XY']
        bfd_moments_table.meta['T_SIGFLX'] = config_data['TEMPLATE']['TMPL_SIGMA_FLUX']
        bfd_moments_table.meta['T_SIGSTP'] = config_data['TEMPLATE']['TMPL_SIGMA_STEP']
        bfd_moments_table.meta['T_SIGMAX'] = config_data['TEMPLATE']['TMPL_SIGMA_MAX']
        bfd_moments_table.meta['T_XYMAX'] = config_data['TEMPLATE']['TMPL_XY_MAX']

    # define pixel scales
    if scale_label in data_stack.stacked_image.header:
        stacked_gal_scale = data_stack.stacked_image.header[scale_label] * 3600
    else:
        stacked_gal_scale = 0.1

    valid_detector = None
    for exposure in data_stack.exposures:
        for detector in exposure.detectors.ravel():
            if detector is not None:
                valid_detector = detector
    if valid_detector is None:
        raise RuntimeError("No valid detectors found")

    if scale_label in valid_detector.header:
        gal_scale = valid_detector.header[scale_label] * 3600
    else:
        gal_scale = 0.1

    if scale_label in data_stack.exposures[0].psf_data_hdulist[2].header:
        psf_scale = data_stack.exposures[0].psf_data_hdulist[2].header[scale_label] * 3600
    else:
        psf_scale = 0.02
    cnt = 0
    # Loop over galaxies and get an estimate for each one
    for row in data_stack.detections_catalogue:

        gal_id = row[detf.ID]
        gal_x_world = row[detf.gal_x_world]
        gal_y_world = row[detf.gal_y_world]

        # Get a stack of the galaxy images
        gal_stamp_stack = data_stack.extract_stamp_stack(x_world=gal_x_world,
                                                         y_world=gal_y_world,
                                                         width=stamp_size,
                                                         x_buffer=x_buffer,
                                                         y_buffer=y_buffer)

        # If there's no data for this galaxy, don't add it to the catalogue at all
        if gal_stamp_stack.is_empty():
            continue

        # Get stacks of the psf images
        bulge_psf_stack, disk_psf_stack = data_stack.extract_psf_stacks(gal_id=gal_id,
                                                                        make_stacked_psf=True,)

        # pull out images, sky variance, and wcs
        stacked_gal_stamp = gal_stamp_stack.stacked_image
        stacked_bulge_psf_stamp = bulge_psf_stack.stacked_image
        stacked_disk_psf_stamp = disk_psf_stack.stacked_image

        stacked_gal_image = stacked_gal_stamp.data
        stacked_psf_image = stacked_bulge_psf_stamp.data
        # FIX? what are units?
        sky_var = float(np.square(stacked_gal_stamp.noisemap.transpose()).mean())  # Galsim requires single float here
        # setup WCS with world and pixel coord of fixed point (nominal gal center)
        uvgal = (0., 0.)  # define center as 0,0 (gal_x_world,gal_y_world)
        xygal = np.array(stacked_gal_stamp.world2pix(gal_x_world, gal_y_world))
        transformation_matrix = stacked_gal_stamp.get_pix2world_transformation(xygal[0], xygal[1]) * 3600.
        duv_dxy = np.array([[transformation_matrix[0, 0], transformation_matrix[0, 1]],
                            [transformation_matrix[1, 0], transformation_matrix[1, 1]]])

        wcs = bfd.WCS(duv_dxy, xyref=xygal, uvref=uvgal)

        bfd_info = get_bfd_info(stacked_gal_image,
                                stacked_psf_image,
                                stacked_gal_scale,
                                psf_scale,
                                sky_var,
                                gal_id,
                                wcs,
                                config_data,
                                seg_map=stacked_gal_stamp.segmentation_map,
                                flag_map=stacked_gal_stamp.mask,
                                bkgd_map=stacked_gal_stamp.background_map,
                                do_bkgd_sub=False)  # bkgd already sub in stack

        if config_data['isTarget'] == True:
            if bfd_info.lost:
                nlost += 1
            else:
                cov_even = []
                for i in range(5):
                    cov_even = np.append(cov_even, bfd_info.cov_even[i, i:5])

                cov_odd = np.append(bfd_info.cov_odd[0, 0:2], bfd_info.cov_odd[1, 1:2])
                xsave= bfd_info.x/3600.+gal_x_world
                ysave = bfd_info.y/3600.+gal_y_world
                # Add this row to the estimates table
                bfd_moments_table.add_row({setf.ID: gal_id,
                                           setf.bfd_moments: bfd_info.m,
                                           setf.x_world: xsave,
                                           setf.y_world: ysave,
                                           setf.bfd_cov_even: cov_even,
                                           setf.bfd_cov_odd: cov_odd})

        else:

            if (bfd_info.template_lost == False) & (cnt < 200):
                for i, tmpl in enumerate(bfd_info.templates):
                    bfd_info.get_template(i)
                    bfd_moments_table.add_row({setf.ID: gal_id,
                                               setf.x_world: 0.0,
                                               setf.y_world: 0.0,
                                               setf.bfd_moments: bfd_info.m,
                                               setf.bfd_deriv_moments_dg1: bfd_info.dm_dg1,
                                               setf.bfd_deriv_moments_dg2: bfd_info.dm_dg2,
                                               setf.bfd_deriv_moments_dmu: bfd_info.dm_dmu,
                                               setf.bfd_2ndderiv_moments_dg1dg1: bfd_info.d2m_dg1dg1,
                                               setf.bfd_2ndderiv_moments_dg1dg2: bfd_info.d2m_dg1dg2,
                                               setf.bfd_2ndderiv_moments_dg2dg2: bfd_info.d2m_dg2dg2,
                                               setf.bfd_2ndderiv_moments_dg1dmu: bfd_info.d2m_dg1dmu,
                                               setf.bfd_2ndderiv_moments_dg2dmu: bfd_info.d2m_dg2dmu,
                                               setf.bfd_2ndderiv_moments_dmudmu: bfd_info.d2m_dmudmu,
                                               setf.bfd_jsuppress: bfd_info.jsuppression,
                                               setf.bfd_template_weight: bfd_info.weight})
        cnt += 1
        if config_data['isTarget'] == True:
            bfd_moments_table.meta['NLOST'] = nlost

    logger.info("Exiting BFD_measure_moments")

    return bfd_moments_table  # No MCMC chains for this method


def bfd_perform_integration(target_file=None, template_file=None):

    logger = getLogger(mv.logger_name)
    logger.info("Entering BFD integration")
    if target_file is None:
        logger.warn("no target file given")

    if template_file is None:
        template_file = getAuxiliaryPath("SHE_BFD_CalculateMoments/templateall.fits")

    sn1 = 5
    sn2 = 35
    sn3 = 100
    nfactor1 = 1.
    nfactor2 = 10.
    nthreads = 1
    call = ["SHE_BFD_BoostTest", "--targetFile", target_file, "--templateFile", template_file, "--selectSN", str(sn1) + "," + str(
        sn2) + "," + str(sn3), "--noiseFactor", str(nfactor1) + "," + str(nfactor2), "--nThreads", str(nthreads), "--seed", str(12345)]
    returncode = subprocess.run(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = np.str(returncode.stdout)
    stderr = np.str(returncode.stderr)
    stdoutsplit = stdout.split("\\n")
    stderrsplit = stderr.split("\\n")
    logger.info("BFD Integration STDOUT INFO:")
    for out in stdoutsplit:
        logger.info(out)
    logger.info("BFD Integration STDERR INFO:")
    for err in stderrsplit:
        logger.info(err)

    logger.debug("Exiting BFD integration")

    return


def bfd_load_training_data(training_data_filename=None, workdir="."):

    logger = getLogger(__name__)
    logger.debug("Entering load_BFD_Training_data")
    # basically just a check that the training data is in the correct format
    # and then returns the filepath+name for later use

    training_data = get_conditional_product(training_data_filename, workdir=workdir)

    if training_data is None:
        logger.warn("No training data provided; using default?")
        template_filename = getAuxiliaryPath("SHE_BFD_CalculateMoments/templateall.fits")
    else:
        # Read in the training data product
        template_filename = training_data.get_filename()

    t = Table.read(template_filename)
    if is_in_format(t, setf) == False:
        logger.warn("template file is not in correct format")

    logger.debug("Exiting bfd load training data")

    return template_filename
