""" @file galsim_estimate_shear.py

    Created 27 Mar 2017

    Provides functions to measure the shape of a galaxy image.
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

__updated__ = "2020-07-10"

from copy import deepcopy
from math import sqrt

from SHE_PPT import flags
from SHE_PPT import mdb
from SHE_PPT.logging import getLogger
from SHE_PPT.magic_values import scale_label, gain_label
from SHE_PPT.she_image import SHEImage
from SHE_PPT.shear_utility import (get_g_from_e, correct_for_wcs_shear_and_rotation, check_data_quality)
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf
from SHE_PPT.table_formats.she_ksb_measurements import initialise_ksb_measurements_table, tf as ksbm_tf
from SHE_PPT.table_formats.she_regauss_measurements import initialise_regauss_measurements_table, tf as regm_tf
from SHE_PPT.utility import run_only_once
import galsim

import numpy as np

stamp_size = 256
x_buffer = -5
y_buffer = -5
get_exposure_estimates = False

default_galaxy_scale = 0.1 / 3600
default_psf_scale = 0.02 / 3600

initialisation_methods = {"KSB": initialise_ksb_measurements_table,
                          "REGAUSS": initialise_regauss_measurements_table}

table_formats = {"KSB": ksbm_tf,
                 "REGAUSS": regm_tf}


class ShearEstimate(object):

    def __init__(self, g1, g2, gerr, re, snr, x, y):
        self.g1 = g1
        self.g2 = g2
        self.gerr = gerr
        self.g1_err = gerr
        self.g2_err = gerr
        self.g1g2_covar = 0
        self.re = re
        self.snr = snr
        self.x = x
        self.y = y
        self.flags = 0


# Set up an instance of ShearEstimate which we'll use for any errors
error_shear_estimate = ShearEstimate(g1=np.NaN,
                                     g2=np.NaN,
                                     gerr=np.inf,
                                     re=np.NaN,
                                     snr=np.NaN,
                                     x=np.NaN,
                                     y=np.NaN)
error_shear_estimate.g1_err = np.inf
error_shear_estimate.g2_err = np.inf
error_shear_estimate.g1g1_covar = np.NaN

snr_cutoff = 15
downweight_error = 0.5
downweight_power = 4


def get_downweight_error(snr):
    return downweight_error / (1 + (snr / snr_cutoff) ** downweight_power)


@run_only_once
def log_no_galaxy_scale():
    logger = getLogger(__name__)
    logger.warning("Cannot find pixel scale in image header. Using default value of " + str(default_galaxy_scale))


@run_only_once
def log_no_psf_scale():
    logger = getLogger(__name__)
    logger.warning("Cannot find pixel scale in PSF header. Using default value of " + str(default_psf_scale))


def get_resampled_image(initial_image, resampled_scale, resampled_nx, resampled_ny, stacked=False):

    if scale_label in initial_image.header:
        in_scale = initial_image.header[scale_label]
    else:
        log_no_galaxy_scale()
        in_scale = default_galaxy_scale

    # Add a default background map if necessary
    initial_image.add_default_background_map()

    if stacked:
        bkg_subtracted_stamp_data = initial_image.data
    else:
        bkg_subtracted_stamp_data = initial_image.data - initial_image.background_map

    window_nx = int(resampled_nx * resampled_scale / in_scale) + 1
    window_ny = int(resampled_ny * resampled_scale / in_scale) + 1

    xm = (bkg_subtracted_stamp_data.shape[0] - window_nx) // 2 + 1
    if xm < 1:
        xm = 1
    xh = xm + window_nx - 1
    if xh > bkg_subtracted_stamp_data.shape[0]:
        xh = bkg_subtracted_stamp_data.shape[0]
    ym = (bkg_subtracted_stamp_data.shape[1] - window_ny) // 2 + 1
    if ym < 1:
        ym = 1
    yh = ym + window_ny - 1
    if yh > bkg_subtracted_stamp_data.shape[1]:
        yh = bkg_subtracted_stamp_data.shape[1]

    subimage = galsim.Image(bkg_subtracted_stamp_data, scale=in_scale).subImage(galsim.BoundsI(xm, xh, ym, yh))

    if (subimage.array == 0).all():
        resampled_image = SHEImage(np.zeros((resampled_nx, resampled_ny)))
    else:
        resampled_gs_image = galsim.Image(resampled_nx, resampled_ny, scale=resampled_scale)

        galsim.InterpolatedImage(subimage).drawImage(resampled_gs_image, use_true_center=True)

        resampled_image = SHEImage(resampled_gs_image.array)

    resampled_image.add_default_header()
    resampled_image.header[scale_label] = resampled_scale

    return resampled_image


def KSB_estimate_shear(data_stack, training_data, calibration_data, workdir, sim_sc4_fix=False, *args, **kwargs):
    """
        Top-level function to run the KSB shear estimation algorithm. 

        The estimation is performed through the implementation of the method in the GalSim toolkit, with a wrapper to
        handle interfacing and to project the galaxy and PSF onto the same pixel scale.

        Parameters
        ----------
        frame_stack : SHEFrameStack
            All the data (image, detections table, psfs, etc.)
        training_data : SHE_CTE.control_training_data.ControlTraining
            All the method parameters
        calibration_data : None
            Currently unused; needed for proper interface
        workdir : string
            Working directory
        sim_sc4_fix : bool
            If set to True, will not correct for astrometric shear (to account for this effect not being present in
            SIM's SC4 simulations)
        *args, **kwargs : Needed in anticipation of future interface changes that might not affect all methods.

        Returns
        -------
        astropy.table.Table
            A table containing one row for each galaxy, in the shear_estimates table format.
    """

    # Not using training or calibration data at this stage
    return GS_estimate_shear(data_stack, training_data=training_data, method="KSB", workdir=workdir,
                             sim_sc4_fix=sim_sc4_fix, *args, **kwargs)


def REGAUSS_estimate_shear(data_stack, training_data, calibration_data, workdir, sim_sc4_fix=False, *args, **kwargs):
    """
        Top-level function to run the REGAUSS shear estimation algorithm. 

        The estimation is performed through the implementation of the method in the GalSim toolkit, with a wrapper to
        handle interfacing and to project the galaxy and PSF onto the same pixel scale.

        Parameters
        ----------
        frame_stack : SHEFrameStack
            All the data (image, detections table, psfs, etc.)
        training_data : SHE_CTE.control_training_data.ControlTraining
            All the method parameters
        calibration_data : None
            Currently unused; needed for proper interface
        workdir : string
            Working directory
        sim_sc4_fix : bool
            If set to True, will not correct for astrometric shear (to account for this effect not being present in
            SIM's SC4 simulations)
        *args, **kwargs : Needed in anticipation of future interface changes that might not affect all methods.

        Returns
        -------
        astropy.table.Table
            A table containing one row for each galaxy, in the shear_estimates table format.
    """

    # Not using training or calibration data at this stage
    return GS_estimate_shear(data_stack, training_data=training_data, method="REGAUSS", workdir=workdir,
                             sim_sc4_fix=sim_sc4_fix, *args, **kwargs)


def get_KSB_shear_estimate(galsim_shear_estimate, scale):

    logger = getLogger(__name__)
    logger.debug("Entering get_KSB_shear_estimate")

    g1 = galsim_shear_estimate.corrected_g1
    g2 = galsim_shear_estimate.corrected_g2
    mag = g1 ** 2 + g2 ** 2

    shear_estimate = ShearEstimate(galsim_shear_estimate.corrected_g1,
                                   galsim_shear_estimate.corrected_g2,
                                   galsim_shear_estimate.corrected_shape_err,
                                   galsim_shear_estimate.moments_sigma * scale,
                                   -1,
                                   galsim_shear_estimate.moments_centroid.x,
                                   galsim_shear_estimate.moments_centroid.y)

    if mag > 1:

        # Flag an error if the shear is too large

        logger.warning("Magnitude of g shear is too large: " + str(mag))
        shear_estimate.flags |= flags.flag_too_large_shear

    logger.debug("Exiting get_KSB_shear_estimate")

    return shear_estimate


def get_REGAUSS_shear_estimate(galsim_shear_estimate, scale):

    logger = getLogger(__name__)
    logger.debug("Entering get_REGAUSS_shear_estimate")

    e1 = galsim_shear_estimate.corrected_e1
    e2 = galsim_shear_estimate.corrected_e2
    mag = e1 ** 2 + e2 ** 2

    shear_estimate = ShearEstimate(np.NaN, np.NaN, np.inf,
                                   galsim_shear_estimate.moments_sigma * scale,
                                   -1,
                                   galsim_shear_estimate.moments_centroid.x,
                                   galsim_shear_estimate.moments_centroid.y)

    if mag > 1:

        # Can't calculate a shear to return if the magnitude is greater than 1,
        # so flag an error

        logger.warning("Magnitude of g shear is too large: " + str(mag))
        shear_estimate.flags |= flags.flag_too_large_shear

    else:

        g1, g2 = get_g_from_e(e1, e2)
        gerr = galsim_shear_estimate.corrected_shape_err * np.sqrt((g1 ** 2 + g2 ** 2) / (e1 ** 2 + e2 ** 2))

        shear_estimate.g1 = g1
        shear_estimate.g2 = g2
        shear_estimate.gerr = gerr
        shear_estimate.g1_err = gerr
        shear_estimate.g2_err = gerr
        shear_estimate.g1g2_covar = 0

    logger.debug("Exiting get_REGAUSS_shear_estimate")

    return shear_estimate


def get_shear_estimate(gal_stamp, psf_stamp, gal_scale, psf_scale, ID, method, stacked=False,
                       sim_sc4_fix=False):

    logger = getLogger(__name__)
    logger.debug("Entering get_shear_estimate")

    # Check that there aren't any obvious issues with the data
    data_quality_flags = check_data_quality(gal_stamp, psf_stamp, stacked=stacked)

    # If we hit any failure flags, return now with an error
    if data_quality_flags & flags.failure_flags:

        shear_estimate = deepcopy(error_shear_estimate)
        shear_estimate.flags |= data_quality_flags

        return shear_estimate
    elif data_quality_flags:
        # For non-failure flags, overwrite with default noisemap and/or segmentation map
        if (data_quality_flags & flags.flag_no_noisemap) or (data_quality_flags & flags.flag_corrupt_noisemap):
            gal_stamp.add_default_noisemap(force=True)
        if (data_quality_flags & flags.flag_no_segmentation_map) or (data_quality_flags &
                                                                     flags.flag_corrupt_segmentation_map):
            gal_stamp.add_default_segmentation_map(force=True)

    # Subtract off the background from the galaxy
    if stacked:
        bkg_subtracted_gal_stamp_data = gal_stamp.data
    else:
        bkg_subtracted_gal_stamp_data = gal_stamp.data - gal_stamp.background_map

    # Estimate the size of the galaxy, so we can figure out how big we need to make the resampled stamp

    gal_sigs = [5.0, 2.0, 0.5, 10.0]

    for gal_sig in gal_sigs:

        try:

            gal_mom = galsim.hsm.FindAdaptiveMom(galsim.Image(bkg_subtracted_gal_stamp_data.transpose(), scale=gal_scale),
                                                 badpix=galsim.Image(
                                                     (gal_stamp.boolmask).astype(np.uint16).transpose(), scale=gal_scale),
                                                 guess_sig=gal_sig,)

            resampled_gal_stamp_size = int(20 * gal_mom.moments_sigma * gal_scale /
                                           psf_scale)    # Calculate the galaxy's S/N
            a_eff = np.pi * (3 * gal_mom.moments_sigma * np.sqrt(2 * np.log(2)))
            gain = mdb.get_mdb_value(mdb.mdb_keys.vis_gain)
            signal_to_noise = (gain * gal_mom.moments_amp / np.sqrt(gain * gal_mom.moments_amp + a_eff *
                                                                    (gain * np.square(gal_stamp.noisemap.transpose()).mean()) ** 2))
            break

        except RuntimeError as e:

            if "HSM Error" in str(e):

                # Flag an error if we're on the last guess_sigma we're trying
                if gal_sig == gal_sigs[-1]:

                    # The galaxy is probably too small in this case, so flag that
                    shear_estimate = deepcopy(error_shear_estimate)
                    shear_estimate.flags |= flags.flag_object_too_small | data_quality_flags

                    # Fill in the x/y with the rough position, at least
                    shear_estimate.x = gal_stamp.shape[0] / 2
                    shear_estimate.y = gal_stamp.shape[1] / 2

                    return shear_estimate
                else:
                    continue
            else:
                raise

    # Get a resampled galaxy stamp
    resampled_gal_stamp = get_resampled_image(gal_stamp, psf_scale, resampled_gal_stamp_size, resampled_gal_stamp_size,
                                              stacked=stacked)

    # Get a resampled badpix map
    supersampled_badpix = SHEImage((gal_stamp.boolmask).astype(float))
    supersampled_badpix.add_default_header()
    if scale_label in gal_stamp.header:
        supersampled_badpix.header[scale_label] = gal_stamp.header[scale_label]
    else:
        supersampled_badpix.header[scale_label] = default_galaxy_scale
    resampled_badpix = get_resampled_image(supersampled_badpix, psf_scale,
                                           resampled_gal_stamp_size, resampled_gal_stamp_size)

    badpix = (resampled_badpix.data > 0.5).astype(np.uint16)  # Galsim requires int array

    # FIXME - What units should sky_var be in?
    sky_var = float(np.square(gal_stamp.noisemap.transpose()).mean())  # Galsim requires single float here

    try:

        galsim_shear_estimate = galsim.hsm.EstimateShear(gal_image=galsim.Image(resampled_gal_stamp.data.transpose(), scale=psf_scale),
                                                         PSF_image=galsim.Image(psf_stamp.data.transpose(),
                                                                                scale=psf_scale),
                                                         badpix=galsim.Image(badpix.transpose(), scale=psf_scale),
                                                         sky_var=float(sky_var),
                                                         shear_est=method,
                                                         guess_sig_gal=gal_sig * gal_scale / psf_scale)

        if method == "KSB":

            shear_estimate = get_KSB_shear_estimate(galsim_shear_estimate, psf_scale)

        elif method == "REGAUSS":

            shear_estimate = get_REGAUSS_shear_estimate(galsim_shear_estimate, psf_scale)

        else:
            raise RuntimeError("Invalid shear estimation method for GalSim: " + str(method))

        # Add the data quality flags to the object
        shear_estimate.flags |= data_quality_flags

        # Set the S/N of the estimate
        shear_estimate.snr = signal_to_noise

        # Adjust the error to apply custom downweighting of low-S/N galaxies
        shear_estimate.gerr = np.sqrt(shear_estimate.gerr ** 2 + get_downweight_error(signal_to_noise))

        # Correct the shear estimate for x and y from resampled stamp
        shear_estimate.x = (gal_stamp.shape[0] / 2 -
                            ((resampled_gal_stamp.shape[0] / 2 - shear_estimate.x) * psf_scale / gal_scale))
        shear_estimate.y = (gal_stamp.shape[1] / 2 -
                            ((resampled_gal_stamp.shape[1] / 2 - shear_estimate.y) * psf_scale / gal_scale))

        # Correct the estimate for WCS shear and rotation
        try:
            correct_for_wcs_shear_and_rotation(shear_estimate, gal_stamp,
                                               sim_sc4_fix=sim_sc4_fix)
        except RuntimeError as e:
            # Report the error
            logger.warning(str(e))
            # Flag that we couldn't correct for the distortion
            shear_estimate.flags |= flags.flag_cannot_correct_distortion

    except RuntimeError as e:

        # If we face any error, set up a proper shear_estimate object to return
        shear_estimate = deepcopy(error_shear_estimate)

        # Set the flag of the shear_estimate object based on the type of error
        if str(e) == "HSM Error: Error: too many iterations in adaptive moments\n":
            # This is normally caused by the PSF being too small
            shear_estimate.flags |= flags.flag_psf_too_small | data_quality_flags
        else:
            raise

        logger.debug(str(e))

    logger.debug("Exiting get_shear_estimate")

    return shear_estimate


def inv_var_stack(a, a_err):

    logger = getLogger(__name__)
    logger.debug("Entering inv_var_stack")

    if a is None:
        return None, None

    a_inv_var = 1 / a_err ** 2

    inv_a_inv_var_sum = 1. / a_inv_var.sum()

    a_m = (a * a_inv_var).sum() * inv_a_inv_var_sum

    a_m_err = sqrt(inv_a_inv_var_sum)

    logger.debug("Exiting inv_var_stack")

    return a_m, a_m_err


def GS_estimate_shear(data_stack, training_data, method, workdir, sim_sc4_fix=False, debug=False):

    logger = getLogger(__name__)
    logger.debug("Entering GS_estimate_shear")

    shear_estimates_table = initialisation_methods[method]()
    tf = table_formats[method]

    if scale_label in data_stack.stacked_image.header:
        stacked_gal_scale = data_stack.stacked_image.header[scale_label]
    else:
        log_no_galaxy_scale()
        stacked_gal_scale = default_galaxy_scale

    if scale_label in data_stack.exposures[0].psf_data_hdulist[2].header:
        psf_scale = data_stack.exposures[0].psf_data_hdulist[2].header[scale_label]
    else:
        log_no_psf_scale()
        psf_scale = default_psf_scale

    # can i index data_stack with galaxy ids ???
    # data_stack.detections_catalogue=data_stack.detections_catalogue[galids]
    row_index = 0

    # Loop over galaxies and get an estimate for each one
    for row in data_stack.detections_catalogue:

        if debug and row_index > 100:
            logger.debug("Debug mode enabled, so exiting GS_estimate_shear early")
            break
        else:
            row_index += 1
            if (row_index - 1) % 100 == 0:
                logger.info("Calculating shear for galaxy " + str(row_index - 1))
            else:
                logger.debug("Calculating shear for galaxy " + str(row_index - 1))

        gal_id = row[mfc_tf.ID]
        gal_x_world = row[mfc_tf.gal_x_world]
        gal_y_world = row[mfc_tf.gal_y_world]

        # Get a stack of the galaxy images
        gal_stamp_stack = data_stack.extract_stamp_stack(x_world=gal_x_world,
                                                         y_world=gal_y_world,
                                                         width=stamp_size,
                                                         x_buffer=x_buffer,
                                                         y_buffer=y_buffer,)

        # Get stacks of the psf images
        bulge_psf_stack, disk_psf_stack = data_stack.extract_psf_stacks(gal_id=gal_id,
                                                                        make_stacked_psf=True,)

        # First, do we have any data at all?
        if gal_stamp_stack.is_empty():
            stack_shear_estimate = deepcopy(error_shear_estimate)
            stack_shear_estimate.flags |= flags.flag_no_data
            stack_x_world, stack_y_world = gal_x_world, gal_y_world
        else:

            # Get the shear estimate from the stacked image

            stacked_gal_stamp = gal_stamp_stack.stacked_image
            stacked_gal_stamp.add_default_header()

            # Make sure the wcs is correct
            stacked_gal_stamp.wcs = data_stack.stacked_image.wcs
            stacked_gal_stamp.galsim_wcs = data_stack.stacked_image.galsim_wcs

            stacked_bulge_psf_stamp = bulge_psf_stack.stacked_image
            stacked_bulge_psf_stamp.add_default_header()

            stacked_disk_psf_stamp = disk_psf_stack.stacked_image
            stacked_disk_psf_stamp.add_default_header()

            # Note the galaxy scale and gain in the stamp's header
            if scale_label in data_stack.stacked_image.header:
                stacked_gal_stamp.header[scale_label] = data_stack.stacked_image.header[scale_label]
            else:
                stacked_gal_stamp.header[scale_label] = default_galaxy_scale
            stacked_gal_stamp.header[gain_label] = mdb.get_mdb_value(mdb.mdb_keys.vis_gain)

            try:
                stack_shear_estimate = get_shear_estimate(stacked_gal_stamp,
                                                          stacked_bulge_psf_stamp,  # FIXME Handle colour gradients
                                                          gal_scale=stacked_gal_scale,
                                                          psf_scale=psf_scale,
                                                          ID=gal_id,
                                                          method=method,
                                                          stacked=True,
                                                          sim_sc4_fix=sim_sc4_fix)
            except RuntimeError as e:
                # For any unidentified errors, flag as such and return appropriate values
                logger.warning("Per-galaxy failsafe block triggered with exception: " + str(e))
                stack_shear_estimate = deepcopy(error_shear_estimate)
                stack_shear_estimate.flags |= flags.flag_unclassified_failure

            # Get world coordinates

            # Use the updated position if we had a successful fit
            if not bool(stack_shear_estimate.flags & flags.failure_flags):
                stack_x_world, stack_y_world = stacked_gal_stamp.pix2world(
                    stack_shear_estimate.x, stack_shear_estimate.y)
            else:
                # Otherwise use the original position
                stack_x_world, stack_y_world = gal_x_world, gal_y_world

            if get_exposure_estimates:

                num_exposures = len(data_stack.exposures)

                # Get estimates for each exposure

                g1s = np.zeros(num_exposures)
                g2s = np.zeros(num_exposures)
                gerrs = np.zeros(num_exposures)
                g1_errs = np.zeros(num_exposures)
                g2_errs = np.zeros(num_exposures)
                g1g2_covars = np.zeros(num_exposures)
                res = np.zeros(num_exposures)
                snrs = np.zeros(num_exposures)
                x_worlds = np.zeros(num_exposures)
                y_worlds = np.zeros(num_exposures)

                for x in range(num_exposures):

                    gal_stamp = gal_stamp_stack.exposures[x]
                    if gal_stamp is None:
                        continue
                    gal_stamp.header[scale_label] = data_stack.stacked_image.header[scale_label]
                    gal_stamp.header[gain_label] = mdb.get_mdb_value(mdb.mdb_keys.vis_gain)

                    # Make sure the wcs is correct
                    gal_stamp.wcs = gal_stamp.parent_image.wcs
                    gal_stamp.galsim_wcs = gal_stamp.parent_image.galsim_wcs

                    bulge_psf_stamp = bulge_psf_stack.exposures[x]
                    disk_psf_stamp = disk_psf_stack.exposures[x]

                    try:
                        shear_estimate = get_shear_estimate(gal_stamp,
                                                            bulge_psf_stamp,  # FIXME Handle colour gradients
                                                            gal_scale=stacked_gal_scale,
                                                            psf_scale=psf_scale,
                                                            ID=gal_id,
                                                            method=method)
                    except RuntimeError as e:
                        # For any unidentified errors, flag as such and return appropriate values
                        shear_estimate = deepcopy(error_shear_estimate)
                        shear_estimate.flags |= flags.flag_unclassified_failure

                    g1s[x] = shear_estimate.g1
                    g2s[x] = shear_estimate.g2
                    gerrs[x] = shear_estimate.gerr
                    g1_errs[x] = shear_estimate.g1_err
                    g2_errs[x] = shear_estimate.g2_err
                    g1g2_covars[x] = shear_estimate.g1g2_covar
                    res[x] = shear_estimate.re
                    snrs[x] = shear_estimate.snr

                    # Use the updated position if we had a successful fit
                    if not bool(stack_shear_estimate.flags & flags.failure_flags):
                        exp_x_world, exp_y_world = gal_stamp.pix2world(shear_estimate.x, shear_estimate.y)
                    else:
                        # Otherwise use the original position
                        exp_x_world, exp_y_world = gal_x_world, gal_y_world

                    x_worlds.append(exp_x_world)
                    y_worlds.append(exp_y_world)

                g1, gerr = inv_var_stack(g1s, gerrs)
                g2, _ = inv_var_stack(g2s, gerrs)
                re, _ = inv_var_stack(res, gerrs)
                snr, _ = inv_var_stack(snrs, gerrs)
                x_world = inv_var_stack(x_worlds, gerrs)
                y_world = inv_var_stack(y_worlds, gerrs)

        # Add this row to the estimates table (for now just using stack values)
        shear_estimates_table.add_row({tf.ID: gal_id,
                                       tf.g1: stack_shear_estimate.g1,
                                       tf.g2: stack_shear_estimate.g2,
                                       tf.g1_err: np.sqrt(stack_shear_estimate.g1_err ** 2 + training_data.e1_var),
                                       tf.g2_err: np.sqrt(stack_shear_estimate.g2_err ** 2 + training_data.e2_var),
                                       tf.g1g2_covar: stack_shear_estimate.g1g2_covar,
                                       tf.fit_class: 2,  # Unknown type, since we can't distinguish stars and galaxies
                                       tf.fit_flags: stack_shear_estimate.flags,
                                       tf.re: stack_shear_estimate.re,
                                       tf.snr: stack_shear_estimate.snr,
                                       tf.ra: stack_x_world,
                                       tf.dec: stack_y_world,
                                       })

    logger.info("Finished estimating shear.")

    logger.debug("Exiting GS_estimate_shear")

    return shear_estimates_table
