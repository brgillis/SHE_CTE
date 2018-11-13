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

from math import sqrt

import galsim

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_PPT.logging import getLogger
from SHE_PPT.magic_values import scale_label, gain_label
from SHE_PPT.she_image import SHEImage
from SHE_PPT.shear_utility import get_g_from_e
from SHE_PPT.table_formats.detections import tf as detf
from SHE_PPT.table_formats.shear_estimates import initialise_shear_estimates_table, tf as setf
from SHE_PPT.utility import run_only_once
import numpy as np


stamp_size = 256
x_buffer = -5
y_buffer = -5
get_exposure_estimates = False

default_galaxy_scale = 0.1
default_psf_scale = 0.02


class ShearEstimate(object):

    def __init__(self, g1, g2, gerr, re, snr, x, y):
        self.g1 = g1
        self.g2 = g2
        self.gerr = gerr
        self.re = re
        self.snr = snr
        self.x = x
        self.y = y


snr_cutoff = 15
downweight_error = 0.5
downweight_power = 4


def downweight_error(snr):
    return downweight_error / (1 + (snr / snr_cutoff)**downweight_power)


@run_only_once
def log_no_galaxy_scale():
    logger = getLogger(__name__)
    logger.warn("Cannot find pixel scale in image header. Using default value of " + str(default_galaxy_scale))


@run_only_once
def log_no_psf_scale():
    logger = getLogger(__name__)
    logger.warn("Cannot find pixel scale in PSF header. Using default value of " + str(default_psf_scale))


def get_resampled_image(initial_image, resampled_scale, resampled_nx, resampled_ny):

    if scale_label in initial_image.header:
        in_scale = initial_image.header[scale_label]
    else:
        log_no_galaxy_scale()
        in_scale = default_galaxy_scale

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

    resampled_image.header[scale_label] = resampled_scale

    return resampled_image


def KSB_estimate_shear(data_stack, training_data, calibration_data, workdir, *args, **kwargs):
    # Not using training or calibration data at this stage
    return GS_estimate_shear(data_stack, training_data=training_data, method="KSB", workdir=workdir, *args, **kwargs)


def REGAUSS_estimate_shear(data_stack, training_data, calibration_data, workdir, *args, **kwargs):
    # Not using training or calibration data at this stage
    return GS_estimate_shear(data_stack, training_data=training_data, method="REGAUSS", workdir=workdir, *args,
                             **kwargs)


def get_KSB_shear_estimate(galsim_shear_estimate, scale):

    logger = getLogger(__name__)
    logger.debug("Entering get_KSB_shear_estimate")

    g1 = galsim_shear_estimate.corrected_g1
    g2 = galsim_shear_estimate.corrected_g2
    mag = g1 ** 2 + g2 ** 2

    if mag > 1:
        raise RuntimeError("HSM Error: Magnitude of g shear is too large: " + str(mag))

    shear_estimate = ShearEstimate(galsim_shear_estimate.corrected_g1,
                                   galsim_shear_estimate.corrected_g2,
                                   galsim_shear_estimate.corrected_shape_err,
                                   galsim_shear_estimate.moments_sigma * scale,
                                   -1,
                                   galsim_shear_estimate.moments_centroid.x,
                                   galsim_shear_estimate.moments_centroid.y)

    logger.debug("Exiting get_KSB_shear_estimate")

    return shear_estimate


def get_REGAUSS_shear_estimate(galsim_shear_estimate, scale):

    logger = getLogger(__name__)
    logger.debug("Entering get_REGAUSS_shear_estimate")

    e1 = galsim_shear_estimate.corrected_e1
    e2 = galsim_shear_estimate.corrected_e2
    mag = e1 ** 2 + e2 ** 2

    if mag > 1:
        raise RuntimeError("HSM Error: Magnitude of e shear is too large: " + str(mag))

    g1, g2 = get_g_from_e(e1, e2)
    gerr = galsim_shear_estimate.corrected_shape_err * np.sqrt((g1 ** 2 + g2 ** 2) / (e1 ** 2 + e2 ** 2))

    shear_estimate = ShearEstimate(g1, g2, gerr,
                                   galsim_shear_estimate.moments_sigma * scale,
                                   -1,
                                   galsim_shear_estimate.moments_centroid.x,
                                   galsim_shear_estimate.moments_centroid.y)

    logger.debug("Exiting get_REGAUSS_shear_estimate")

    return shear_estimate


def get_shear_estimate(gal_stamp, psf_stamp, gal_scale, psf_scale, ID, method):

    logger = getLogger(__name__)
    logger.debug("Entering get_shear_estimate")

    # Subtract off the background from the galaxy
    bkg_subtracted_gal_stamp_data = gal_stamp.data - gal_stamp.background_map

    # Estimate the size of the galaxy, so we can figure out how big we need to make the resampled stamp

    gal_sigs = [5.0, 2.0, 0.5, 10.0]

    for gal_sig in gal_sigs:

        try:
            gal_mom = galsim.hsm.FindAdaptiveMom(galsim.Image(bkg_subtracted_gal_stamp_data.transpose(), scale=psf_scale),
                                                 badpix=galsim.Image(
                                                     (gal_stamp.boolmask).astype(np.uint16).transpose(), scale=gal_scale),
                                                 guess_sig=gal_sig,)

            resampled_gal_stamp_size = int(5 * gal_mom.moments_sigma * gal_scale / psf_scale)
            break
        except RuntimeError as e:
            if ("HSM Error" not in str(e)):
                raise
            elif gal_sig == 10.0:
                # If it fails, it's probably because the galaxy is small, so a small size will suffice
                resampled_gal_stamp_size = 50
            else:
                continue

    # Calculate the galaxy's S/N
    a_eff = np.pi * (3 * gal_mom.moments_sigma * np.sqrt(2 * np.log(2)))
    gain = gal_stamp.header[gain_label]
    signal_to_noise = (gain * gal_mom.moments_amp / np.sqrt(gain * gal_mom.moments_amp + a_eff *
                                                            (gain * np.square(gal_stamp.noisemap.transpose()).mean())**2))

    # Get a resampled galaxy stamp
    resampled_gal_stamp = get_resampled_image(gal_stamp, psf_scale, resampled_gal_stamp_size, resampled_gal_stamp_size)

    # Get a resampled badpix map
    supersampled_badpix = SHEImage((gal_stamp.boolmask).astype(float))
    supersampled_badpix.header[scale_label] = gal_stamp.header[scale_label]
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

        # Correct the shear estimate for x and y from resampled stamp
        shear_estimate.x = (gal_stamp.shape[0] / 2 -
                            ((resampled_gal_stamp.shape[0] / 2 - shear_estimate.x) * psf_scale / gal_scale))
        shear_estimate.y = (gal_stamp.shape[1] / 2 -
                            ((resampled_gal_stamp.shape[1] / 2 - shear_estimate.y) * psf_scale / gal_scale))

        # Set the proper snr for the estimate, and use it to downweight as appropriate
        shear_estimate.snr = signal_to_noise
        shear_estimate.gerr = np.sqrt(shear_estimate.gerr**2 + downweight_error(signal_to_noise))

    except RuntimeError as e:
        if ("HSM Error" not in str(e)):
            raise

        logger.debug(str(e))
        shear_estimate = ShearEstimate(np.NaN,
                                       np.NaN,
                                       np.inf,
                                       np.NaN,
                                       np.NaN,
                                       gal_stamp.data.shape[0],
                                       gal_stamp.data.shape[1])

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


def GS_estimate_shear(data_stack, training_data, method, workdir, debug=False):

    logger = getLogger(__name__)
    logger.debug("Entering GS_estimate_shear")

    shear_estimates_table = initialise_shear_estimates_table()

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

        gal_id = row[detf.ID]
        gal_x_world = row[detf.gal_x_world]
        gal_y_world = row[detf.gal_y_world]

        # Get a stack of the galaxy images
        gal_stamp_stack = data_stack.extract_stamp_stack(x_world=gal_x_world,
                                                         y_world=gal_y_world,
                                                         width=stamp_size,
                                                         x_buffer=x_buffer,
                                                         y_buffer=y_buffer,)

        # If there's no data for this galaxy, don't add it to the catalogue at all
        if gal_stamp_stack.is_empty():
            continue

        # Get stacks of the psf images
        bulge_psf_stack, disk_psf_stack = data_stack.extract_psf_stacks(gal_id=gal_id,
                                                                        make_stacked_psf=True,)

        # Get the shear estimate from the stacked image

        stacked_gal_stamp = gal_stamp_stack.stacked_image
        stacked_bulge_psf_stamp = bulge_psf_stack.stacked_image
        stacked_disk_psf_stamp = disk_psf_stack.stacked_image

        # Note the galaxy scale in the stamp's header
        stacked_gal_stamp.header[scale_label] = data_stack.stacked_image.header[scale_label]

        shear_estimate = get_shear_estimate(stacked_gal_stamp,
                                            stacked_bulge_psf_stamp,  # FIXME Handle colour gradients
                                            gal_scale=stacked_gal_scale,
                                            psf_scale=psf_scale,
                                            ID=gal_id,
                                            method=method)

        stack_g_pix = np.matrix([[shear_estimate.g1], [shear_estimate.g2]])
        stack_re = shear_estimate.re
        stack_snr = shear_estimate.snr

        # Get world coordinates

        stack_x_world, stack_y_world = stacked_gal_stamp.pix2world(shear_estimate.x, shear_estimate.y)

        # Need to convert g1/g2 and errors to -ra/dec coordinates
        stack_rotation_matrix = stacked_gal_stamp.get_pix2world_rotation(shear_estimate.x, shear_estimate.y)
        stack_double_rotation_matrix = stack_rotation_matrix @ stack_rotation_matrix  # 2x2 so it's commutative
        stack_g_world = stack_double_rotation_matrix @ stack_g_pix

        stack_covar_pix = np.matrix([[shear_estimate.gerr, 0], [0, shear_estimate.gerr]])
        stack_covar_world = stack_double_rotation_matrix @ stack_covar_pix @ stack_double_rotation_matrix.transpose()

        if get_exposure_estimates:

            # Get estimates for each exposure

            g1s = []
            g2s = []
            gerrs = []
            res = []
            snrs = []
            x_worlds = []
            y_worlds = []

            for x in range(len(data_stack.exposures)):

                gal_stamp = gal_stamp_stack.exposures[x]
                if gal_stamp is None:
                    continue
                gal_stamp.header[scale_label] = data_stack.stacked_image.header[scale_label]
                bulge_psf_stamp = bulge_psf_stack.exposures[x]
                disk_psf_stamp = disk_psf_stack.exposures[x]

                shear_estimate = get_shear_estimate(gal_stamp,
                                                    bulge_psf_stamp,  # FIXME Handle colour gradients
                                                    gal_scale=stacked_gal_scale,
                                                    psf_scale=psf_scale,
                                                    ID=gal_id,
                                                    method=method)

                g_pix = np.matrix([[shear_estimate.g1], [shear_estimate.g2]])

                # Need to convert g1/g2 and errors to -ra/dec coordinates
                rotation_matrix = gal_stamp.get_pix2world_rotation(shear_estimate.x, shear_estimate.y)
                g_world = rotation_matrix @ (rotation_matrix @ g_pix)

                g1s.append(g_world[0])
                g2s.append(g_world[1])
                gerrs.append(shear_estimate.gerr)
                res.append(shear_estimate.re)
                snrs.append(shear_estimate.snr)

                x_world, y_world = gal_stamp.pix2world(shear_estimate.x, shear_estimate.y)
                x_worlds.append(x_world)
                y_worlds.append(y_world)

            g1s = np.array(g1s)
            g2s = np.array(g2s)
            gerrs = np.array(gerrs)
            res = np.array(res)
            snrs = np.array(snrs)
            x_worlds = np.array(x_worlds)
            y_worlds = np.array(y_worlds)

            g1, gerr = inv_var_stack(g1s, gerrs)
            g2, _ = inv_var_stack(g2s, gerrs)
            re, _ = inv_var_stack(res, gerrs)
            snr, _ = inv_var_stack(snrs, gerrs)
            x_world = inv_var_stack(x_worlds, gerrs)
            y_world = inv_var_stack(y_worlds, gerrs)

        # Add this row to the estimates table (for now just using stack values)
        shear_estimates_table.add_row({setf.ID: gal_id,
                                       setf.g1: stack_g_world[0],
                                       setf.g2: stack_g_world[1],
                                       setf.g1_err: np.sqrt(stack_covar_world[0, 0] ** 2 + training_data.e1_var),
                                       setf.g2_err: np.sqrt(stack_covar_world[1, 1] ** 2 + training_data.e2_var),
                                       setf.g1g2_covar: stack_covar_world[0, 1],
                                       setf.re: stack_re,
                                       setf.snr: stack_snr,
                                       setf.x_world: stack_x_world,
                                       setf.y_world: stack_y_world,
                                       })

    logger.info("Finished estimating shear.")

    logger.debug("Exiting GS_estimate_shear")

    return shear_estimates_table
