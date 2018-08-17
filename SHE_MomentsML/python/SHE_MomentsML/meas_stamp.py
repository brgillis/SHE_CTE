#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

"""
File: python/SHE_MomentsML/meas_stamp.py

Created on: 09/08/17
Author: Malte Tewes
"""

import logging

import galsim

from SHE_MomentsML import utils_table
import numpy as np


logger = logging.getLogger(__name__)


def prep_galsim_adamom(cat, prefix="adamom_"):
    """Adds columns required for galsim_adamom

    """
    utils_table.add_cols(cat, ["flag"], dtype=int, prefix=prefix, masked=False)
    utils_table.add_cols(cat, ["flux", "x", "y", "g1", "g2", "sigma", "rho4", "logflux", "g", "theta"],
                         dtype=float, prefix=prefix, masked=True)


def galsim_adamom(stamp, src, prefix="adamom_"):
    """Measures galsim adamom features on a stamp, writing results in the table row src

    Parameters
    ----------
    stamp : SHEImage
        The image on which to measure moments
    src : Table row
        Row of an astropy Table in which to write results
    prefix : str


    FLAGS
    -----
        0 : OK
        1 : 
        2 : Centroid far from center
        3 : Negative flux
        4 : Did not converge

    """
    bkg_subtracted_stamp_data = stamp.data - stamp.background_map
    
    gs_stamp = galsim.image.Image(bkg_subtracted_stamp_data, make_const=True)
    try:
        try:  # First we try defaults:
            res = galsim.hsm.FindAdaptiveMom(gs_stamp)
        except:  # We change a bit the settings:
            logger.debug("HSM defaults failed, retrying with larger sigma...")
            hsmparams = galsim.hsm.HSMParams(max_mom2_iter=1000)
            res = galsim.hsm.FindAdaptiveMom(gs_stamp, guess_sig=15.0, hsmparams=hsmparams)
    except:  # If this also fails, we give up:
        logger.debug("HSM second attempt also failed")
        src[prefix + "flag"] = 4
        return

    # Note how we invert x and y and change signs of g1 !
    src[prefix + "flux"] = res.moments_amp
    src[prefix + "x"] = res.moments_centroid.y - 0.5  # Compensating for GalSim's default origin
    src[prefix + "y"] = res.moments_centroid.x - 0.5  # Center of first pixel is at (0.5, 0.5), not (1, 1)
    src[prefix + "g1"] = -1.0 * res.observed_shape.g1
    src[prefix + "g2"] = res.observed_shape.g2
    src[prefix + "sigma"] = res.moments_sigma
    src[prefix + "rho4"] = res.moments_rho4

    # We shift x and y by the offset
    src[prefix + "x"] += stamp.offset[0]
    src[prefix + "y"] += stamp.offset[1]

    # We add some derived quantities, as they turn out to be frequently useful
    src[prefix + "logflux"] = np.ma.log10(src[prefix + "flux"])  # Negative values get masked.
    src[prefix + "g"] = np.ma.hypot(src[prefix + "g1"], src[prefix + "g2"])
    src[prefix + "theta"] = 0.5 * np.arctan2(src[prefix + "g2"], src[prefix + "g1"])

    # If we made it to this point, we check that the centroid is roughly ok:
    stamp_center_x = stamp.shape[0] / 2.0
    stamp_center_y = stamp.shape[1] / 2.0
    if np.hypot(stamp_center_x - src[prefix + "x"], stamp_center_y - src[prefix + "y"]) > 10.0:
        src[prefix + "flag"] = 2

    if src[prefix + "flux"] < 0:
        src[prefix + "flag"] = 3


def mad(nparray):
    """The Median Absolute Deviation

    Multiply this by 1.4826 to convert into an estimate of the Gaussian std.
    """
    return np.median(np.fabs(nparray - np.median(nparray)))


def prep_skystats(cat, prefix="skystats_"):
    """Adds columns required for skystats

    """
    utils_table.add_cols(cat, ["flag"], dtype=int, prefix=prefix, masked=False)
    utils_table.add_cols(cat, ["std", "smad", "mean", "med", "stampsum"], dtype=float, prefix=prefix, masked=True)


def skystats(stamp, src, prefix="skystats_"):
    """Measures some statistics of pixels along the edge of an array

    Note that "smad" (s for scaled) is already rescaled by 1.4826 to be comparable with std.
    """
    bkg_subtracted_stamp_data = stamp.data - stamp.background_map
    a = bkg_subtracted_stamp_data  # a numpy array
    edgepixels = np.concatenate([
        a[0, 1:],  # left
        a[-1, 1:],  # right
        a[:, 0],  # bottom
        a[1:-1, -1]  # top
    ])
    assert len(edgepixels) == 2 * (a.shape[0] - 1) + 2 * (a.shape[0] - 1)

    src[prefix + "stampsum"] = np.sum(a)
    src[prefix + "std"] = np.std(edgepixels)
    src[prefix + "smad"] = 1.4826 * mad(edgepixels)
    src[prefix + "mean"] = np.mean(edgepixels)
    src[prefix + "med"] = np.median(edgepixels)


def prep_snr(cat, prefix="adamom_"):
    """Adds columns required for snr

    """
    utils_table.add_cols(cat, ["snr"], dtype=float, prefix=prefix, masked=True)


def snr(stamp, src, adamom_prefix="adamom_", skystats_prefix="skystats_", prefix="adamom_", gain=None):
    """Computes a S/N from previously measured features by galsim_adamom and skystats.

    We use the standard MomentsML S/N definition, where the noise contributes in a circular aperture of radius
    three times the measured half-light radius.
    """
    if gain is None:
        gain = stamp.header["gain"]

    sourcefluxes = src[adamom_prefix + "flux"] * gain
    skynoisefluxes = (src[skystats_prefix + "smad"] * gain) ** 2  # per pixel

    areas = np.pi * (3.0 * 1.1774 * src[adamom_prefix + "sigma"]) ** 2

    noises = np.sqrt(sourcefluxes + areas * skynoisefluxes)

    src[prefix + "snr"] = sourcefluxes / noises
