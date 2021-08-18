""" @file estimate_shear_test.py

    Created 1 Sep 2017

    Unit tests for the control shear estimation methods.
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

from copy import deepcopy
import gc

from SHE_PPT import flags
from SHE_PPT import mdb
from SHE_PPT.constants.fits import SCALE_LABEL, GAIN_LABEL
from SHE_PPT.file_io import find_file
from SHE_PPT.she_image import SHEImage
from astropy.io import fits
import galsim
import pytest

from SHE_CTE_ShearEstimation.galsim_estimate_shear import (get_resampled_image, inv_var_stack,
                                                           get_shear_estimate)
import numpy as np


class TestCase:
    """


    """

    @pytest.fixture(autouse=True)
    def setup(self):
        """ Set up a default galaxy stamp and PSF stamp for testing.
        """

        mdb.init(mdb_files=find_file("WEB/SHE_PPT_8_7/sample_mdb-SC8.xml"))

        self.sky_var = 0
        self.bkg_level = 1000
        self.psf_pixel_scale = 0.02
        self.gal_pixel_scale = 0.10

        self.xs = 100
        self.ys = 100

        self.psf_xs = 250
        self.psf_ys = 250

        self.g1 = 0
        self.g2 = 0

        self.gal_ID = 4

        # Set up the galaxy profile we'll be using
        self.base_gal = galsim.Sersic(n=1, half_light_radius=0.5)

        # Set up the psf we'll be using and a subsampled image of it
        self.psf = galsim.Airy(lam_over_diam=0.2)

        self.ss_psf_image = galsim.Image(self.psf_xs, self.psf_ys, scale=self.psf_pixel_scale)
        self.psf.drawImage(self.ss_psf_image, use_true_center=False)

        self.bkg_image = galsim.Image(self.xs, self.ys, scale=self.gal_pixel_scale) + self.bkg_level

        self.psf_stamp = SHEImage(self.ss_psf_image.array.transpose())
        self.psf_stamp.add_default_header()
        self.psf_stamp.header[SCALE_LABEL] = self.ss_psf_image.scale

        # Draw the default galaxy
        self.observed_gal = galsim.Convolve([self.base_gal.shear(g1=self.g1, g2=self.g2), self.psf])
        self.observed_gal_image = galsim.Image(self.xs, self.ys, scale=self.gal_pixel_scale)
        self.observed_gal.drawImage(self.observed_gal_image, use_true_center=False)

        self.observed_gal_image += self.bkg_image

        self.gal_stamp = SHEImage(self.observed_gal_image.array.transpose(),
                                  mask=np.zeros_like(self.observed_gal_image.array.transpose(), dtype=np.int8),
                                  segmentation_map=self.gal_ID * np.ones_like(
                                  self.observed_gal_image.array.transpose(), dtype=np.int8),
                                  background_map=self.bkg_image.array.transpose(),
                                  noisemap=0.0001 * np.ones_like(
                                      self.observed_gal_image.array.transpose(), dtype=float),
                                  header=fits.Header())
        self.gal_stamp.add_default_header()
        self.gal_stamp.header[SCALE_LABEL] = self.observed_gal_image.scale
        self.gal_stamp.header[GAIN_LABEL] = 1.0

        return

    def test_get_resampled_image(self):

        for ss_scale in (0.2, 0.4):
            for rb_scale in (0.8, 1.2):

                ss_factor = int(round(rb_scale / ss_scale))

                # Make a mock subsampled image
                ss_data = np.zeros((3 * ss_factor, 3 * ss_factor))
                for i in range(3):
                    for j in range(3):
                        ss_data[ss_factor * i:ss_factor * i + ss_factor,
                                ss_factor * j:ss_factor * j + ss_factor] += i + 3 * j + 1

                ss_data /= ss_data.sum()

                ss_image = SHEImage(ss_data)
                ss_image.add_default_header()
                ss_image.header[SCALE_LABEL] = 1. / ss_factor

                # Try rebinning it
                rb_image = get_resampled_image(
                    ss_image, 1., ss_image.shape[0] / ss_factor, ss_image.shape[1] / ss_factor)
                rb_data = rb_image.data
                rb_data /= rb_data.sum()

                # Check the sections are close
                for i in range(3):
                    for j in range(3):
                        assert np.isclose(ss_data[ss_factor * i:ss_factor * i + ss_factor, ss_factor * j:ss_factor * j + ss_factor].sum(),
                                          rb_data[i:i + 1, j:j + 1].sum(),
                                          rtol=0.2)

    def test_inv_var_stack(self):

        a = np.array([1.0, 2.0, 2.0])
        a_err = np.array([1.0, 0.5, 2.0])

        # Test expected result

        a_m, a_m_err = inv_var_stack(a, a_err)

        assert np.isclose(a_m, 1.8095238095238095)
        assert np.isclose(a_m_err, 0.4364357804719847)

        # Test addition

        ap1_m, ap1_m_err = inv_var_stack(a + 1, a_err)

        assert np.isclose(ap1_m, a_m + 1)
        assert np.isclose(ap1_m_err, a_m_err)

        # Test multiplication

        at2_m, at2_m_err = inv_var_stack(a * 2, a_err)

        assert np.isclose(at2_m, a_m * 2)
        assert np.isclose(at2_m_err, a_m_err)

        # Test changing the error

        aet2_m, aet2_m_err = inv_var_stack(a, a_err * 2)

        assert np.isclose(aet2_m, a_m)
        assert np.isclose(aet2_m_err, a_m_err * 2)

    def test_get_shear_estimate(self):

        for method in "KSB", "REGAUSS":
            for g1, g2 in ((0., 0.),
                           (0.1, 0.),
                           (0., -0.1)):

                # Draw the galaxy
                observed_gal = galsim.Convolve([self.base_gal.shear(g1=g1, g2=g2), self.psf])
                observed_gal_image = galsim.Image(self.xs, self.ys, scale=self.gal_pixel_scale)
                observed_gal.drawImage(observed_gal_image, use_true_center=False)

                observed_gal_image += self.bkg_image

                gal_stamp = SHEImage(observed_gal_image.array.transpose(),
                                     mask=np.zeros_like(observed_gal_image.array.transpose(), dtype=np.int8),
                                     segmentation_map=self.gal_ID * np.ones_like(
                                         observed_gal_image.array.transpose(), dtype=np.int8),
                                     background_map=self.bkg_image.array.transpose(),
                                     noisemap=0.0001 * np.ones_like(
                                         observed_gal_image.array.transpose(), dtype=float),
                                     header=fits.Header())
                gal_stamp.add_default_header()
                gal_stamp.header[SCALE_LABEL] = self.observed_gal_image.scale
                gal_stamp.header[GAIN_LABEL] = 1.0

                # Get the shear estimate
                shear_estimate = get_shear_estimate(gal_stamp,
                                                    self.psf_stamp,
                                                    gal_scale=self.gal_pixel_scale,
                                                    psf_scale=self.psf_pixel_scale,
                                                    method=method,
                                                    ID=1)
                est_g1, est_g2 = shear_estimate.g1, shear_estimate.g2

                if shear_estimate.flags & flags.failure_flags:
                    raise RuntimeError("Error in shear estimate: " + bin(shear_estimate.flags))

                assert np.isclose(est_g1, g1, rtol=0.2, atol=0.01)
                assert np.isclose(est_g2, g2, rtol=0.2, atol=0.01)

    def test_shear_estimate_flags(self):
        """Test proper flagging by forcing failure conditions and then checking output flags.
        """

        # Test both methods equivalently
        for method in ("KSB", "REGAUSS"):

            for (attr, corrupt_flag, missing_flag) in (("_data", flags.flag_corrupt_science_image,
                                                        flags.flag_no_science_image),
                                                       ("mask", flags.flag_corrupt_mask,
                                                        flags.flag_no_mask),
                                                       ("background_map", flags.flag_corrupt_background_map,
                                                        flags.flag_no_background_map),
                                                       ("noisemap", flags.flag_corrupt_noisemap,
                                                        flags.flag_no_noisemap),
                                                       ("segmentation_map", flags.flag_corrupt_segmentation_map,
                                                        flags.flag_no_segmentation_map),
                                                       ):

                # Test corrupt failure
                gal_stamp = deepcopy(self.gal_stamp)
                for bad_val in (np.nan, np.inf, -2):
                    image = getattr(gal_stamp, attr)
                    try:
                        image[self.xs // 2, self.ys // 2] = bad_val
                    except ValueError as e:
                        # For the integer maps, ignore failures on assigning NaN
                        if "cannot convert float NaN to integer" in str(e):
                            continue
                        else:
                            raise
                    except OverflowError as e:
                        # For the integer maps, ignore failures on assigning inf
                        if "cannot convert float infinity to integer" in str(e):
                            continue
                        else:
                            raise

                    shear_estimate = get_shear_estimate(gal_stamp,
                                                        self.psf_stamp,
                                                        gal_scale=self.gal_pixel_scale,
                                                        psf_scale=self.psf_pixel_scale,
                                                        method=method,
                                                        ID=1)

                    if not shear_estimate.flags & corrupt_flag:
                        raise RuntimeError("Failed to raise flag " + bin(corrupt_flag) + " when expected " +
                                           "for method " + method + ".")

                # Test missing failure
                setattr(gal_stamp, attr, None)
                shear_estimate = get_shear_estimate(gal_stamp,
                                                    self.psf_stamp,
                                                    gal_scale=self.gal_pixel_scale,
                                                    psf_scale=self.psf_pixel_scale,
                                                    method=method,
                                                    ID=1)

                if not shear_estimate.flags & missing_flag:
                    raise RuntimeError("Failed to raise flag " + bin(missing_flag) + " when expected " +
                                       "for method " + method + ".")

                del gal_stamp
                gc.collect()

        return
