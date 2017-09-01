""" @file test_estimate_shear.py

    Created 1 Sep 2017

    Unit tests for the control shear estimation methods.

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

import pytest
import numpy as np

import galsim

from SHE_PPT.she_image import SHEImage
from SHE_PPT.magic_values import scale_label

from SHE_CTE_ShearEstimation.estimate_shear import (get_resampled_image,inv_var_stack)

class TestCase:
    """


    """


    def test_get_resampled_image(self):
        
        for ss_scale in (0.2,0.4):
            for rb_scale in (0.8,1.2):
        
                ss_factor = int(round(rb_scale/ss_scale))
                
                # Make a mock subsampled image
                ss_data = np.zeros((3*ss_factor,3*ss_factor))
                for i in range(3):
                    for j in range(3):
                        ss_data[ ss_factor*i:ss_factor*i+ss_factor, ss_factor*j:ss_factor*j+ss_factor ] += i + 3*j + 1
                        
                ss_data /= ss_data.sum()
                        
                ss_image = SHEImage(ss_data)
                ss_image.header[scale_label] = 1./ss_factor
                
                # Try rebinning it
                rb_image = get_resampled_image(ss_image,1.)
                rb_data = rb_image.data
                rb_data /= rb_data.sum()
                
                # Check the sections are close
                for i in range(3):
                    for j in range(3):
                        assert np.isclose( ss_data[ ss_factor*i:ss_factor*i+ss_factor, ss_factor*j:ss_factor*j+ss_factor ].sum(),
                                           rb_data[ i:i+1, j:j+1 ].sum(),
                                           rtol=0.2 )
        
    def test_inv_var_stack(self):
        
        a = np.array([1.0, 2.0, 2.0])
        a_err = np.array([1.0, 0.5, 2.0])
        
        # Test expected result
        
        a_m, a_m_err = inv_var_stack(a,a_err)
        
        assert np.isclose( a_m, 1.8095238095238095)
        assert np.isclose( a_m_err, 0.19047619047619047)
        
        # Test addition
        
        ap1_m, ap1_m_err = inv_var_stack(a+1,a_err)
        
        assert np.isclose( ap1_m, a_m+1 )
        assert np.isclose( ap1_m_err, a_m_err )
        
        # Test multiplication
        
        at2_m, at2_m_err = inv_var_stack(a*2,a_err)
        
        assert np.isclose( at2_m, a_m*2 )
        assert np.isclose( at2_m_err, a_m_err )
        
        # Test changing the error
        
        aet2_m, aet2_m_err = inv_var_stack(a,a_err*2)
        
        assert np.isclose( aet2_m, a_m )
        assert np.isclose( aet2_m_err, a_m_err*2 )
        
    def test_get_shear_estimate(self):
        
        sky_var = 0
        
        # Set up the galaxy profile we'll be using
        base_gal = galsim.Sersic(n=1, half_light_radius=2.0)
        
        # Set up the psf we'll be using and a subsampled image of it
        psf = galsim.Airy(lam_over_diam=0.085)
        
        ss_psf_image = galsim.Image(512,512,scale=0.02)
        psf.drawImage(ss_psf_image,use_true_center=False)
        
        psf_stamp = SHEImage(ss_psf_image.array.transpose())
        psf_stamp.header[scale_label] = ss_psf_image.scale
        
        for method in "KSB", "REGAUSS":
            for g1, g2 in ((0.,0.),
                           (0.1,0.),
                           (0.,-0.1)):
                
                # Draw the galaxy
                observed_gal = galsim.Convolve([base_gal,psf])
                observed_gal_image = galsim.Image(100,100,scale=0.10)
                observed_gal.drawImage(observed_gal_image)
                
                gal_stamp = SHEImage(observed_gal_image.array.transpose())
                gal_stamp.header[scale_label] = observed_gal_image.scale
                
                # Get the shear estimate
                est_g1, est_g2, _gerr = get_shear_estimate(gal_stamp, psf_stamp, sky_var)
                
                assert np.isclose( est_g1, g1, rtol=0.2, atol=0.01 )
