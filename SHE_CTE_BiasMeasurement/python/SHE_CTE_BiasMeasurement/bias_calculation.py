""" @file measurement_extraction.py

    Created 10 Apr 2017

    Function to calculate bias from a table of shear measurements
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

import numpy as np

from SHE_PPT.table_formats.shear_estimates import tf as setf

from SHE_CTE_BiasMeasurement import magic_values as mv

class BiasMeasurement(object):
    def __init__(self):
        self.m1 = None
        self.m2 = None
        self.m1_err = None
        self.m2_err = None
        self.c1 = None
        self.c2 = None
        self.c1_err = None
        self.c2_err = None
        self.m1c1_covar = None
        self.m2c2_covar = None
        
    def get_m(self):
        return np.sqrt(self.m1**2 + self.m2**2)
        
    def get_c(self):
        return np.sqrt(self.c1**2 + self.c2**2)
    
    def get_m_err(self):
        m1s_err = 2*self.m1_err*np.abs(self.m1)
        m2s_err = 2*self.m2_err*np.abs(self.m2)
        
        ms_err = np.sqrt(m1s_err**2+m2s_err**2)
        m = self.get_m()
        
        ms_ferr = ms_err / m**2
        m_ferr = ms_ferr / 2
        
        m_err = m_ferr * m
        
        return m_err
    
    def get_c_err(self):
        c1s_err = 2*self.c1_err*np.abs(self.c1)
        c2s_err = 2*self.c2_err*np.abs(self.c2)
        
        cs_err = np.sqrt(c1s_err**2+c2s_err**2)
        c = self.get_c()
        
        cs_ferr = cs_err / c**2
        c_ferr = cs_ferr / 2
        
        c_err = c_ferr * c
        
        return c_err
    
    def get_mc_covar(self):
        """
        @brief Estimate the mc covariance
        @details Wild guess using harmonic means. If we actually need this I'll look into it more.
        """
        
        sign = np.sign(self.m1c1_covar*self.m2c2_covar)
        
        mc_covar = sign*np.sqrt(self.get_m_err()**2/(self.m1_err*self.m2_err) *
                           self.get_c_err()**2/(self.c1_err*self.c2_err) *
                           np.abs(self.m1c1_covar*self.m2c2_covar) )
        
        return mc_covar
        
def compress_measurements(real_values,measurements,measurement_errors):
    """
    @brief
        Compress measurements when shape noise cancellation was used, to combine
        measurements made on the same input values
        
    @param real_values <np.ndarray>
    @param measurements <np.ndarray>
    @param measurement_errors <np.ndarray>
    
    @return compressed_real_values <np.ndarray>,
            compressed_measurements <np.ndarray>,
            compressed_measurement_errors <np.ndarray>,
    """
    
    compressed_real_values = []
    compressed_measurements = []
    compressed_measurement_errors = []
    
    N = len(real_values)
    i = 0
    
    while i < N:
        real_value = real_values[i]
        temp_measurements = []
        temp_measurement_errors = []
        
        while real_values[i]==real_value:
            # Check that the measurement is good
            new_measurement = measurements[i]
            new_measurement_error = measurement_errors[i]
            if (new_measurement > -2 and new_measurement < 2 and 
                    new_measurement_error > 0 and new_measurement_error < 1e99):
                temp_measurements.append(new_measurement)
                temp_measurement_errors.append(new_measurement_error)
            i += 1
            if i >= N:
                break
            
        temp_measurements = np.array(temp_measurements)
        temp_measurement_errors = np.array(temp_measurement_errors)
            
        temp_measurement_weights = temp_measurement_errors**-2
        
        total_weight = np.sum(temp_measurement_weights)
        
        if total_weight==0:
            continue
        
        mean_measurement = np.sum(temp_measurements*temp_measurement_weights)/total_weight
        mean_measurement_error = np.sqrt(1/np.sum(temp_measurement_weights))
        
        compressed_real_values.append(real_value)
        compressed_measurements.append(mean_measurement)
        compressed_measurement_errors.append(mean_measurement_error)
        
    return (np.array(compressed_real_values),
            np.array(compressed_measurements),
            np.array(compressed_measurement_errors))
    
def linregress_with_errors(x, y, y_err):
    """
    @brief
        Perform a linear regression with errors on the y values
    @details
        This uses a direct translation of GSL's gsl_fit_wlinear function, using
        inverse-variance weighting
        
    @param x <np.ndarray>
    @param y <np.ndarray>
    @param y_err <np.ndarray>
    
    @return slope <float>,
            intercept <float>,
            slope_err <float>,
            intercept_err <float>,
            slope_intercept_covar <float>
    """
    
    y_weights = y_err**-2
    total_weight = y_weights.sum()
    
    x_weighted_mean = np.average(x,weights=y_weights)
    y_weighted_mean = np.average(y,weights=y_weights)
    
    dx = x-x_weighted_mean
    dy = y-y_weighted_mean
    
    dx2_weighted_mean = np.average(dx**2,weights=y_weights)
    dxdy_weighted_mean = np.average(dx*dy,weights=y_weights)
    
    slope = dxdy_weighted_mean / dx2_weighted_mean
    intercept = y_weighted_mean - x_weighted_mean*slope
    
    slope_err = np.sqrt(1./(total_weight*dx2_weighted_mean))
    intercept_err = np.sqrt((1.0 + x_weighted_mean**2 / dx2_weighted_mean) / total_weight)
    slope_intercept_covar = -x_weighted_mean / (total_weight*dx2_weighted_mean)

    return slope, intercept, slope_err, intercept_err, slope_intercept_covar
        
def regress_shear_measurements(real_values, measurements, measurement_errors):
    """
    @brief
        Perform a linear regression on a set of shear measurements
        
    @param real_values <np.ndarray>
    @param measurements <np.ndarray>
    @param measurement_errors <np.ndarray>
    
    @return m <float>,
            c <float>,
            m_err <float>,
            c_err <float>,
            mc_covar <float>
    """
    
    # Compress measurements to properly account for shape noise cancellation
    (compressed_real_values,
     compressed_measurements,
     compressed_measurement_errors) = compress_measurements(real_values,measurements,measurement_errors)
     
    slope, c, m_err, c_err, mc_covar = linregress_with_errors(compressed_real_values,
                                                              compressed_measurements,
                                                              compressed_measurement_errors)
    
    m = slope - 1
    
    return m, m_err, c, c_err, mc_covar
        
def calculate_bias(all_shear_measurements):
    """
    @brief
    
        Calculate bias from a table of all shear measurements.
    @param all_shear_measurements <astropy.table.Table>
        Table of all shear measurements
        
    @return <BiasMeasurement>
    """
    
    bias_measurement = BiasMeasurement()
    
    # Get bias for both index 1 and 2 independently
    
    (bias_measurement.m1, bias_measurement.m1_err,
     bias_measurement.c1, bias_measurement.c1_err,
     bias_measurement.m1c1_covar) = (
            regress_shear_measurements(all_shear_measurements[mv.fits_table_sim_g1_label],
                                       all_shear_measurements[setf.g1],
                                       all_shear_measurements[setf.g1_err]) )
    
    (bias_measurement.m2, bias_measurement.m2_err,
     bias_measurement.c2, bias_measurement.c2_err,
     bias_measurement.m2c2_covar) = (
            regress_shear_measurements(all_shear_measurements[mv.fits_table_sim_g2_label],
                                       all_shear_measurements[setf.g2],
                                       all_shear_measurements[setf.g2_err]) )
     
    return bias_measurement 
        