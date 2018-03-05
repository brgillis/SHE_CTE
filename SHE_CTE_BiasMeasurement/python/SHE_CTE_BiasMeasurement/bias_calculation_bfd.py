""" @file bias_calculation.py

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

from SHE_PPT.math import linregress_with_errors
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
 hdu_targ=fits.open("targettest_m_0.3_c_0.01_gdist.fits")
hdu_temp=fits.open("templatetest.fits")
hdu_targ_PQR = fits.open("targettestPQR_m_0.3_c_0.01_gdist.fits")
g1t=hdu_targ[2].data['g1_true']
g2t=hdu_targ[2].data['g2_true']

P=hdu_targ_PQR[1].data['PQR'][:,0]
Q1=hdu_targ_PQR[1].data['PQR'][:,1]
Q2=hdu_targ_PQR[1].data['PQR'][:,2]
R11=hdu_targ_PQR[1].data['PQR'][:,3]
R12=hdu_targ_PQR[1].data['PQR'][:,4]
R22=hdu_targ_PQR[1].data['PQR'][:,5]

term1a = Q1/P
term1b = Q2/P
term2a = (Q1 * Q1)/P**2 - R11/P
term2b = (Q1 * Q2)/P**2 - R12/P
term2c = (Q2 * Q1)/P**2 - R12/P
term2d = (Q2 * Q2)/P**2 - R22/P

sum1= np.sum(g1t*term1a+g2t*term1b)
sum2= np.sum(term1a + term1b)
sum3= np.sum(g1t*(g1t*term2a+g2t*term2c)+g2t*(g1t*term2b+g2t*term2d))
sum4= np.sum(g1t*term2a+g2t*term2c + g1t*term2b+g2t*term2d)
sum5= np.sum(g1t*(term2a + term2c) + g2t*(term2b+term2d))
sum6= np.sum(term2a+term2b+term2c+term2d)


bb=np.matrix([[sum1],[sum2]])
AA=np.matrix([[sum3,sum4],[sum5,sum6]])

AaAAinv=np.linalg.inv(AA)


print("m = %s +- %s" %((AAinv*bb)[0,0]-1.0,np.sqrt(AAinv[0,0])))
print("c = %s +- %s" %((AAinv*bb)[1,0],np.sqrt(AAinv[1,1])))
   
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
        
