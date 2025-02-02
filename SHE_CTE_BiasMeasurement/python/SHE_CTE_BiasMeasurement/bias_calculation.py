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


class CombinedBiasMeasurement(object):

    def __init__(self, g1_bias_measurement, g2_bias_measurement):

        if g1_bias_measurement is not None:
            self.m1 = g1_bias_measurement.m
            self.m1_err = g1_bias_measurement.m_err
            self.c1 = g1_bias_measurement.c
            self.c1_err = g1_bias_measurement.c_err
            self.m1c1_covar = g1_bias_measurement.mc_covar
        else:
            self.m1 = None
            self.m1_err = None
            self.c1 = None
            self.c1_err = None
            self.m1c1_covar = None

        if g2_bias_measurement is not None:
            self.m2 = g2_bias_measurement.m
            self.m2_err = g2_bias_measurement.m_err
            self.c2 = g2_bias_measurement.c
            self.c2_err = g2_bias_measurement.c_err
            self.m2c2_covar = g2_bias_measurement.mc_covar
        else:
            self.m2 = None
            self.m2_err = None
            self.c2 = None
            self.c2_err = None
            self.m2c2_covar = None

    def get_m(self):
        return np.sqrt(self.m1**2 + self.m2**2)

    def get_c(self):
        return np.sqrt(self.c1**2 + self.c2**2)

    def get_m_err(self):
        m1s_err = 2 * self.m1_err * np.abs(self.m1)
        m2s_err = 2 * self.m2_err * np.abs(self.m2)

        ms_err = np.sqrt(m1s_err**2 + m2s_err**2)
        m = self.get_m()

        ms_ferr = ms_err / m**2
        m_ferr = ms_ferr / 2

        m_err = m_ferr * m

        return m_err

    def get_c_err(self):
        c1s_err = 2 * self.c1_err * np.abs(self.c1)
        c2s_err = 2 * self.c2_err * np.abs(self.c2)

        cs_err = np.sqrt(c1s_err**2 + c2s_err**2)
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

        sign = np.sign(self.m1c1_covar * self.m2c2_covar)

        mc_covar = sign * np.sqrt(self.get_m_err()**2 / (self.m1_err * self.m2_err) *
                                  self.get_c_err()**2 / (self.c1_err * self.c2_err) *
                                  np.abs(self.m1c1_covar * self.m2c2_covar))

        return mc_covar
