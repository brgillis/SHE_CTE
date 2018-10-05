""" @file measure_statistics_test.py

    Created 13 September 2018

    Unit tests for measuring shear bias statistics.
"""

__updated__ = "2018-09-13"

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
import os
import pytest
from pytest import approx

from SHE_PPT import products
from SHE_PPT import table_formats

from SHE_PPT.table_formats.details import tf as datf
from SHE_PPT.table_formats.shear_estimates import tf as setf
from SHE_PPT.table_formats.bfd_moments import tf as bfdf

from SHE_PPT.file_io import write_xml_product, read_xml_product
from SHE_PPT.math import BiasMeasurements, LinregressResults, linregress_with_errors

from SHE_CTE_BiasMeasurement.measure_statistics import measure_statistics_from_args
from SHE_CTE_BiasMeasurement.bfd_statistics_calculation import calculate_bfd_shear_bias_statistics

class TestBFDMeasureStatistics:
    """Tests for Running BFD Statistics Measurements.
    Actual tests of method are in SHE_BFD
    """

    @classmethod
    def setup_class(cls):

        # Set up some mock data for the test
        cls.details = table_formats.details.initialise_details_table()
        cls.bfd_estimates = table_formats.bfd_moments.initialise_bfd_moments_table(optional_columns=[bfdf.bfd_pqr])
        # outputs that expect edit
        cls.ex_b1 = -0.0033333
        cls.ex_b2 = -0.1666666
        cls.ex_b3 = -0.00625
        cls.ex_b4 =  0.2083333

        cls.ex_A11 = -0.0003555
        cls.ex_A12 = -0.0177777
        cls.ex_A13 = 9.58333e-5
        cls.ex_A14 = -0.0031944
        cls.ex_A22 = -0.8888888
        cls.ex_A23 = 0.00479166
        cls.ex_A24 = -0.1597222
        cls.ex_A33 = -0.0003734375
        cls.ex_A34 = 0.0124479166
        cls.ex_A44 = -0.414930555
        # input shears
        g1_true = np.array([0.02])
        g2_true = np.array([-0.03])
        pqr=np.zeros((1,6))
        pqr[0,:]=np.array([2.4,-0.4,0.5,2.2,0.3,1.1])

        for i in range(len(g1_true)):
            cls.details.add_row(vals={datf.ID: i,
                                      datf.group_ID: i,
                                      datf.g1:0.,
                                      datf.g2:0.})
            cls.bfd_estimates.add_row(vals={bfdf.ID: i,
                                            bfdf.bfd_pqr:np.zeros((1,6))})


        # Save details data
        cls.details[datf.g1] = g1_true
        cls.details[datf.g2] = g2_true

        # Save shear_estimates data
        cls.bfd_estimates[bfdf.bfd_pqr] = pqr

        return

    def test_calculate_shear_bias_statistics(self):
        """Try using the calculate_shear_bias_statistics function and check the results.
        """
        bfd_stats,extra = calculate_bfd_shear_bias_statistics(self.bfd_estimates, self.details)
        
        assert extra == None
        assert bfd_stats.b1 == approx(self.ex_b1,abs=1e-5)
        assert bfd_stats.b2 == approx(self.ex_b2,abs=1e-5)
        assert bfd_stats.b3 == approx(self.ex_b3,abs=1e-4)
        assert bfd_stats.b4 == approx(self.ex_b4,abs=1e-5)

        assert bfd_stats.A11 == approx(self.ex_A11,abs=1e-5)
        assert bfd_stats.A12 == approx(self.ex_A12,abs=1e-5)
        assert bfd_stats.A13 == approx(self.ex_A13,abs=1e-5)
        assert bfd_stats.A14 == approx(self.ex_A14,abs=1e-5)
        assert bfd_stats.A22 == approx(self.ex_A22,abs=1e-5)
        assert bfd_stats.A23 == approx(self.ex_A23,abs=1e-5)
        assert bfd_stats.A24 == approx(self.ex_A24,abs=1e-5)
        assert bfd_stats.A33 == approx(self.ex_A33,abs=1e-5)
        assert bfd_stats.A34 == approx(self.ex_A34,abs=1e-5)
        assert bfd_stats.A44 == approx(self.ex_A44,abs=1e-5)
        return
