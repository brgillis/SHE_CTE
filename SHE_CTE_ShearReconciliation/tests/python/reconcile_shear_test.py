""" @file reconcile_shear_test.py

    Created 17 Aug 2020

    Unit tests for the shear reconciliation task.
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

from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf, initialise_mer_final_catalog
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksbm_tf, initialise_ksb_measurements_table
from SHE_PPT.table_formats.she_lensmc_measurements import tf as lmcm_tf, initialise_lensmc_measurements_table
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_formats.she_momentsml_measurements import tf as mmlm_tf, initialise_momentsml_measurements_table
from SHE_PPT.table_formats.she_regauss_measurements import tf as regm_tf, initialise_regauss_measurements_table
from SHE_PPT.table_utility import is_in_format, add_row
import pytest

import numpy as np

sem_names = ("KSB",
             "LensMC",
             "MomentsML",
             "REGAUSS")

sem_tfs = {"KSB":ksbm_tf,
           "LensMC":lmcm_tf,
           "MomentsML":mmlm_tf,
           "REGAUSS":regm_tf}

sem_initialisers = {"KSB":initialise_ksb_measurements_table,
                    "LensMC":initialise_lensmc_measurements_table,
                    "MomentsML":initialise_momentsml_measurements_table,
                    "REGAUSS":initialise_regauss_measurements_table}


class TestCase:
    """
    """

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):

        self.workdir = tmpdir.strpath

        # Set up a mock catalogue with 20 objects
        self.mer_final_catalog = initialise_mer_final_catalog(size=20)
        self.mer_final_catalog[mfc_tf.ID] = np.arange(20, dtype=int)

        # Set up "true" values for g1 and g2, which each other table will be biased from
        true_g1 = np.linspace(-0.8, 0.8, num=20)
        true_g2 = np.where(true_g1 < 0, 0.8 + true_g1, -0.8 + true_g1)

        # We'll set up some mock tables from each method, using the same values for each method
        self.sem_tables = {}
        for sem in sem_names:

            tf = sem_tfs[sem]

            tables = []

            # In brief - Table 1 and 2 cancel out for invvar weighting, but not shape weighting. Table 3 dominates,
            # Table 4 has NaN errors and 0 estimates, Table 5 has NaN estimates and errors and Inf weight,
            # Table 6 has impossible (negative) values for weight and errors
            for i_min, i_max, g1_offset, g2_offset, g1_err, g2_err, weight in ((0, 9, 0.01, -0.01, 0.01, 0.01, 1.0),
                                                                               (6, 12, -0.04, 0.04, 0.02, 0.02, 0.5),
                                                                               (8, 19, 0.1, 0.2, 0.001, 0.001, 100),
                                                                               (0, 19, 0, 0, np.nan, np.nan, np.nan),
                                                                               (0, 19, np.nan, np.nan, np.nan, np.nan, np.inf),
                                                                               (0, 19, 0, 0, -1, -np.inf, -1),):

                l = i_max - i_min + 1
                t = sem_initialisers(size=l)
                t[tf.ID] = np.arange(l) + i_min
                t[tf.g1] = true_g1[i_min:i_max + 1] + g1_offset
                t[tf.g2] = true_g2[i_min:i_max + 1] + g2_offset
                t[tf.g1] = np.ones(l) * g1_err
                t[tf.g2] = np.ones(l) * g2_err
                tables.append(t)

            self.sem_tables[sem] = tables

        return

    def test_reconcile_best(self):
        pass

    def test_reconcile_shape_weight(self):
        pass

    def test_reconcile_invvar(self):
        pass
