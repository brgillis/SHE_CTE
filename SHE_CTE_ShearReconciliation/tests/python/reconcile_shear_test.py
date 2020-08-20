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

from SHE_CTE_ShearReconciliation.reconcile_shear import reconcile_tables
from SHE_CTE_ShearReconciliation.reconciliation_functions import (reconcile_best,
                                                                  reconcile_shape_weight,
                                                                  reconcile_invvar)
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


def assert_rows_equal(t1, t2, i, sem_tf):
    """ Check that two rows match in shear parameters
    """

    r1 = t1.loc[i]
    r2 = t2.loc[i]
    assert r1[sem_tf.g1] == r2[sem_tf.g1], "Row " + str(i) + " doesn't match expected value for g1."
    assert r1[sem_tf.g2] == r2[sem_tf.g2], "Row " + str(i) + " doesn't match expected value for g2."
    assert r1[sem_tf.g1_err] == r2[sem_tf.g1_err], "Row " + str(i) + " doesn't match expected value for g1_err."
    assert r1[sem_tf.g2_err] == r2[sem_tf.g2_err], "Row " + str(i) + " doesn't match expected value for g2_err."
    assert r1[sem_tf.weight] == r2[sem_tf.weight], "Row " + str(i) + " doesn't match expected value for weight."

    return


class TestCase:
    """
    """

    @classmethod
    def setup_class(cls):

        # Set up a mock catalogue with 20 objects, but say only the first 19 are in the tile
        cls.object_ids_in_tile = frozenset(np.arange(19, dtype=int))

        # Set up "true" values for g1 and g2, which each other table will be biased from
        true_g1 = np.linspace(-0.8, 0.8, num=20, dtype=">f4")
        true_g2 = np.where(true_g1 < 0, 0.8 + true_g1, -0.8 + true_g1)

        # We'll set up some mock tables from each method, using the same values for each method
        cls.sem_table_lists = {}
        for sem in sem_names:

            tf = sem_tfs[sem]

            tables = []

            # In brief - Table 1 and 2 cancel out for invvar weighting, but not shape weighting. Table 3 dominates,
            # Table 4 has NaN errors and 0 estimates, Table 5 has NaN estimates and errors and Inf weight,
            # Table 6 has impossible (negative) values for weight and errors
            for i_min, i_max, g1_offset, g2_offset, g1_err, g2_err, weight in ((0, 9, 0.01, -0.01, 0.01, 0.01, 15.974440894568689),
                                                                               (6, 12, -0.04, 0.04, 0.02, 0.02, 15.89825119236884),
                                                                               (8, 19, 0.1, 0.2, 0.001, 0.001, 99.9900009999),
                                                                               (0, 19, 0, 0, np.nan, np.nan, np.nan),
                                                                               (0, 19, np.nan, np.nan, np.nan, np.nan, np.inf),
                                                                               (0, 19, 0, 0, -1, -np.inf, -1),):

                l = i_max - i_min + 1
                t = sem_initialisers[sem]()
                for _ in range(l):
                    t.add_row()

                t[tf.ID] = np.arange(l) + i_min
                t[tf.g1] = true_g1[i_min:i_max + 1] + g1_offset
                t[tf.g2] = true_g2[i_min:i_max + 1] + g2_offset
                t[tf.g1_err] = np.ones(l, dtype=">f4") * g1_err
                t[tf.g2_err] = np.ones(l, dtype=">f4") * g2_err
                t[tf.weight] = np.ones(l, dtype=">f4") * weight
                t.add_index(tf.ID)
                tables.append(t)

            cls.sem_table_lists[sem] = tables

        return

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):

        self.workdir = tmpdir.strpath
        return

    def test_reconcile_best(self):

        for sem in sem_names:

            sem_table_list = self.sem_table_lists[sem]
            reconciled_catalog = reconcile_tables(shear_estimates_tables=sem_table_list,
                                                  shear_estimation_method=sem,
                                                  object_ids_in_tile=self.object_ids_in_tile,
                                                  reconciliation_function=reconcile_best,
                                                  workdir=self.workdir)

            assert(len(reconciled_catalog) == 19)

            sem_tf = sem_tfs[sem]

            reconciled_catalog.add_index(sem_tf.ID)

            # Row 1 should exactly match results from table 0
            assert_rows_equal(reconciled_catalog, sem_table_list[0], 1, sem_tf)

            # Also for 6, which overlaps with table 1, but table 0 has higher weight
            assert_rows_equal(reconciled_catalog, sem_table_list[0], 6, sem_tf)

            # Row 8 should be from table 2, which has the highest weight
            assert_rows_equal(reconciled_catalog, sem_table_list[2], 8, sem_tf)

        return

    def test_reconcile_shape_weight(self):

        for sem in sem_names:

            sem_table_list = self.sem_table_lists[sem]
            reconciled_catalog = reconcile_tables(shear_estimates_tables=sem_table_list,
                                                  shear_estimation_method=sem,
                                                  object_ids_in_tile=self.object_ids_in_tile,
                                                  reconciliation_function=reconcile_shape_weight,
                                                  workdir=self.workdir)

            assert(len(reconciled_catalog) == 19)

            sem_tf = sem_tfs[sem]

            reconciled_catalog.add_index(sem_tf.ID)

            # TODO - add tests of results here

            # Row 1 should exactly match results from table 0, since there's no other data for it
            assert_rows_equal(reconciled_catalog, sem_table_list[0], 1, sem_tf)

            # Row 18 should exactly match results from table 2, since there's no other data for it
            assert_rows_equal(reconciled_catalog, sem_table_list[2], 18, sem_tf)

        return

    def test_reconcile_invvar(self):

        for sem in sem_names:

            sem_table_list = self.sem_table_lists[sem]
            reconciled_catalog = reconcile_tables(shear_estimates_tables=sem_table_list,
                                                  shear_estimation_method=sem,
                                                  object_ids_in_tile=self.object_ids_in_tile,
                                                  reconciliation_function=reconcile_shape_weight,
                                                  workdir=self.workdir)

            assert(len(reconciled_catalog) == 19)

            sem_tf = sem_tfs[sem]

            reconciled_catalog.add_index(sem_tf.ID)

            # Row 1 should exactly match results from table 0, since there's no other data for it
            assert_rows_equal(reconciled_catalog, sem_table_list[0], 1, sem_tf)

            # Row 18 should exactly match results from table 2, since there's no other data for it
            assert_rows_equal(reconciled_catalog, sem_table_list[2], 18, sem_tf)

        return
