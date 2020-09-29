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

import os

from SHE_PPT import products
from SHE_PPT.file_io import write_xml_product, read_xml_product, get_allowed_filename, write_listfile
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf, initialise_mer_final_catalog
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksbm_tf, initialise_ksb_measurements_table
from SHE_PPT.table_formats.she_lensmc_chains import tf as lmcc_tf, initialise_lensmc_chains_table, len_chain
from SHE_PPT.table_formats.she_lensmc_measurements import tf as lmcm_tf, initialise_lensmc_measurements_table
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_formats.she_momentsml_measurements import tf as mmlm_tf, initialise_momentsml_measurements_table
from SHE_PPT.table_formats.she_regauss_measurements import tf as regm_tf, initialise_regauss_measurements_table
from SHE_PPT.table_utility import is_in_format, add_row
from astropy.table import Table
import pytest

import SHE_CTE
from SHE_CTE_ShearReconciliation.chains_reconciliation_functions import (reconcile_chains_best,
                                                                         reconcile_chains_shape_weight,
                                                                         reconcile_chains_invvar,
                                                                         reconcile_chains_keep)
from SHE_CTE_ShearReconciliation.reconcile_shear import reconcile_shear_from_args, reconcile_tables, reconcile_chains
from SHE_CTE_ShearReconciliation.reconciliation_functions import (reconcile_best,
                                                                  reconcile_shape_weight,
                                                                  reconcile_invvar)
import numpy as np


sem_names = ("KSB",
             "LensMC",
             "MomentsML",
             "REGAUSS")

sem_tfs = {"KSB": ksbm_tf,
           "LensMC": lmcm_tf,
           "MomentsML": mmlm_tf,
           "REGAUSS": regm_tf}

sem_initialisers = {"KSB": initialise_ksb_measurements_table,
                    "LensMC": initialise_lensmc_measurements_table,
                    "MomentsML": initialise_momentsml_measurements_table,
                    "REGAUSS": initialise_regauss_measurements_table}


class Args(object):

    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])


def assert_rows_equal(t1, t2, i, sem_tf):
    """ Check that two rows match in shear parameters
    """

    r1 = t1.loc[i]
    r2 = t2.loc[i]
    assert np.isclose(r1[sem_tf.g1], r2[sem_tf.g1]), "Row " + str(i) + " doesn't match expected value for g1."
    assert np.isclose(r1[sem_tf.g2], r2[sem_tf.g2]), "Row " + str(i) + " doesn't match expected value for g2."
    assert np.isclose(r1[sem_tf.g1_err], r2[sem_tf.g1_err]), "Row " + \
        str(i) + " doesn't match expected value for g1_err."
    assert np.isclose(r1[sem_tf.g2_err], r2[sem_tf.g2_err]), "Row " + \
        str(i) + " doesn't match expected value for g2_err."
    assert np.isclose(r1[sem_tf.weight], r2[sem_tf.weight]), "Row " + \
        str(i) + " doesn't match expected value for weight."

    return


def assert_chains_rows_equal(t1, t2, i):
    """ Check that two chains rows match in shear parameters
    """

    r1 = t1.loc[i]
    r2 = t2.loc[i]
    assert np.isclose(r1[lmcc_tf.g1], r2[lmcc_tf.g1]).all(), "Row " + str(i) + " doesn't match expected value for g1."
    assert np.isclose(r1[lmcc_tf.g2], r2[lmcc_tf.g2]).all(), "Row " + str(i) + " doesn't match expected value for g2."
    assert np.isclose(r1[lmcc_tf.shape_noise], r2[lmcc_tf.shape_noise]), "Row " + \
        str(i) + " doesn't match expected value for shape noise."
    assert np.isclose(r1[lmcc_tf.weight], r2[lmcc_tf.weight]), "Row " + \
        str(i) + " doesn't match expected value for weight."

    return


class TestReconcileShear:
    """
    """

    @classmethod
    def setup_class(cls):

        # Set up a mock catalogue with 20 objects, but say only the first 19 are in the tile
        cls.ids = np.arange(19, dtype=int)
        cls.object_ids_in_tile = frozenset(cls.ids)

        # Set up "true" values for g1 and g2, which each other table will be biased from
        cls.true_g1 = np.linspace(-0.8, 0.8, num=20, dtype=">f4")
        cls.true_g2 = np.where(cls.true_g1 < 0, 0.8 + cls.true_g1, -0.8 + cls.true_g1)

        cls.shape_noise = 0.25
        cls.max_weight = 1 / cls.shape_noise ** 2

        # We'll set up some mock tables from each method, using the same values for each method
        cls.sem_table_lists = {}
        cls.chains_table_list = []
        for sem in sem_names:

            tf = sem_tfs[sem]

            make_chains = sem == "LensMC"

            tables = []

            # In brief - Table 1 and 2 cancel out for invvar weighting, but not shape weighting. Table 3 dominates,
            # Table 4 has NaN errors and 0 estimates, Table 5 has NaN estimates and errors and Inf weight,
            # Table 6 has impossible (negative) values for weight and errors
            for i_min, i_max, g1_offset, g2_offset, g1_err, g2_err, weight in ((0, 9, 0.01, -0.01, 0.01, 0.01, None),
                                                                               (6, 12, -0.04, 0.04, 0.02, 0.02, None),
                                                                               (8, 19, 0.1, 0.2, 0.001, 0.001, None),
                                                                               (0, 19, 0, 0, np.nan, np.nan, np.nan),
                                                                               (0, 19, np.nan, np.nan,
                                                                                np.nan, np.nan, np.inf),
                                                                               (0, 19, 0, 0, -1, -np.inf, -1),):

                if weight is None:
                    weight = 1 / (0.5 * (g1_err ** 2 + g2_err ** 2) + cls.shape_noise ** 2)
                l = i_max - i_min + 1

                # Create the measurements table

                # t = sem_initialisers[sem](optional_columns=(tf.shape_weight, tf.e_var, tf.shape_noise))
                t = sem_initialisers[sem](optional_columns=(tf.e_var, tf.shape_noise))
                for _ in range(l):
                    t.add_row()

                t[tf.ID] = np.arange(l) + i_min
                t[tf.g1] = cls.true_g1[i_min:i_max + 1] + g1_offset
                t[tf.g2] = cls.true_g2[i_min:i_max + 1] + g2_offset
                t[tf.g1_err] = np.ones(l, dtype=">f4") * g1_err
                t[tf.g2_err] = np.ones(l, dtype=">f4") * g2_err
                t[tf.weight] = np.ones(l, dtype=">f4") * weight
                # t[tf.shape_weight] = np.ones(l, dtype=">f4") / (g1_err**2 + g2_err**2)
                t[tf.e_var] = np.ones(l, dtype=">f4") * (g1_err**2 + g2_err**2)
                t[tf.shape_noise] = np.ones(l, dtype=">f4") * cls.shape_noise
                t.add_index(tf.ID)
                tables.append(t)

                if make_chains:
                    # Create the chains table

                    tc = initialise_lensmc_chains_table(optional_columns=(lmcc_tf.e_var, lmcc_tf.shape_noise))
                    for _ in range(l):
                        tc.add_row()

                    # Use a stable set of deviates which average to zero for consistency
                    deviates_first_half = np.random.standard_normal((l, len_chain // 2))
                    deviates_second_half = -deviates_first_half
                    deviates = np.concatenate((deviates_first_half, deviates_second_half), axis=1)

                    tc[lmcc_tf.ID] = np.arange(l) + i_min
                    tc[lmcc_tf.g1] = (np.array((cls.true_g1[i_min:i_max + 1],)).transpose() + g1_offset +
                                      deviates * g1_err).astype(lmcc_tf.dtypes[lmcc_tf.g1])
                    tc[lmcc_tf.g2] = (np.array((cls.true_g2[i_min:i_max + 1],)).transpose() + g2_offset +
                                      deviates * g2_err).astype(lmcc_tf.dtypes[lmcc_tf.g1])
                    tc[lmcc_tf.weight] = np.ones(l, dtype=lmcc_tf.dtypes[lmcc_tf.weight]) * weight
                    # tc[lmcc_tf.shape_weight] = np.ones(l, dtype=">f4") / (g1_err**2 + g2_err**2)
                    tc[lmcc_tf.e_var] = np.ones(l, dtype=lmcc_tf.dtypes[lmcc_tf.e_var]) * (g1_err**2 + g2_err**2)
                    tc[lmcc_tf.shape_noise] = np.ones(l, dtype=lmcc_tf.dtypes[lmcc_tf.shape_noise]) * cls.shape_noise
                    tc.add_index(lmcc_tf.ID)
                    cls.chains_table_list.append(tc)

            cls.sem_table_lists[sem] = tables

        return

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):

        self.workdir = tmpdir.strpath
        os.makedirs(self.workdir + "/data")

        return

    @pytest.mark.skip(reason="Not testing at present to save time.")
    def test_reconcile_best(self):

        sem = "LensMC"

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

    @pytest.mark.skip(reason="Not available until DM update")
    def test_reconcile_shape_weight(self):

        sem = "LensMC"

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

        # Check combination of tables 0 and 1 is sensible
        test_row = reconciled_catalog.loc[6]
        assert np.isclose(test_row[sem_tf.g1], self.true_g1[6])
        assert np.isclose(test_row[sem_tf.g2], self.true_g2[6])

        # Weight should be less than the max weight, but higher than at least one individual weight
        assert test_row[sem_tf.weight] < self.max_weight
        assert (test_row[sem_tf.weight] > sem_table_list[0].loc[6][sem_tf.weight] or
                test_row[sem_tf.weight] > sem_table_list[1].loc[6][sem_tf.weight])

        return

    @pytest.mark.skip(reason="Not testing at present to save time.")
    def test_reconcile_invvar(self):

        for sem in sem_names:

            sem_table_list = self.sem_table_lists[sem]
            reconciled_catalog = reconcile_tables(shear_estimates_tables=sem_table_list,
                                                  shear_estimation_method=sem,
                                                  object_ids_in_tile=self.object_ids_in_tile,
                                                  reconciliation_function=reconcile_invvar,
                                                  workdir=self.workdir)

            assert(len(reconciled_catalog) == 19)

            sem_tf = sem_tfs[sem]

            reconciled_catalog.add_index(sem_tf.ID)

            # Row 1 should exactly match results from table 0, since there's no other data for it
            assert_rows_equal(reconciled_catalog, sem_table_list[0], 1, sem_tf)

            # Row 18 should exactly match results from table 2, since there's no other data for it
            assert_rows_equal(reconciled_catalog, sem_table_list[2], 18, sem_tf)

            # Check combination of tables 0 and 1 is sensible
            test_row = reconciled_catalog.loc[6]
            assert np.isclose(test_row[sem_tf.g1], self.true_g1[6])
            assert np.isclose(test_row[sem_tf.g2], self.true_g2[6])

            # Weight should be less than the max weight, but higher than either individual weight
            assert test_row[sem_tf.weight] < self.max_weight
            assert test_row[sem_tf.weight] > sem_table_list[0].loc[6][sem_tf.weight]
            assert test_row[sem_tf.weight] > sem_table_list[1].loc[6][sem_tf.weight]

        return

    @pytest.mark.skip(reason="Not testing at present to save time.")
    def test_reconcile_chains_best(self):

        reconciled_chains = reconcile_chains(chains_tables=self.chains_table_list,
                                             object_ids_in_tile=self.object_ids_in_tile,
                                             chains_reconciliation_function=reconcile_chains_best,
                                             workdir=self.workdir)

        assert(len(reconciled_chains) == 19)

        reconciled_chains.add_index(lmcc_tf.ID)

        # Row 1 should exactly match results from table 0
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 1)

        # Also for 6, which overlaps with table 1, but table 0 has higher weight
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 6)

        # Row 8 should be from table 2, which has the highest weight
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[2], 8)

        return

    @pytest.mark.skip(reason="Not available until DM update")
    def test_reconcile_chains_shape_weight(self):

        reconciled_chains = reconcile_chains(chains_tables=self.chains_table_list,
                                             object_ids_in_tile=self.object_ids_in_tile,
                                             chains_reconciliation_function=reconcile_chains_shape_weight,
                                             workdir=self.workdir)

        assert(len(reconciled_chains) == 19)

        reconciled_chains.add_index(lmcc_tf.ID)

        # Row 1 should exactly match results from table 0, since there's no other data for it
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 1)

        # Row 18 should exactly match results from table 2, since there's no other data for it
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[2], 18)

        # Check combination of tables 0 and 1 is sensible
        test_row = reconciled_chains.loc[6]
        assert np.isclose(test_row[lmcc_tf.g1].mean(), self.true_g1[6])
        assert np.isclose(test_row[lmcc_tf.g2].mean(), self.true_g2[6])

        # Weight should be less than the max weight, but higher than at least one individual weight
        assert test_row[lmcc_tf.weight] < self.max_weight
        assert (test_row[lmcc_tf.weight] > self.chains_table_list[0].loc[6][lmcc_tf.weight] or
                test_row[lmcc_tf.weight] > self.chains_table_list[1].loc[6][lmcc_tf.weight])

        return

    def test_reconcile_chains_invvar(self):

        reconciled_chains = reconcile_chains(chains_tables=self.chains_table_list,
                                             object_ids_in_tile=self.object_ids_in_tile,
                                             chains_reconciliation_function=reconcile_chains_invvar,
                                             workdir=self.workdir)

        assert(len(reconciled_chains) == 19)

        reconciled_chains.add_index(lmcc_tf.ID)

        # Row 1 should exactly match results from table 0, since there's no other data for it
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 1)

        # Row 18 should exactly match results from table 2, since there's no other data for it
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[2], 18)

        # Check combination of tables 0 and 1 is sensible
        test_row = reconciled_chains.loc[6]
        assert np.isclose(test_row[lmcc_tf.g1].mean(), self.true_g1[6])
        assert np.isclose(test_row[lmcc_tf.g2].mean(), self.true_g2[6])

        # Weight should be less than the max weight, but higher than at least one individual weight
        assert test_row[lmcc_tf.weight] < self.max_weight
        assert (test_row[lmcc_tf.weight] > self.chains_table_list[0].loc[6][lmcc_tf.weight] or
                test_row[lmcc_tf.weight] > self.chains_table_list[1].loc[6][lmcc_tf.weight])

        return

    @pytest.mark.skip(reason="Not testing at present to save time.")
    def test_interface(self):
        """Run through the full exectuable, to test the interface.
        """

        # Start by creating data products with the input data

        # Create a data table and product for the MER Final Catalog
        mer_final_catalog_table = initialise_mer_final_catalog()
        for id in self.ids:
            mer_final_catalog_table.add_row()
            mer_final_catalog_table[-1][mfc_tf.ID] = id
        mer_final_catalog_data_filename = get_allowed_filename("MFC", "TEST", version=SHE_CTE.__version__)
        mer_final_catalog_table.write(os.path.join(self.workdir, mer_final_catalog_data_filename))

        mer_final_catalog_product = products.mer_final_catalog.create_dpd_mer_final_catalog(
            mer_final_catalog_data_filename)
        mer_final_catalog_product_filename = get_allowed_filename(
            "MFC-P", "TEST", version=SHE_CTE.__version__, subdir="")
        write_xml_product(mer_final_catalog_product, mer_final_catalog_product_filename,
                          workdir=self.workdir, allow_pickled=False)

        # Create a data product for each input shear estimates table
        sem_product_filename_list = []
        for i in range(len(self.sem_table_lists["KSB"])):
            sem_product = products.she_validated_measurements.create_she_validated_measurements_product()
            for sem in sem_names:
                sem_table_filename = get_allowed_filename("SEM-" + sem.upper(), str(i), version=SHE_CTE.__version__,)
                self.sem_table_lists[sem][i].write(os.path.join(self.workdir, sem_table_filename))
                sem_product.set_method_filename(sem, sem_table_filename)

            sem_product_filename = get_allowed_filename("SEM-P", str(i), version=SHE_CTE.__version__, subdir="",)
            write_xml_product(sem_product, sem_product_filename, workdir=self.workdir)
            sem_product_filename_list.append(sem_product_filename)

        sem_listfile_filename = get_allowed_filename("SEM-L", "TEST", version=SHE_CTE.__version__, subdir="",)
        write_listfile(os.path.join(self.workdir, sem_listfile_filename), sem_product_filename_list)

        # Set up arguments to call the main reconciliation function
        srm_product_filename = get_allowed_filename("SRM-P", "TEST", version=SHE_CTE.__version__, subdir="",)
        args = Args(profile=False,
                    dry_run=False,
                    debug=False,
                    she_validated_measurements_listfile=sem_listfile_filename,
                    mer_final_catalog=mer_final_catalog_product_filename,
                    she_reconciliation_config="None",
                    method=None,
                    she_reconciled_measurements=srm_product_filename,
                    workdir=self.workdir,
                    logdir="logs")

        # Call the program, then check the results
        reconcile_shear_from_args(args)

        srm_product = read_xml_product(srm_product_filename, workdir=self.workdir, allow_pickled=False)

        for sem in sem_names:
            sem_table_filename = srm_product.get_method_filename(sem)
            loaded_sem_table = Table.read(os.path.join(self.workdir, sem_table_filename))
            loaded_sem_table.add_index(sem_tfs[sem].ID)

            # Just a quick test on results, since we do detailed tests elsewhere
            assert_rows_equal(loaded_sem_table, self.sem_table_lists[sem][0], 1, sem_tfs[sem])

        return
