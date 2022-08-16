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
from copy import deepcopy

import numpy as np
import pytest
from astropy.table import Row, Table

import SHE_CTE
from SHE_CTE_ShearReconciliation.ReconcileShear import D_REC_SHEAR_CONFIG_DEFAULTS, D_REC_SHEAR_CONFIG_TYPES
from SHE_CTE_ShearReconciliation.chains_reconciliation_functions import (reconcile_chains_best, reconcile_chains_invvar,
                                                                         reconcile_chains_keep,
                                                                         reconcile_chains_shape_weight, )
from SHE_CTE_ShearReconciliation.reconcile_shear import reconcile_chains, reconcile_shear_from_args, reconcile_tables
from SHE_CTE_ShearReconciliation.reconciliation_functions import (reconcile_best, reconcile_invvar,
                                                                  reconcile_shape_weight, )
from SHE_PPT import products
from SHE_PPT.constants.shear_estimation_methods import D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS, ShearEstimationMethods
from SHE_PPT.file_io import get_allowed_filename, read_xml_product, write_listfile, write_xml_product
from SHE_PPT.pipeline_utility import ReconciliationConfigKeys, read_config
from SHE_PPT.table_formats.mer_final_catalog import initialise_mer_final_catalog, tf as mfc_tf
from SHE_PPT.table_formats.she_lensmc_chains import initialise_lensmc_chains_table, len_chain, tf as lmcc_tf


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


def assert_chains_rows_equal(t1, t2, i, which_of_kept = None):
    """ Check that two chains rows match in shear parameters
    """

    r1 = t1.loc[i]
    if which_of_kept is not None:
        if (isinstance(r1, Row) and which_of_kept != 0) or which_of_kept > len(r1):
            raise ValueError("Fewer rows kept (" + str(len(r1)) +
                             ") than requested index (" + str(which_of_kept) + ").")
        if not isinstance(r1, Row):
            r1 = r1[which_of_kept]
    r2 = t2.loc[i]
    assert np.isclose(r1[lmcc_tf.g1], r2[lmcc_tf.g1]).all(), "Row " + str(i) + " doesn't match expected value for g1."
    assert np.isclose(r1[lmcc_tf.g2], r2[lmcc_tf.g2]).all(), "Row " + str(i) + " doesn't match expected value for g2."
    assert np.isclose(r1[lmcc_tf.shape_noise], r2[lmcc_tf.shape_noise]), "Row " + \
                                                                         str(i) + " doesn't match expected value for " \
                                                                                  "shape noise."
    assert np.isclose(r1[lmcc_tf.weight], r2[lmcc_tf.weight]), "Row " + \
                                                               str(i) + " doesn't match expected value for weight."


class TestReconcileShear:
    """
    """

    @classmethod
    def setup_class(cls):

        # Set up a mock catalogue with 20 objects, but say only the first 19 are in the tile
        cls.ids = np.arange(19, dtype = int)
        cls.object_ids_in_tile = frozenset(cls.ids)

        # Set up "true" values for g1 and g2, which each other table will be biased from
        cls.true_g1 = np.linspace(-0.8, 0.8, num = 20, dtype = ">f4")
        cls.true_g2 = np.where(cls.true_g1 < 0, 0.8 + cls.true_g1, -0.8 + cls.true_g1)

        cls.shape_noise = 0.25
        cls.max_weight = 1 / cls.shape_noise ** 2

        # We'll set up some mock tables from each method, using the same values for each method
        cls.sem_table_lists = {}
        cls.chains_table_list = []
        for sem in ShearEstimationMethods:

            tf = D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS[sem]

            make_chains = sem == ShearEstimationMethods.LENSMC

            tables = []

            # In brief - Table 1 and 2 cancel out for invvar weighting, but not shape weighting. Table 3 dominates,
            # Table 4 has NaN errors and 0 estimates, Table 5 has NaN estimates and errors and Inf weight,
            # Table 6 has impossible (negative) values for weight and errors
            for i_min, i_max, g1_offset, g2_offset, e1_err, e2_err, weight in ((0, 9, 0.01, -0.01, 0.01, 0.01, None),
                                                                               (6, 12, -0.04, 0.04, 0.02, 0.02, None),
                                                                               (8, 19, 0.1, 0.2, 0.001, 0.001, None),
                                                                               (0, 19, 0, 0, np.nan, np.nan, np.nan),
                                                                               (0, 19, np.nan, np.nan,
                                                                                np.nan, np.nan, np.inf),
                                                                               (0, 19, 0, 0, -1, -np.inf, -1),):

                if weight is None:
                    weight = 1 / (0.5 * (e1_err ** 2 + e2_err ** 2) + cls.shape_noise ** 2)
                g1_err = np.sqrt(e1_err ** 2 + cls.shape_noise ** 2)
                g2_err = np.sqrt(e2_err ** 2 + cls.shape_noise ** 2)

                l = i_max - i_min + 1

                # Create the measurements table

                t = tf.init_table(optional_columns = (tf.shape_weight, tf.e_var, tf.shape_noise))
                for _ in range(l):
                    t.add_row()

                t[tf.ID] = np.arange(l) + i_min
                t[tf.g1] = cls.true_g1[i_min:i_max + 1] + g1_offset
                t[tf.g2] = cls.true_g2[i_min:i_max + 1] + g2_offset
                t[tf.e1_err] = np.ones(l, dtype = tf.dtypes[tf.g1_err]) * e1_err
                t[tf.e2_err] = np.ones(l, dtype = tf.dtypes[tf.g2_err]) * e2_err
                t[tf.g1_err] = np.ones(l, dtype = tf.dtypes[tf.g1_err]) * g1_err
                t[tf.g2_err] = np.ones(l, dtype = tf.dtypes[tf.g2_err]) * g2_err
                t[tf.weight] = np.ones(l, dtype = tf.dtypes[tf.weight]) * weight
                t[tf.shape_weight] = np.ones(l, dtype = tf.dtypes[tf.shape_weight]) / (e1_err ** 2 + e2_err ** 2)
                t[tf.e_var] = np.ones(l, dtype = tf.dtypes[tf.e_var]) * (e1_err ** 2 + e2_err ** 2)
                t[tf.shape_noise] = np.ones(l, dtype = tf.dtypes[tf.shape_noise]) * cls.shape_noise
                t.add_index(tf.ID)
                tables.append(t)

                if make_chains:
                    # Create the chains table

                    tc = initialise_lensmc_chains_table(optional_columns = (lmcc_tf.e_var, lmcc_tf.shape_noise))
                    for _ in range(l):
                        tc.add_row()

                    # Use a stable set of deviates which average to zero for consistency
                    deviates_first_half = np.random.standard_normal((l, len_chain // 2))
                    deviates_second_half = -deviates_first_half
                    deviates = np.concatenate((deviates_first_half, deviates_second_half), axis = 1)

                    tc[lmcc_tf.ID] = np.arange(l) + i_min
                    tc[lmcc_tf.g1] = (np.array((cls.true_g1[i_min:i_max + 1],)).transpose() + g1_offset +
                                      deviates * g1_err).astype(lmcc_tf.dtypes[lmcc_tf.g1])
                    tc[lmcc_tf.g2] = (np.array((cls.true_g2[i_min:i_max + 1],)).transpose() + g2_offset +
                                      deviates * g2_err).astype(lmcc_tf.dtypes[lmcc_tf.g1])
                    tc[lmcc_tf.weight] = np.ones(l, dtype = lmcc_tf.dtypes[lmcc_tf.weight]) * weight
                    tc[lmcc_tf.shape_weight] = np.ones(
                        l, dtype = lmcc_tf.dtypes[lmcc_tf.shape_weight]) * 2 / (e1_err ** 2 + e2_err ** 2)
                    tc[lmcc_tf.e_var] = np.ones(l, dtype = lmcc_tf.dtypes[lmcc_tf.e_var]) * (e1_err ** 2 + e2_err ** 2)
                    tc[lmcc_tf.shape_noise] = np.ones(l, dtype = lmcc_tf.dtypes[lmcc_tf.shape_noise]) * cls.shape_noise
                    tc.add_index(lmcc_tf.ID)
                    cls.chains_table_list.append(tc)

            cls.sem_table_lists[sem] = tables

    @classmethod
    def teardown_class(cls):

        return

    @pytest.fixture(autouse = True)
    def setup(self, tmpdir):

        self.workdir = tmpdir.strpath
        os.makedirs(self.workdir + "/data")

    def test_reconcile_best(self):

        sem = ShearEstimationMethods.LENSMC

        sem_table_list = self.sem_table_lists[sem]
        reconciled_catalog = reconcile_tables(shear_estimates_tables = sem_table_list,
                                              shear_estimation_method = sem,
                                              object_ids_in_tile = self.object_ids_in_tile,
                                              reconciliation_function = reconcile_best,
                                              workdir = self.workdir)

        assert (len(reconciled_catalog) == 19)

        sem_tf = D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS[sem]

        reconciled_catalog.add_index(sem_tf.ID)

        # Row 1 should exactly match results from table 0
        assert_rows_equal(reconciled_catalog, sem_table_list[0], 1, sem_tf)

        # Also for 6, which overlaps with table 1, but table 0 has higher weight
        assert_rows_equal(reconciled_catalog, sem_table_list[0], 6, sem_tf)

        # Row 8 should be from table 2, which has the highest weight
        assert_rows_equal(reconciled_catalog, sem_table_list[2], 8, sem_tf)

    @pytest.mark.skip(reason = "Not available until DM update")
    def test_reconcile_shape_weight(self):

        sem = ShearEstimationMethods.LENSMC

        sem_table_list = self.sem_table_lists[sem]
        reconciled_catalog = reconcile_tables(shear_estimates_tables = sem_table_list,
                                              shear_estimation_method = sem,
                                              object_ids_in_tile = self.object_ids_in_tile,
                                              reconciliation_function = reconcile_shape_weight,
                                              workdir = self.workdir)

        assert (len(reconciled_catalog) == 19)

        sem_tf = D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS[sem]

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

    def test_reconcile_invvar(self):

        for sem in ShearEstimationMethods:

            sem_tf = D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS[sem]

            sem_table_list = self.sem_table_lists[sem]

            reconciled_catalog = reconcile_tables(shear_estimates_tables = sem_table_list,
                                                  shear_estimation_method = sem,
                                                  object_ids_in_tile = self.object_ids_in_tile,
                                                  reconciliation_function = reconcile_invvar,
                                                  workdir = self.workdir)

            assert (len(reconciled_catalog) == 19)

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

            # Try some extra tables with varying columns, just for KSB to save time
            if sem == "KSB":
                t_alt1 = deepcopy(sem_table_list[0])
                t_alt1.remove_column(sem_tf.shape_weight)
                sem_table_list_alt = [t_alt1, t_alt1]
                _ = reconcile_tables(shear_estimates_tables = sem_table_list_alt,
                                     shear_estimation_method = sem,
                                     object_ids_in_tile = self.object_ids_in_tile,
                                     reconciliation_function = reconcile_invvar,
                                     workdir = self.workdir)

                t_alt2 = deepcopy(t_alt1)
                t_alt2.remove_column(sem_tf.e_var)
                sem_table_list_alt = [t_alt2, t_alt2]
                _ = reconcile_tables(shear_estimates_tables = sem_table_list_alt,
                                     shear_estimation_method = sem,
                                     object_ids_in_tile = self.object_ids_in_tile,
                                     reconciliation_function = reconcile_invvar,
                                     workdir = self.workdir)

                t_alt3 = deepcopy(t_alt2)
                t_alt3[sem_tf.e1_err] = 0.
                sem_table_list_alt = [t_alt3, t_alt3]
                _ = reconcile_tables(shear_estimates_tables = sem_table_list_alt,
                                     shear_estimation_method = sem,
                                     object_ids_in_tile = self.object_ids_in_tile,
                                     reconciliation_function = reconcile_invvar,
                                     workdir = self.workdir)

    def test_reconcile_chains_best(self):

        reconciled_chains = reconcile_chains(chains_tables = self.chains_table_list,
                                             object_ids_in_tile = self.object_ids_in_tile,
                                             chains_reconciliation_function = reconcile_chains_best,
                                             workdir = self.workdir)

        assert (len(reconciled_chains) == 19)

        reconciled_chains.add_index(lmcc_tf.ID)

        # Row 1 should exactly match results from table 0
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 1)

        # Also for 6, which overlaps with table 1, but table 0 has higher weight
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 6)

        # Row 8 should be from table 2, which has the highest weight
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[2], 8)

    @pytest.mark.skip(reason = "Not available until DM update")
    def test_reconcile_chains_shape_weight(self):

        reconciled_chains = reconcile_chains(chains_tables = self.chains_table_list,
                                             object_ids_in_tile = self.object_ids_in_tile,
                                             chains_reconciliation_function = reconcile_chains_shape_weight,
                                             workdir = self.workdir)

        assert (len(reconciled_chains) == 19)

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

    def test_reconcile_chains_invvar(self):

        reconciled_chains = reconcile_chains(chains_tables = self.chains_table_list,
                                             object_ids_in_tile = self.object_ids_in_tile,
                                             chains_reconciliation_function = reconcile_chains_invvar,
                                             workdir = self.workdir)

        assert (len(reconciled_chains) == 19)

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

    def test_reconcile_chains_keep(self):

        reconciled_chains = reconcile_chains(chains_tables = self.chains_table_list,
                                             object_ids_in_tile = self.object_ids_in_tile,
                                             chains_reconciliation_function = reconcile_chains_keep,
                                             workdir = self.workdir)

        assert (len(reconciled_chains) == 85)

        reconciled_chains.add_index(lmcc_tf.ID)

        # Row 1 should exactly match results from table 0, since there's no other data for it
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 1, 0)

        # Row 18 should exactly match results from table 2, since there's no other data for it
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[2], 18, 0)

        # Check combination of tables 0 and 1 is sensible
        test_row_0 = reconciled_chains.loc[6][0]
        test_row_1 = reconciled_chains.loc[6][1]
        assert test_row_0[lmcc_tf.g1].mean() > self.true_g1[6]
        assert test_row_0[lmcc_tf.g2].mean() < self.true_g2[6]
        assert test_row_1[lmcc_tf.g1].mean() < self.true_g1[6]
        assert test_row_1[lmcc_tf.g2].mean() > self.true_g2[6]

        # Test that the weights add up properly
        assert test_row_0[lmcc_tf.weight] + test_row_1[lmcc_tf.weight] < self.max_weight

    def test_interface(self):
        """Run through the full executable, to test the interface.
        """

        # Start by creating data products with the input data

        # Create a data table and product for the MER Final Catalog
        mer_final_catalog_table = initialise_mer_final_catalog()
        for object_id in self.ids:
            mer_final_catalog_table.add_row()
            mer_final_catalog_table[-1][mfc_tf.ID] = object_id
        mer_final_catalog_data_filename = get_allowed_filename("MFC", "TEST", version = SHE_CTE.__version__)
        mer_final_catalog_table.write(os.path.join(self.workdir, mer_final_catalog_data_filename))

        mer_final_catalog_product = products.mer_final_catalog.create_dpd_mer_final_catalog(
            mer_final_catalog_data_filename)
        mer_final_catalog_product_filename = get_allowed_filename(
            "MFC-P", "TEST", version = SHE_CTE.__version__, subdir = "")
        write_xml_product(mer_final_catalog_product, mer_final_catalog_product_filename,
                          workdir = self.workdir, allow_pickled = False)

        # Create a data product for each input shear estimates table, with overlapping pointing IDs
        sem_product_filename_list = []
        num_products = len(self.sem_table_lists[ShearEstimationMethods.LENSMC])

        for i in range(num_products):

            # Make a separate list of mock Pointing IDs for each product - each will have an ID of 100 (overlapping),
            # and the rest will be the index of this product
            l_pointing_ids = [100, i]

            sem_product = products.she_validated_measurements.create_she_validated_measurements_product()
            sem_product.Data.PointingIdList = l_pointing_ids

            for sem in ShearEstimationMethods:
                sem_table_filename = get_allowed_filename(
                    "SEM-TEST-" + sem.name, str(i), version = SHE_CTE.__version__, )
                self.sem_table_lists[sem][i].write(os.path.join(self.workdir, sem_table_filename))
                sem_product.set_method_filename(sem, sem_table_filename)

            sem_product_filename = get_allowed_filename(
                "SEM-TEST-P", str(i), version = SHE_CTE.__version__, subdir = "", extension = ".xml", )
            write_xml_product(sem_product, sem_product_filename, workdir = self.workdir)
            sem_product_filename_list.append(sem_product_filename)

        sem_listfile_filename = get_allowed_filename(
            "SEM-TEST-L", "0", version = SHE_CTE.__version__, subdir = "", extension = ".json", )
        write_listfile(os.path.join(self.workdir, sem_listfile_filename), sem_product_filename_list)

        # Create a data product for each input chains table, and a listfile of their names
        chains_product_filename_list = []

        for i in range(len(self.chains_table_list)):
            chains_product = products.she_lensmc_chains.create_dpd_she_lensmc_chains()

            chains_table_filename = get_allowed_filename("CHAINS-TEST", str(i), version = SHE_CTE.__version__, )
            self.chains_table_list[i].write(os.path.join(self.workdir, chains_table_filename))
            chains_product.set_filename(chains_table_filename)

            chains_product_filename = get_allowed_filename(
                "CHAINS-TEST-P", str(i), version = SHE_CTE.__version__, subdir = "", extension = ".xml")
            write_xml_product(chains_product, chains_product_filename, workdir = self.workdir)
            chains_product_filename_list.append(chains_product_filename)

        chains_listfile_filename = get_allowed_filename(
            "CHAINS-TEST-L", "0", version = SHE_CTE.__version__, subdir = "", extension = ".json", )
        write_listfile(os.path.join(self.workdir, chains_listfile_filename), chains_product_filename_list)

        # Get desired filenames for output products
        srm_product_filename = get_allowed_filename(
            "SRM-P", "TEST", version = SHE_CTE.__version__, subdir = "", extension = ".xml", )
        rec_chains_product_filename = get_allowed_filename(
            "REC-CHAINS-TEST-P", "0", version = SHE_CTE.__version__, subdir = "", extension = ".xml", )

        # Set up arguments to call the main reconciliation function
        args = Args(profile = False,
                    dry_run = False,
                    debug = False,
                    she_validated_measurements_listfile = sem_listfile_filename,
                    she_lensmc_chains_listfile = chains_listfile_filename,
                    mer_final_catalog = mer_final_catalog_product_filename,
                    she_reconciliation_config = read_config(None, config_keys = ReconciliationConfigKeys,
                                                            d_defaults = D_REC_SHEAR_CONFIG_DEFAULTS,
                                                            d_types = D_REC_SHEAR_CONFIG_TYPES),
                    method = None,
                    chains_method = None,
                    she_reconciled_measurements = srm_product_filename,
                    she_reconciled_lensmc_chains = rec_chains_product_filename,
                    workdir = self.workdir,
                    logdir = "logs")

        # Call the program
        reconcile_shear_from_args(args)

        # Check the reconciled measurements results
        srm_product = read_xml_product(srm_product_filename, workdir = self.workdir, allow_pickled = False)

        # Check the Pointing ID list is as expected - one of each product index, plus 100
        assert srm_product.Data.PointingIdList == list(range(num_products)) + [100]

        for sem in ShearEstimationMethods:
            sem_table_filename = srm_product.get_method_filename(sem)
            loaded_sem_table = Table.read(os.path.join(self.workdir, sem_table_filename))
            loaded_sem_table.add_index(D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS[sem].ID)

            # Just a quick test on results, since we do detailed tests elsewhere
            assert_rows_equal(loaded_sem_table, self.sem_table_lists[sem]
            [0], 1, D_SHEAR_ESTIMATION_METHOD_TABLE_FORMATS[sem])

        # Check the reconciled chains results
        rec_chains_product = read_xml_product(rec_chains_product_filename, workdir = self.workdir,
                                              allow_pickled = False)

        rec_chains_table_filename = rec_chains_product.get_filename()
        reconciled_chains = Table.read(os.path.join(self.workdir, rec_chains_table_filename))
        reconciled_chains.add_index(lmcc_tf.ID)

        # Just a quick test on results, since we do detailed tests elsewhere
        assert_chains_rows_equal(reconciled_chains, self.chains_table_list[0], 1, 0)
