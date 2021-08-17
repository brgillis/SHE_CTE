""" @file estimate_shear_interface_test.py

    Created 6 August 2018

    Unit tests for the control shear estimation methods.
"""

__updated__ = "2021-08-17"

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
import time

import SHE_LensMC.she_measure_shear
from SHE_MomentsML.estimate_shear import estimate_shear as ML_estimate_shear
from SHE_PPT import mdb
from SHE_PPT.file_io import read_xml_product, find_file, read_listfile
from SHE_PPT.logging import getLogger
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksbm_tf
from SHE_PPT.table_formats.she_lensmc_chains import tf as lmcc_tf, len_chain
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_formats.she_regauss_measurements import tf as regm_tf
import pytest

from ElementsServices.DataSync import DataSync
from SHE_CTE_ShearEstimation.control_training_data import load_control_training_data
from SHE_CTE_ShearEstimation.estimate_shears import fill_measurements_table_meta
from SHE_CTE_ShearEstimation.galsim_estimate_shear import (KSB_estimate_shear, REGAUSS_estimate_shear)
import numpy as np


test_data_location = "/tmp"

data_images_filename = "vis_calibrated_frames.json"
segmentation_images_filename = "she_exposure_reprojected_segmentation_maps.json"
stacked_image_filename = "vis_stacked_image.xml"
stacked_segmentation_image_filename = "she_stack_reprojected_segmentation_map.xml"
psf_images_and_tables_filename = "she_psf_model_images.json"
detections_tables_filename = "mer_final_catalogs.json"
ksb_training_filename = "test_ksb_training.xml"
regauss_training_filename = "test_regauss_training.xml"

expected_observation_id = -1
expected_observation_time = "2020-06-10 15:00:36.660000+00:00"
expected_pointing_id_list = "-1 -1 -1 -1"
expected_tile_id = "1"


class TestCase:
    """


    """

    @classmethod
    def setup_class(cls):

        # Download the MDB from WebDAV
        sync_mdb = DataSync("testdata/sync.conf", "testdata/test_mdb.txt")
        sync_mdb.download()
        mdb_filename = sync_mdb.absolutePath("SHE_PPT_8_7/sample_mdb-SC8.xml")

        mdb.init(mdb_filename)

        # Download the training data files from WebDAV
        sync_training = DataSync("testdata/sync.conf", "testdata/test_control_training_data.txt")
        sync_training.download()

        # Download the data stack files from WebDAV
        sync_datastack = DataSync("testdata/sync.conf", "testdata/test_data_stack.txt")
        sync_datastack.download()
        qualified_data_images_filename = sync_datastack.absolutePath("SHE_PPT_8_7/vis_calibrated_frames.json")
        assert os.path.isfile(qualified_data_images_filename), f"Cannot find file: {qualified_data_images_filename}"

        # Get the workdir based on where the data images listfile is
        cls.workdir = os.path.split(qualified_data_images_filename)[0]
        cls.logdir = os.path.join(cls.workdir, "logs")

        # Read in the test data
        cls.data_stack = SHEFrameStack.read(exposure_listfile_filename=data_images_filename,
                                            seg_listfile_filename=segmentation_images_filename,
                                            stacked_image_product_filename=stacked_image_filename,
                                            stacked_seg_product_filename=stacked_segmentation_image_filename,
                                            psf_listfile_filename=psf_images_and_tables_filename,
                                            detections_listfile_filename=detections_tables_filename,
                                            workdir=cls.workdir,
                                            clean_detections=True,
                                            memmap=True,
                                            mode='denywrite')

    @classmethod
    def teardown_class(cls):

        return

    def test_ksb(self):
        """Test that the interface for the KSB method works properly.
        """

        ksb_training_data = load_control_training_data(find_file(ksb_training_filename, self.workdir),
                                                       workdir=self.workdir)

        ksb_cat = KSB_estimate_shear(self.data_stack,
                                     training_data=ksb_training_data,
                                     calibration_data=None,
                                     workdir=self.workdir)

        # Check that we have valid data
        for row in ksb_cat:
            for colname in (ksbm_tf.g1, ksbm_tf.g2, ksbm_tf.g1_err, ksbm_tf.g2_err):
                g = row[colname]
                if not (g > -1 and g < 1):
                    raise AssertionError("Bad value for " + colname + ": " + str(g))

    def test_regauss(self):
        """Test that the interface for the REGAUSS method works properly.
        """

        regauss_training_data = load_control_training_data(find_file(regauss_training_filename, self.workdir),
                                                           workdir=self.workdir)

        regauss_cat = REGAUSS_estimate_shear(self.data_stack,
                                             training_data=regauss_training_data,
                                             calibration_data=None,
                                             workdir=self.workdir)

        # Check that we have valid data
        for row in regauss_cat:
            for colname in (regm_tf.g1, regm_tf.g2, regm_tf.g1_err, regm_tf.g2_err):
                g = row[colname]
                if not (g > -1 and g < 1):
                    raise AssertionError("Bad value for " + colname + ": " + str(g))

    def test_return_chains(self):

        ksb_training_data = load_control_training_data(find_file(ksb_training_filename, self.workdir),
                                                       workdir=self.workdir)

        # Seed the random number generator for consistent test results
        np.random.seed(1234)

        ksb_cat, ksb_chains = KSB_estimate_shear(self.data_stack,
                                                 training_data=ksb_training_data,
                                                 calibration_data=None,
                                                 workdir=self.workdir,
                                                 return_chains=True)

        # Check that we have valid data
        assert len(ksb_cat) == len(ksb_chains)

        for i in range(len(ksb_cat)):

            estimates_row = ksb_cat[i]
            chains_row = ksb_chains[i]

            for param in ("g1", "g2", "ra", "dec"):

                vals = chains_row[getattr(lmcc_tf, param)]
                ex_mean = estimates_row[getattr(ksbm_tf, param)]
                ex_stddev = estimates_row[getattr(ksbm_tf, param + "_err")]

                assert len(vals) == len_chain
                assert len(vals) == ksb_chains.meta[lmcc_tf.m.len_chain]

                mean_val = np.mean(vals)
                stddev_val = np.std(vals)

                assert np.isclose(mean_val, ex_mean, rtol=0, atol=1e-8 + 2 * ex_stddev /
                                  np.sqrt(len_chain)), "Error matching mean for parameter " + param
                assert np.isclose(stddev_val, ex_stddev, rtol=0, atol=1e-8 + 2 * ex_stddev /
                                  np.sqrt(len_chain)), "Error matching stddev for parameter " + param

    def test_fill_measurements_table_meta(self):
        """Test that the fill_measurements_table_meta utility function works as expected.
        """
        mer_final_catalog_products = []
        for mer_final_catalog_filename in read_listfile(os.path.join(self.workdir, detections_tables_filename)):
            mer_final_catalog_products.append(read_xml_product(os.path.join(self.workdir, mer_final_catalog_filename)))

        vis_calibrated_frame_products = []
        for vis_calibrated_frame_filename in read_listfile(os.path.join(self.workdir, data_images_filename)):
            vis_calibrated_frame_products.append(read_xml_product(
                os.path.join(self.workdir, vis_calibrated_frame_filename)))

        t = ksbm_tf.init_table()

        fill_measurements_table_meta(t=t,
                                     mer_final_catalog_products=mer_final_catalog_products,
                                     vis_calibrated_frame_products=vis_calibrated_frame_products)

        assert t.meta[sm_tf.m.observation_id] == expected_observation_id
        assert t.meta[sm_tf.m.observation_time] == expected_observation_time
        assert t.meta[sm_tf.m.pointing_id] == expected_pointing_id_list
        assert t.meta[sm_tf.m.tile_id] == expected_tile_id
