""" @file estimate_shear_interface_test.py

    Created 6 August 2018

    Unit tests for the control shear estimation methods.
"""

__updated__ = "2020-07-10"

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

from SHE_PPT import mdb
from SHE_PPT.file_io import read_pickled_product, find_file
from SHE_PPT.logging import getLogger
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksbm_tf
from SHE_PPT.table_formats.she_regauss_measurements import tf as regm_tf
import pytest

from ElementsServices.DataSync import DataSync
from SHE_CTE_ShearEstimation.control_training_data import load_control_training_data
from SHE_CTE_ShearEstimation.galsim_estimate_shear import (KSB_estimate_shear, REGAUSS_estimate_shear)
import SHE_LensMC.she_measure_shear
# from SHE_MomentsML.estimate_shear import estimate_shear as ML_estimate_shear # TODO: Uncomment when MomentsML is updated to EDEN 2.1
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


class TestCase:
    """


    """

    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):

        # Download the MDB from WebDAV
        sync_mdb = DataSync("testdata/sync.conf", "testdata/test_mdb.txt")
        sync_mdb.download()
        mdb_filename = sync_mdb.absolutePath("SHE_CTE_8_1/sample_mdb.xml")

        mdb.init(mdb_filename)

        # Download the training data files from WebDAV
        sync_training = DataSync("testdata/sync.conf", "testdata/test_control_training_data.txt")
        sync_training.download()

        # Download the data stack files from WebDAV
        sync_datastack = DataSync("testdata/sync.conf", "testdata/test_data_stack.txt")
        sync_datastack.download()
        qualified_data_images_filename = sync_datastack.absolutePath("SHE_CTE_8_1/vis_calibrated_frames.json")
        assert os.path.isfile(qualified_data_images_filename), f"Cannot find file: {qualified_data_images_filename}"

        # Get the workdir based on where the data images listfile is
        self.workdir = os.path.split(qualified_data_images_filename)[0]
        self.logdir = os.path.join(self.workdir, "logs")

        # Read in the test data
        self.data_stack = SHEFrameStack.read(exposure_listfile_filename=data_images_filename,
                                             seg_listfile_filename=segmentation_images_filename,
                                             stacked_image_product_filename=stacked_image_filename,
                                             stacked_seg_product_filename=stacked_segmentation_image_filename,
                                             psf_listfile_filename=psf_images_and_tables_filename,
                                             detections_listfile_filename=detections_tables_filename,
                                             workdir=self.workdir,
                                             clean_detections=True,
                                             memmap=True,
                                             mode='denywrite')

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
                    raise Exception("Bad value for " + colname + ": " + str(g))

        return

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
                    raise Exception("Bad value for " + colname + ": " + str(g))

        return

    @pytest.mark.skip(reason="Too extensive for a unit test - should be turned into a smoke test")
    def test_momentsml(self):
        """Test that the interface for the MomentsML method works properly.
        """

        # Read in the test data
        she_frame = read_pickled_product(find_file(she_frame_location))
        original_she_frame = deepcopy(she_frame)

        momentsml_cat = ML_estimate_shear(she_frame,
                                          training_data=None,
                                          calibration_data=None,
                                          workdir=self.workdir)

        # Check that we have valid data
        for row in momentsml_cat:
            for colname in (sm_tf.g1, sm_tf.g2, sm_tf.g1_err, sm_tf.g2_err):
                g = row[colname]
                if not (g > -1 and g < 1):
                    raise Exception("Bad value for " + colname + ": " + str(g))

        # Check that the input data isn't changed by the method
        assert she_frame == original_she_frame

        return
