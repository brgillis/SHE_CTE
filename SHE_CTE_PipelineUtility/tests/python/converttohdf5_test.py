#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

"""
:file: tests/python/converttohdf5_test.py

:date: 23/03/2023
:author: Gordon Gibb
"""

import pytest
import os

from ST_DM_DmUtils.DmUtils import read_product_metadata

from SHE_PPT.testing import (
    generate_mock_vis_images,
    generate_mock_reprojected_segmentation_maps,
)

from SHE_PPT.she_io.vis_exposures import VisExposureHDF5, VisExposureAstropyFITS

from SHE_CTE_PipelineUtility.ConvertToHDF5 import mainMethod, defineSpecificProgramOptions

NUM_OBJECTS = 2
OBJSIZE = 2.5


@pytest.fixture
def workdir(tmp_path):
    """Creates a workdir (via pytest's tmpdir fixture) and a datadir"""
    workdir = tmp_path

    # make the datadir
    datadir = os.path.join(workdir, "data")
    os.mkdir(datadir)

    return workdir


@pytest.fixture
def input_products(workdir):
    """Creates the data products needed as an input"""
    exposure_product, sky_coords, pixel_coords, detectors, wcs_list = generate_mock_vis_images.create_exposure(
        workdir=workdir, n_objs_per_det=NUM_OBJECTS, objsize=OBJSIZE, seed=2
    )

    object_ids = range(NUM_OBJECTS)
    seg_map_product = generate_mock_reprojected_segmentation_maps.create_reprojected_segmentation_map(
        object_ids, pixel_coords, detectors, wcs_list, objsize=OBJSIZE, workdir=workdir
    )

    return exposure_product, seg_map_product


class TestConvertToHDF5(object):
    def test_executable(self, workdir, input_products):
        vis_prod, seg_prod = input_products

        output_hdf5 = "out.h5"

        parser = defineSpecificProgramOptions()
        argstring = [
            f"--vis_frame_prod={vis_prod}",
            f"--remapped_seg_prod={seg_prod}",
            f"--output_hdf5={output_hdf5}",
            f"--workdir={workdir}",
            f"--chunksize=50",
        ]
        args = parser.parse_args(argstring)

        # run the executable
        mainMethod(args)

        assert os.path.exists(os.path.join(workdir, output_hdf5)), f"Output file {output_hdf5} does not exist"

        # Now ensure the output is correct (identical data to the input)

        # Parse input data products to get the FITS files
        datadir = os.path.join(workdir, "data")

        vis_dpd = read_product_metadata(os.path.join(workdir, vis_prod))
        seg_dpd = read_product_metadata(os.path.join(workdir, seg_prod))

        det_fits = os.path.join(datadir, vis_dpd.Data.DataStorage.DataContainer.FileName)
        wgt_fits = os.path.join(datadir, vis_dpd.Data.WeightStorage.DataContainer.FileName)
        bkg_fits = os.path.join(datadir, vis_dpd.Data.BackgroundStorage.DataContainer.FileName)
        seg_fits = os.path.join(datadir, seg_dpd.Data.DataStorage.DataContainer.FileName)

        # Create VisExposure objects from the input FITS files, and the output HDF5 file.
        vis_fits = VisExposureAstropyFITS(det_file=det_fits, wgt_file=wgt_fits, bkg_file=bkg_fits, seg_file=seg_fits)
        vis_hdf5 = VisExposureHDF5(os.path.join(workdir, output_hdf5))

        n_vis = len(vis_fits)
        n_hdf5 = len(vis_hdf5)
        assert (
            n_vis == n_hdf5
        ), f"Input FITS exposure has different number of detectors than output HDF5 ({n_vis} vs {n_hdf5})"

        # compare the detector data
        # The detector objects have an __eq__ dunder method defined that checks the data contents are identical
        for i in range(n_vis):
            det_fits = vis_fits[i]
            det_hdf5 = vis_hdf5[i]
            assert det_fits == det_hdf5, "Detector from the input FITS differs from the detector from the output HDF5"

