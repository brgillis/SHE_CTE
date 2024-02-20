""" @file reorder_listfile_test.py

    Created 20 Feb 2024

    Tests for SHE_CTE_ReorderListfile
"""

__updated__ = "2024-02-20"

# Copyright (C) 2012-2022 Euclid Science Ground Segment
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


import pytest
import json

import numpy as np

from ST_DM_FilenameProvider.FilenameProvider import FileNameProvider
from ST_DM_DmUtils import DmUtils

from SHE_PPT.products import she_psf_field_parameters, vis_calibrated_frame, she_exposure_segmentation_map

from SHE_CTE_PipelineUtility.ReorderListfile import mainMethod, defineSpecificProgramOptions, rgetattr


POINTINGS = [1, 2, 3, 4, 5]
OBSERVATIONS = [1, 1, 2, 3, 3]


@pytest.fixture
def workdir(tmp_path):
    workdir = tmp_path
    (tmp_path / "data").mkdir()
    return workdir


@pytest.fixture
def vis_products(workdir):
    """Creates vis products with observation and pointing ids"""

    files = []
    for pnt, obs in zip(POINTINGS, OBSERVATIONS):
        filename = FileNameProvider().get_allowed_filename(
            type_name="VIS-FRAME",
            instance_id=f"{obs}-{pnt}",
            extension=".xml",
            release="0.0",
            processing_function="SHE",
        )

        p = vis_calibrated_frame.create_dpd_vis_calibrated_frame()
        p.Data.ObservationSequence.ObservationId = obs
        p.Data.ObservationSequence.PointingId = pnt

        DmUtils.save_product_metadata(p, workdir / filename)
        files.append(filename)

    listfile_filename = "vis_products.json"

    with open(workdir / listfile_filename, "w") as f:
        json.dump(files, f)

    return listfile_filename


@pytest.fixture
def seg_products(workdir):
    """Creates seg products with reversed obs and pointing ids"""

    files = []
    for pnt, obs in zip(reversed(POINTINGS), reversed(OBSERVATIONS)):
        filename = FileNameProvider().get_allowed_filename(
            type_name="SEG",
            instance_id=f"{obs}-{pnt}",
            extension=".xml",
            release="0.0",
            processing_function="SHE",
        )

        p = she_exposure_segmentation_map.create_dpd_she_exposure_segmentation_map()
        p.Data.ObservationId = obs
        p.Data.PointingId = pnt

        DmUtils.save_product_metadata(p, workdir / filename)
        files.append(filename)

    listfile_filename = "seg_products.json"

    with open(workdir / listfile_filename, "w") as f:
        json.dump(files, f)

    return listfile_filename


@pytest.fixture
def psf_products(workdir):
    """Creates a listfile of PSF field parameters products containing only the unique obs ids"""

    unique_obs, inds = np.unique(OBSERVATIONS, return_index=True)
    unique_pointings = [POINTINGS[i] for i in inds]

    files = []
    for pnt, obs in zip(reversed(unique_pointings), reversed(list(unique_obs))):
        filename = FileNameProvider().get_allowed_filename(
            type_name="PSF-FIELD-PARAMS",
            instance_id=f"{obs}-{pnt}",
            extension=".xml",
            release="0.0",
            processing_function="SHE",
        )

        p = she_psf_field_parameters.create_dpd_she_psf_field_parameters()
        p.Data.ObservationId = obs
        p.Data.PointingId = pnt

        DmUtils.save_product_metadata(p, workdir / filename)
        files.append(filename)

    listfile_filename = "psf_products.json"

    with open(workdir / listfile_filename, "w") as f:
        json.dump(files, f)

    return listfile_filename


def test_basic_reordering(workdir, vis_products, seg_products):
    """Checks that we can reorder segmentation map products with the VIS products"""

    output_products = "output.json"

    parser = defineSpecificProgramOptions()

    argstring = [
        f"--workdir={workdir}",
        f"--reference_listfile={vis_products}",
        f"--input_listfile={seg_products}",
        f"--output_listfile={output_products}",
        "--reference_path=Data.ObservationSequence.PointingId",
        "--input_path=Data.PointingId",
    ]

    args = parser.parse_args(argstring)

    mainMethod(args)

    with open(workdir / output_products) as f:
        prods = [DmUtils.read_product_metadata(workdir / p) for p in json.load(f)]

    pointings = [rgetattr(p, "Data.PointingId") for p in prods]
    observations = [rgetattr(p, "Data.ObservationId") for p in prods]

    assert pointings == POINTINGS
    assert observations == OBSERVATIONS


def test_invalid_reordering_reference(workdir, vis_products, seg_products):
    """Test that the reordering fails if there are multiple identical entries in the reference listfile"""

    output_products = "output.json"

    parser = defineSpecificProgramOptions()

    argstring = [
        f"--workdir={workdir}",
        f"--reference_listfile={vis_products}",
        f"--input_listfile={seg_products}",
        f"--output_listfile={output_products}",
        "--reference_path=Data.ObservationSequence.ObservationId",
        "--input_path=Data.ObservationId",
    ]

    args = parser.parse_args(argstring)

    with pytest.raises(ValueError, match="There are duplicate reference values"):
        mainMethod(args)


def test_invalid_reordering_input(workdir, vis_products, seg_products):
    """Tests that the code correctly throws an error if the input products have duplicate values"""

    output_products = "output.json"

    parser = defineSpecificProgramOptions()

    # This is being called in a way we wouldn't want to call, matching obs id to pointing id, but the code should
    # flag that the seg product's obs ids are not unique
    argstring = [
        f"--workdir={workdir}",
        f"--reference_listfile={vis_products}",
        f"--input_listfile={seg_products}",
        f"--output_listfile={output_products}",
        "--reference_path=Data.ObservationSequence.PointingId",
        "--input_path=Data.ObservationId",
    ]

    args = parser.parse_args(argstring)

    with pytest.raises(ValueError, match="Some of the input values are duplicates"):
        mainMethod(args)


def test_duplicate_reordering(workdir, vis_products, psf_products):
    """Checks that the executable can reorder listfiles with duplicate reference values"""

    output_products = "output.json"

    parser = defineSpecificProgramOptions()

    argstring = [
        f"--workdir={workdir}",
        f"--reference_listfile={vis_products}",
        f"--input_listfile={psf_products}",
        f"--output_listfile={output_products}",
        "--reference_path=Data.ObservationSequence.ObservationId",
        "--input_path=Data.ObservationId",
        "--allow_duplicates",
    ]

    args = parser.parse_args(argstring)

    mainMethod(args)

    with open(workdir / output_products) as f:
        prods = [DmUtils.read_product_metadata(workdir / p) for p in json.load(f)]

    assert len(prods) == len(OBSERVATIONS)

    observations = [rgetattr(p, "Data.ObservationId") for p in prods]

    assert observations == OBSERVATIONS
