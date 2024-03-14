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
:file: tests/python/prunebatchexposures_test.py

:date: 04/03/2024
:author: Gordon Gibb
"""

import pytest
import json

import numpy as np

from ST_DM_DmUtils.DmUtils import read_product_metadata

from SHE_PPT.testing import (
    generate_mock_vis_images,
    generate_mock_mer_catalogues,
)

from SHE_PPT.she_io.hdf5 import convert_to_hdf5

from SHE_CTE_PipelineUtility.PruneBatchExposures import mainMethod, defineSpecificProgramOptions

NUM_EXPOSURES = 3
NUM_DETECTORS = NUM_EXPOSURES
NUM_OBJECTS_PER_DET = 10


@pytest.fixture
def workdir(tmp_path):
    """Creates a workdir (via pytest's tmpdir fixture) and a datadir"""

    # make the datadir
    datadir = tmp_path / "data"
    datadir.mkdir()

    return tmp_path


@pytest.fixture
def vis_and_mer_listfiles(workdir):
    """Creates the data products needed as an input"""

    # Create NUM_EXPOSURES exposures, with NUM_OBJECTS_PER_DET objects per detector.
    # The first exposure contains 1 detector, the second contains 2 etc...
    # This means that each subsequent exposure contains more detectors, and accordingly more objects.
    # The last exposure contains all the objects present in all previous exposures, so we can use the last
    # exposure's objects to create the catalogue of all objects.
    vis_products = []
    seg_products = []
    psf_products = []
    for i in range(NUM_EXPOSURES):
        exposure_product, sky_coords, pixel_coords, detectors, wcs_list = generate_mock_vis_images.create_exposure(
            workdir=workdir,
            n_objs_per_det=NUM_OBJECTS_PER_DET,
            n_detectors=i + 1,
            objsize=2,
        )
        vis_products.append(exposure_product)
        # NOTE: The executable does not try to read the segmap or psf products, so we can produce a dummy filename here
        seg_products.append(f"segmap_{i}.xml")
        psf_products.append(f"psf_{i}.xml")

    vis_listfile = "vis.json"
    with open(workdir / vis_listfile, "w") as f:
        json.dump(vis_products, f)

    seg_listfile = "seg.json"
    with open(workdir / seg_listfile, "w") as f:
        json.dump(seg_products, f)

    psf_listfile = "psf.json"
    with open(workdir / psf_listfile, "w") as f:
        json.dump(psf_products, f)

    # Create batches of objects such that each batch constitutes all objects in one detector.
    # This means that batch n will contain objects only in exposures n and onwards.
    mer_products = []
    for batch_coords in np.array_split(sky_coords, NUM_DETECTORS):
        mer_prod, _ = generate_mock_mer_catalogues.create_catalogue(obj_coords=batch_coords, workdir=workdir)
        mer_products.append(mer_prod)

    mer_listfile = "mer.json"
    with open(workdir / mer_listfile, "w") as f:
        json.dump(mer_products, f)

    return vis_listfile, seg_listfile, psf_listfile, mer_listfile


@pytest.fixture
def hdf5_listfile(workdir, vis_and_mer_listfiles):
    """Creates a listfile of hdf5 files"""

    vis_listfile, _, _, _ = vis_and_mer_listfiles

    datadir = workdir / "data"

    with open(workdir / vis_listfile) as f:
        vis_prods = json.load(f)

    hdf5_files = []
    for i, vis_prod in enumerate(vis_prods):
        dpd = read_product_metadata(workdir / vis_prod)
        det_file = datadir / dpd.Data.DataStorage.DataContainer.FileName
        bkg_file = datadir / dpd.Data.BackgroundStorage.DataContainer.FileName
        wgt_file = datadir / dpd.Data.WeightStorage.DataContainer.FileName
        # NOTE: seg_file is needed to create the HDF5, but is not needed for this test,
        # so just use the similarly structured BKG file as a placeholder
        seg_file = bkg_file

        output_filename = f"{i}.h5"
        hdf5_files.append(output_filename)

        convert_to_hdf5(det_file, bkg_file, wgt_file, seg_file, workdir / output_filename)

    hdf5_listfile = "hdf5_files.json"

    with open(workdir / hdf5_listfile, "w") as f:
        json.dump(hdf5_files, f)

    return hdf5_listfile


def test_executable_vis_only(workdir, vis_and_mer_listfiles):
    """Tests that the executable correctly generates listfiles for each batch (VIS inputs only)"""
    vis_listfile, seg_listfile, psf_listfile, mer_listfile = vis_and_mer_listfiles

    output_vis_listfile = "output_vis_listfile.json"
    output_seg_listfile = "output_seg_listfile.json"
    output_psf_listfile = "output_psf_listfile.json"

    argstring = [
        f"--workdir={workdir}",
        f"--vis_exposure_listfile={vis_listfile}",
        f"--segmap_exposure_listfile={seg_listfile}",
        f"--psf_exposure_listfile={psf_listfile}",
        f"--batch_mer_listfile={mer_listfile}",
        f"--pruned_batch_vis_listfile={output_vis_listfile}",
        f"--pruned_batch_psf_listfile={output_psf_listfile}",
        f"--pruned_batch_seg_listfile={output_seg_listfile}"
    ]

    parser = defineSpecificProgramOptions()

    args = parser.parse_args(argstring)

    mainMethod(args)

    with open(workdir / output_vis_listfile) as f:
        output_vis_batches = json.load(f)

    with open(workdir / output_psf_listfile) as f:
        output_psf_batches = json.load(f)

    with open(workdir / output_seg_listfile) as f:
        output_seg_batches = json.load(f)

    with open(workdir / vis_listfile) as f:
        vis_prods = json.load(f)

    with open(workdir / seg_listfile) as f:
        seg_prods = json.load(f)

    with open(workdir / psf_listfile) as f:
        psf_prods = json.load(f)

    # Should be one batch per exposure
    assert len(output_vis_batches) == NUM_DETECTORS
    assert len(output_seg_batches) == NUM_DETECTORS
    assert len(output_psf_batches) == NUM_DETECTORS

    for i, (vis_batch, seg_batch, psf_batch) in enumerate(
        zip(output_vis_batches, output_seg_batches, output_psf_batches)
    ):
        # First batch should contain all exposures, second should contain 2 onwards, third 3 onwards etc...
        assert vis_batch == vis_prods[i:]
        assert seg_batch == seg_prods[i:]
        assert psf_batch == psf_prods[i:]


def test_executable_vis_and_hdf5(workdir, vis_and_mer_listfiles, hdf5_listfile):
    """Tests that the executable correctly generates listfiles for each batch (VIS and HDF5 inputs)"""
    vis_listfile, seg_listfile, psf_listfile, mer_listfile = vis_and_mer_listfiles

    output_vis_listfile = "output_vis_listfile.json"
    output_seg_listfile = "output_seg_listfile.json"
    output_psf_listfile = "output_psf_listfile.json"
    output_hdf5_listfile = "output_hdf5_listfile.json"

    argstring = [
        f"--workdir={workdir}",
        f"--vis_exposure_listfile={vis_listfile}",
        f"--segmap_exposure_listfile={seg_listfile}",
        f"--psf_exposure_listfile={psf_listfile}",
        f"--hdf5_listfile={hdf5_listfile}",
        f"--batch_mer_listfile={mer_listfile}",
        f"--pruned_batch_vis_listfile={output_vis_listfile}",
        f"--pruned_batch_seg_listfile={output_seg_listfile}",
        f"--pruned_batch_psf_listfile={output_psf_listfile}",
        f"--pruned_batch_hdf5_listfile={output_hdf5_listfile}",
    ]

    parser = defineSpecificProgramOptions()

    args = parser.parse_args(argstring)

    mainMethod(args)

    with open(workdir / output_vis_listfile) as f:
        output_vis_batches = json.load(f)

    with open(workdir / output_seg_listfile) as f:
        output_seg_batches = json.load(f)

    with open(workdir / output_psf_listfile) as f:
        output_psf_batches = json.load(f)

    with open(workdir / output_hdf5_listfile) as f:
        output_hdf5_batches = json.load(f)

    with open(workdir / vis_listfile) as f:
        vis_prods = json.load(f)

    with open(workdir / seg_listfile) as f:
        seg_prods = json.load(f)

    with open(workdir / psf_listfile) as f:
        psf_prods = json.load(f)

    with open(workdir / hdf5_listfile) as f:
        hdf5_files = json.load(f)

    # Should be one batch per exposure
    assert len(output_vis_batches) == NUM_DETECTORS
    assert len(output_hdf5_batches) == NUM_DETECTORS
    assert len(output_seg_batches) == NUM_DETECTORS
    assert len(output_psf_batches) == NUM_DETECTORS

    for i, (vis_batch, seg_batch, psf_batch, hdf5_batch) in enumerate(
        zip(output_vis_batches, output_seg_batches, output_psf_batches, output_hdf5_batches)
    ):
        # First batch should contain all exposures, second should contain 2 onwards, third 3 onwards etc...
        assert vis_batch == vis_prods[i:]
        assert seg_batch == seg_prods[i:]
        assert hdf5_batch == hdf5_files[i:]
        assert psf_batch == psf_prods[i:]


def test_inconsistent_inputs(workdir, vis_and_mer_listfiles, hdf5_listfile):
    """Tests that the executable raises exceptions if inconsistent combinations of HDF5 arguments are provided"""
    vis_listfile, seg_listfile, psf_listfile, mer_listfile = vis_and_mer_listfiles

    output_vis_listfile = "output_vis_listfile.json"
    output_seg_listfile = "output_seg_listfile.json"
    output_hdf5_listfile = "output_hdf5_listfile.json"
    output_psf_listfile = "output_psf_listfile.json"

    parser = defineSpecificProgramOptions()

    # Test if only input HDF5 is provided
    argstring = [
        f"--workdir={workdir}",
        f"--vis_exposure_listfile={vis_listfile}",
        f"--segmap_exposure_listfile={seg_listfile}",
        f"--psf_exposure_listfile={psf_listfile}",
        f"--hdf5_listfile={hdf5_listfile}",
        f"--batch_mer_listfile={mer_listfile}",
        f"--pruned_batch_vis_listfile={output_vis_listfile}",
        f"--pruned_batch_seg_listfile={output_seg_listfile}",
        f"--pruned_batch_psf_listfile={output_psf_listfile}",
    ]

    args = parser.parse_args(argstring)

    with pytest.raises(ValueError, match="Inconsistent arguments"):
        mainMethod(args)

    # test if only output HDF5 is provided
    argstring = [
        f"--workdir={workdir}",
        f"--vis_exposure_listfile={vis_listfile}",
        f"--segmap_exposure_listfile={seg_listfile}",
        f"--psf_exposure_listfile={psf_listfile}",
        f"--batch_mer_listfile={mer_listfile}",
        f"--pruned_batch_vis_listfile={output_vis_listfile}",
        f"--pruned_batch_seg_listfile={output_seg_listfile}",
        f"--pruned_batch_hdf5_listfile={output_hdf5_listfile}",
        f"--pruned_batch_psf_listfile={output_psf_listfile}",
    ]
    args = parser.parse_args(argstring)

    with pytest.raises(ValueError, match="Inconsistent arguments"):
        mainMethod(args)
