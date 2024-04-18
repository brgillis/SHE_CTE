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
:file: tests/python/detector_batches_test.py

:date: 2024/04/05
:author: @rrollins

"""

from astropy.io import fits
from astropy.table import Table
import pytest

from SHE_PPT.file_io import read_listfile
from SHE_PPT.testing.generate_mock_vis_images import create_exposure
from SHE_PPT.testing.generate_mock_mer_catalogues import create_catalogue
from ST_DM_DmUtils.DmUtils import read_product_metadata

from SHE_CTE_PipelineUtility.DetectorBatches import mainMethod, defineSpecificProgramOptions


N_DET = 2
N_OBJ = 5
USE_QUAD = False


@pytest.fixture
def input_products(tmp_path):
    vis_filenme, sky_coords, _, _, _ = create_exposure(
        n_detectors=N_DET, workdir=tmp_path, n_objs_per_det=N_OBJ, use_quadrant=USE_QUAD
    )
    mer_filename, object_ids = create_catalogue(obj_coords=sky_coords, workdir=tmp_path)
    return vis_filenme, mer_filename, object_ids


def test_detector_batches(tmp_path, input_products):

    vis_filename, mer_filename, object_ids = input_products
    vis_batch_file = 'vis_batch.json'
    mer_batch_file = 'mer_batch.json'
    she_objid_file = 'she_objid.json'
    batch_size = 3
    max_batches = 3

    parser = defineSpecificProgramOptions()
    argstring = [
        f"--workdir={tmp_path}",
        f"--vis_calibrated_frame={vis_filename}",
        f"--mer_final_catalog={mer_filename}",
        f"--vis_detector_batches={vis_batch_file}",
        f"--mer_catalog_batches={mer_batch_file}",
        f"--she_objectid_lists={she_objid_file}",
        f"--batch_size={batch_size}",
        f"--max_batches={max_batches}",
    ]

    args = parser.parse_args(argstring)
    mainMethod(args)

    expected_ccd_ids = ['1-1', '1-1', '1-2']
    expected_batch_ids = [0, 1, 0]
    expected_batch_sizes = [3, 2, 3]
    expected_objids = object_ids[:sum(expected_batch_sizes)]
    mer_objids = list()
    she_objids = list()

    vis_batch_list = read_listfile(vis_batch_file, workdir=tmp_path)
    assert len(vis_batch_list) == max_batches
    for batch, ccd_id in zip(vis_batch_list, expected_ccd_ids):
        dpd = read_product_metadata(tmp_path / batch)
        with fits.open(tmp_path / 'data' / dpd.Data.DataStorage.DataContainer.FileName) as hdu_list:
            assert len(hdu_list) == 4
            assert hdu_list[1].name == f'{ccd_id}.SCI'
            assert hdu_list[2].name == f'{ccd_id}.RMS'
            assert hdu_list[3].name == f'{ccd_id}.FLG'
        with fits.open(tmp_path / 'data' / dpd.Data.BackgroundStorage.DataContainer.FileName) as hdu_list:
            assert len(hdu_list) == 2
            assert hdu_list[1].name == f'{ccd_id}'
        with fits.open(tmp_path / 'data' / dpd.Data.WeightStorage.DataContainer.FileName) as hdu_list:
            assert len(hdu_list) == 2
            assert hdu_list[1].name == f'{ccd_id}'

    mer_batch_list = read_listfile(mer_batch_file, workdir=tmp_path)
    assert len(mer_batch_list) == max_batches
    for batch, size in zip(mer_batch_list, expected_batch_sizes):
        dpd = read_product_metadata(tmp_path / batch)
        table = Table.read(tmp_path / 'data' / dpd.Data.DataStorage.DataContainer.FileName)
        assert len(table) == size
        mer_objids += list(table['OBJECT_ID'])
    assert mer_objids == expected_objids

    she_objid_list = read_listfile(she_objid_file, workdir=tmp_path)
    assert len(she_objid_list) == max_batches
    for batch, bid, size in zip(she_objid_list, expected_batch_ids, expected_batch_sizes):
        dpd = read_product_metadata(tmp_path / batch)
        assert (dpd.Data.BatchIndex) == bid
        assert len(dpd.Data.ObjectIdList) == size
        she_objids += dpd.Data.ObjectIdList
    assert she_objids == expected_objids
