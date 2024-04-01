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
:file: tests/python/utils_test.py

:date: 2024/01/23
:author: @rrollins

"""

from astropy.io import fits
from astropy.table import Table
from hypothesis import given
from hypothesis.strategies import integers, iterables
from itertools import islice, tee
import numpy as np
from pathlib import Path
import pytest
import sys

from ST_DataModelBindings.dpd.she.raw.intermediategeneral_stub import DpdSheIntermediateGeneral
from ST_DataModelBindings.pro.she_stub import sheIntermediateGeneral
from ST_DM_DmUtils.DmUtils import get_product_name, read_product_metadata
from ST_DM_HeaderProvider.GenericHeaderProvider import GenericHeaderGenerator

from SHE_CTE import SHE_CTE_RELEASE_STRING
from SHE_CTE_PipelineUtility import utils


class NotIterable:
    pass


@given(iterables(integers()), integers(min_value=1, max_value=sys.maxsize))
def test_batched_integers(iterable, batch_size):
    i1, i2 = tee(iterable)
    for batch in utils.batched(i1, batch_size):
        assert batch == tuple(islice(i2, batch_size))


def test_ceiling_division():

    # Integer Cases
    assert utils.ceiling_division(-1, 2) == 0
    assert utils.ceiling_division(0, 2) == 0
    assert utils.ceiling_division(1, 2) == 1
    assert utils.ceiling_division(2, 2) == 1
    assert utils.ceiling_division(3, 2) == 2

    # Float Cases
    assert utils.ceiling_division(-1.0, 1.5) == 0
    assert utils.ceiling_division(0.0, 1.5) == 0
    assert utils.ceiling_division(1.0, 1.5) == 1
    assert utils.ceiling_division(2.0, 1.5) == 2
    assert utils.ceiling_division(3.0, 1.5) == 2


@given(integers(min_value=1))
def test_batched_not_iterable(batch_size):
    with pytest.raises(TypeError, match='object is not iterable'):
        batched = utils.batched(NotIterable(), batch_size)
        next(batched)


@given(integers(max_value=0))
def test_batched_invalid_batch_size(batch_size):
    with pytest.raises(ValueError, match='batch_size must be at least one'):
        batched = utils.batched(None, batch_size)
        next(batched)


@pytest.fixture
def she_intermediate_general_data_product():
    dpd = DpdSheIntermediateGeneral()
    product_name = get_product_name(dpd)
    hg = GenericHeaderGenerator(product_type=product_name)
    hg.set_tag_value('SoftwareName', 'SHE_CTE')
    hg.set_tag_value('SoftwareRelease', SHE_CTE_RELEASE_STRING)
    hg.set_tag_value('Curator', 'SHE')
    dpd.Header = hg.generate()
    dpd.Data = sheIntermediateGeneral()
    dpd.Data.DataStorage = []
    return dpd


def test_write_data_product(she_intermediate_general_data_product, tmp_path):

    # Call write_data_product for simple data product
    type_name = 'TN'
    instance_id = 'IID'
    filename = utils.write_data_product(type_name, instance_id, she_intermediate_general_data_product, tmp_path)

    # Check filename is valid
    assert (tmp_path / filename).is_file()
    assert type_name in filename
    assert instance_id in filename
    assert SHE_CTE_RELEASE_STRING in filename
    assert 'SHE' in filename
    assert Path(filename).suffix == '.xml'

    # Check file contains a valid data product with expected contents
    dpd = read_product_metadata(tmp_path / filename)
    assert dpd.validateBinding()
    assert dpd.Header.SoftwareName == 'SHE_CTE'
    assert dpd.Header.SoftwareRelease == SHE_CTE_RELEASE_STRING
    assert dpd.Header.Curator == 'SHE'


def test_write_fits_hdus(tmp_path):

    # Call write_data_product for simple data product
    type_name = 'TN'
    instance_id = 'IID'
    hdu1 = fits.ImageHDU()
    hdu2 = fits.TableHDU()
    hdu1.header['TEST1'] = 'TEST1'
    hdu2.header['TEST2'] = 'TEST2'
    filename = utils.write_fits_hdus(type_name, instance_id, (hdu1, hdu2), tmp_path)

    # Check filename is valid
    assert (tmp_path / filename).is_file()
    assert type_name in filename
    assert instance_id in filename
    assert SHE_CTE_RELEASE_STRING in filename
    assert 'SHE' in filename
    assert Path(filename).suffix == '.fits'

    # Check file is a valid FITS file with expected contents
    with fits.open(tmp_path / filename) as hdulist:
        assert isinstance(hdulist[0], fits.PrimaryHDU)
        assert isinstance(hdulist[1], fits.ImageHDU)
        assert isinstance(hdulist[2], fits.TableHDU)
        assert hdulist[1].header['TEST1'] == 'TEST1'
        assert hdulist[2].header['TEST2'] == 'TEST2'


def test_write_fits_table(tmp_path):

    # Call write_fits_table for random astropy table
    type_name = 'TN'
    instance_id = 'IID'
    data = np.random.uniform(size=(2, 3))
    input_table = Table(data=data)
    filename = utils.write_fits_table(type_name, instance_id, input_table, tmp_path)

    # Check filename is valid
    assert (tmp_path / filename).is_file()
    assert type_name in filename
    assert instance_id in filename
    assert SHE_CTE_RELEASE_STRING in filename
    assert 'SHE' in filename
    assert Path(filename).suffix == '.fits'

    # Check file is a valid FITS file with expected contents
    fits_table = Table.read(tmp_path / filename)
    np.testing.assert_array_equal(input_table, fits_table)
