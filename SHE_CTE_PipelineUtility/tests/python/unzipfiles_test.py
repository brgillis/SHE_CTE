#
# Copyright (C) 2012-2023 Euclid Science Ground Segment
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
:file: tests/python/unzipfiles_test.py

:date: 31/08/2023
:author: Gordon Gibb
"""

import pytest
import subprocess
import pathlib

from ST_DM_DmUtils import DmUtils

from ST_DM_HeaderProvider.GenericHeaderProvider import GenericHeaderGenerator
from ST_DM_FilenameProvider.FilenameProvider import FileNameProvider
from ST_DataModelBindings.dpd.she.raw.intermediategeneral_stub import DpdSheIntermediateGeneral
from ST_DataModelBindings.pro.she_stub import sheIntermediateGeneral, sheGenericFile

from SHE_CTE import SHE_CTE_RELEASE_STRING
from SHE_CTE_PipelineUtility.UnzipFiles import convert_product, gunzip_file


@pytest.fixture
def workdir(tmp_path):
    workdir = tmp_path
    datadir = workdir / "data"
    datadir.mkdir()

    return workdir


@pytest.fixture
def compressed_file(workdir):
    """Create a file and its gzipped version"""

    unzipped_filename = "test_file.txt"
    unzipped_path = workdir / "data" / unzipped_filename
    with open(unzipped_path, "w") as f:
        f.write("Lorem ipsum dolor sit amet...")

    zipped_filename = unzipped_filename + ".gz"
    zipped_path = workdir / "data" / zipped_filename

    # Run the gzip command. We do not want to overwrite the input file.
    # Sadly the "-k" option is not available on some machines, so instead we redirect
    # the compressed output to stdout, capture that and write it to the output file
    result = subprocess.run(("gzip", "-c", f"{unzipped_path}"), capture_output=True)

    if result.returncode != 0:
        raise Exception(f"Command failed: {result.stderr}")
    with open(zipped_path, "wb") as f:
        f.write(result.stdout)

    return unzipped_filename, zipped_filename


def test_unzipping(workdir, compressed_file):
    """Check that the unzip function works"""

    datadir = workdir / "data"

    unzipped_file, gzipped_file = compressed_file

    reference_file = datadir / unzipped_file
    input_file = datadir / gzipped_file
    output_file = datadir / "uncompressed_file"

    gunzip_file(input_file, output_file)

    with open(reference_file) as f:
        reference_contents = f.read()

    with open(output_file) as f:
        output_contents = f.read()

    assert reference_contents == output_contents


def test_convert_product(workdir, compressed_file):
    """Tests that we can gunzip all gzipped files in a data product"""

    unzipped_file, gzipped_file = compressed_file

    # create input data product

    # filenames to store in the product
    data_filenames = [f"{gzipped_file}", "somefile.txt", "anotherfile.txt.gz"]

    # filenames we expect in the output
    expected_filenames = [f"{unzipped_file}", "somefile.txt", "anotherfile.txt"]

    # NOTE: convert_product should work with any data product, so this test is not limited to
    # dpdSheIntermediateGeneral, however this is a simple product (that we control) which can
    # contain an arbitrary number of data files.
    dpd = DpdSheIntermediateGeneral()

    product_name = DmUtils.get_product_name(dpd)

    hg = GenericHeaderGenerator(product_type=product_name)
    hg.set_tag_value("SoftwareName", "SHE_PSFToolkitFitting")
    hg.set_tag_value("SoftwareRelease", SHE_CTE_RELEASE_STRING)
    hg.set_tag_value("Curator", "SHE")

    dpd.Header = hg.generate()

    dpd.Data = sheIntermediateGeneral()

    dpd.Data.DataStorage = []
    for data_filename in data_filenames:
        shefile = sheGenericFile.Factory()
        shefile.DataContainer = DmUtils.create_data_container(data_filename, "PROPOSED")
        dpd.Data.DataStorage.append(shefile)

    input_product = "input.xml"
    DmUtils.save_product_metadata(dpd, workdir / input_product)

    output_product = "output.xml"

    # Run the function and check the output product has the expected filenames
    # NOTE: No need to check the decompression as this was done in test_unzipping

    convert_product(workdir, input_product, output_product)

    output_dpd = DmUtils.read_product_metadata(workdir / output_product)

    output_filenames = [ds.DataContainer.FileName for ds in output_dpd.Data.DataStorage]

    for output_fname, expected_fname in zip(output_filenames, expected_filenames):
        # NOTE: as we add a UUID to the stem of the unzipped file's name, it will not match exactly
        output_stem = pathlib.Path(output_fname).stem
        expected_stem = pathlib.Path(expected_fname).stem
        print(expected_fname, output_fname)
        assert expected_stem in output_stem

