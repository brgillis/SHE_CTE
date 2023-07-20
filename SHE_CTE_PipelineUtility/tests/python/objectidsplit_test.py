""" @file objectidsplit_smoketest.py

    Created 21 March 2022

    Smoke test for SHE_CTE_ObjectIdSplit
"""

__updated__ = "2022-03-21"

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

import os
import uuid
import json
import shutil
import math

import pytest

from astropy.table import Table

from ST_DM_DmUtils.DmUtils import read_product_metadata

from SHE_PPT.testing import generate_mock_mer_catalogues, generate_mock_vis_images
from SHE_PPT.file_io import read_listfile, write_listfile
from SHE_PPT.products.mer_final_catalog import dpdMerFinalCatalog
from SHE_PPT.products.she_object_id_list import dpdSheObjectIdList
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf
from SHE_PPT.table_utility import is_in_format


from SHE_CTE_PipelineUtility.ObjectIdSplit import defineSpecificProgramOptions, mainMethod

base_workdir="/tmp/SHE_CTE_tests/"

num_objects = 45
batch_size = 10



@pytest.fixture
def input_files(tmpdir):
    """Produces the input data needed by ObjectIdSplit"""

    workdir = tmpdir

    #create the vis images and mer catalogue
    exposure_product, obj_coords, _, _, _ = generate_mock_vis_images.create_exposure(workdir=workdir, n_objs_per_det=num_objects, objsize=2)
    catalogue_product, _ = generate_mock_mer_catalogues.create_catalogue(obj_coords=obj_coords, workdir=workdir)
    
    #create the contents of the listfiles for the vis and mer products
    exposure_list = [exposure_product for i in range(4)]
    catalogue_list = [catalogue_product]

    data_images = "data_images.json"
    mer_final_catalog_tables = "mer_final_catalog_tables.json"
    mer_final_catalog_product = catalogue_product
    pipeline_config = "pipeline_config.txt"

    #write them to disk
    write_listfile(os.path.join(workdir,data_images), exposure_list)
    write_listfile(os.path.join(workdir,mer_final_catalog_tables), catalogue_list)
    
    #create the pipeline config, specifying the batch size
    with open(os.path.join(workdir,pipeline_config),"w") as f:
        f.write("SHE_CTE_ObjectIdSplit_batch_size=%d\n"%batch_size)

    return data_images, mer_final_catalog_tables, mer_final_catalog_product, pipeline_config


class TestCase:
    """
       Tests for the SHE_CTE_ObjectIdSplit executable

    """

    def _test_objectidsplit(self, argstring):
        """Runs ObjectIdSplit with a set of input arguments, and checks its outputs are consistent"""

        parser = defineSpecificProgramOptions()

        args = parser.parse_args(argstring)
        
        #run the executable
        mainMethod(args)

        #check that the correct number of batches have been produced (for the mer catalogues and object lists)
        n_batches_expected = math.ceil(num_objects/batch_size)
        
        workdir = args.workdir
        batch_mer_catalogs_listfile = args.batch_mer_catalogs
        object_ids_prod = args.object_ids

        catalog_list = read_listfile(os.path.join(workdir,batch_mer_catalogs_listfile))
        assert len(catalog_list) >= n_batches_expected, (
            "Unexpected number of mer_final_catalog batches created. Got %d expected >= %d"%(len(catalog_list), n_batches_expected))

        object_list = read_listfile(os.path.join(workdir,object_ids_prod))
        assert len(object_list) >= n_batches_expected, (
            "Unexpected number of object_id batches created. Got %d, expected >= %d"%(len(object_list), n_batches_expected))

        #check the output final catalogs are the correct type, and their tables are in the correct format
        for dpd_xml in catalog_list:
            dpd = read_product_metadata(os.path.join(workdir,dpd_xml))

            #ensure the read in data product is of the correct type
            assert type(dpd) == dpdMerFinalCatalog, "Expected output to be dpdMerFinalCatalog, but got %s"%(type(dpd).__name__)
            
            #now read in the table, and ensure it is the correct format
            table_filename = dpd.Data.DataStorage.DataContainer.FileName
            qualified_table_filename = os.path.join(workdir,"data",table_filename)
            table = Table.read(qualified_table_filename)
            assert is_in_format(table, mfc_tf), "MER final catalog table is not in the correct format"

        # check that the output object id lists are of the correct type
        output_objs = []
        for product in object_list:
            dpd = read_product_metadata(os.path.join(workdir,product))

            assert type(dpd) == dpdSheObjectIdList, "Expected output to be dpdSheObjectIdList, but got %s"%(type(dpd).__name__)
            
            #extract the object ids in this product for use below
            ids = dpd.get_id_list()
            output_objs += ids


        #check we have the correct number of objects
        assert len(output_objs) == num_objects, "Expected a total of %d objects, but got %d"%(num_objects, len(output_objs))
        
        #check all the objects are unique
        all_objs = set(output_objs)
        assert len(all_objs) == num_objects, "Not all object_ids are unique"

        

    def test_objectidsplit_listfile(self, tmpdir, input_files):
        """Test ObjectIdSplit with a MER final catalogue listfile as input"""

        workdir = tmpdir

        data_images, mer_final_catalog_tables, _, pipeline_config = input_files
        
        argstring = [
            f"--workdir={workdir}",
            f"--data_images={data_images}",
            f"--mer_final_catalog_tables={mer_final_catalog_tables}",
            "--object_ids=object_ids0.json",
            "--batch_mer_catalogs=batch_mer_catalogs0.json",
            f"--pipeline_config={pipeline_config}",
        ]

        self._test_objectidsplit(argstring)

    def test_objectidsplit_product(self, tmpdir, input_files):
        """Test ObjectIdSplit with a MER final catalogue product as input"""
        
        workdir = tmpdir

        data_images, _, mer_final_catalog_product, pipeline_config = input_files
        
        argstring = [
            f"--workdir={workdir}",
            f"--data_images={data_images}",
            f"--mer_final_catalog_tables={mer_final_catalog_product}",
            "--object_ids=object_ids1.json",
            "--batch_mer_catalogs=batch_mer_catalogs1.json",
            f"--pipeline_config={pipeline_config}",
        ]

        self._test_objectidsplit(argstring)


    def test_objectidsplit_unknown(self, tmpdir, input_files):
        """Test ObjectIdSplit with a MER final catalogue file with an unexpected file extension"""
        
        workdir = tmpdir

        data_images, _, mer_final_catalog_product, pipeline_config = input_files
        
        argstring = [
            f"--workdir={workdir}",
            f"--data_images={data_images}",
            "--mer_final_catalog_tables=some_file.xyz",
            "--object_ids=object_ids1.json",
            "--batch_mer_catalogs=batch_mer_catalogs1.json",
            f"--pipeline_config={pipeline_config}",
        ]

        with pytest.raises(ValueError):
            self._test_objectidsplit(argstring)
        
