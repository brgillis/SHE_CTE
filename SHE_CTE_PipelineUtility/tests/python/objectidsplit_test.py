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

from astropy.table import Table

from SHE_PPT.testing import generate_mock_mer_catalogues, generate_mock_vis_images
from SHE_PPT.file_io import read_listfile, read_xml_product, write_listfile
from SHE_PPT.products.mer_final_catalog import dpdMerFinalCatalog
from SHE_PPT.products.she_object_id_list import dpdSheObjectIdList
from SHE_PPT.table_formats.mer_final_catalog import tf as mfc_tf
from SHE_PPT.table_utility import is_in_format


from SHE_CTE_PipelineUtility.ObjectIdSplit import defineSpecificProgramOptions, mainMethod

base_workdir="/tmp/SHE_CTE_tests/"

num_objects = 45
batch_size = 10




class TestCase:
    """
       Tests for the SHE_CTE_ObjectIdSplit executable

    """

    @classmethod
    def setup_class(cls):
        
        #get a UUID for the workdir for this testcase
        workdir = str(uuid.uuid4())
        
        #create this workdir
        cls.workdir = os.path.join(base_workdir,workdir)
        os.makedirs(cls.workdir)

        #create the vis images and mer catalogue
        exposure_product, obj_coords = generate_mock_vis_images.create_exposure(workdir=cls.workdir, n_objs_per_det=num_objects)
        catalogue_product = generate_mock_mer_catalogues.create_catalogue(obj_coords=obj_coords, workdir=cls.workdir)
        
        #create the contents of the listfiles for the vis and mer products
        exposure_list = [exposure_product for i in range(4)]
        catalogue_list = [catalogue_product]

        #write them to disk
        cls.data_images = os.path.join(cls.workdir,"data_images.json")
        write_listfile(cls.data_images, exposure_list)
       
        cls.mer_final_catalog_tables = os.path.join(cls.workdir,"mer_final_catalog_tables.json")
        write_listfile(cls.mer_final_catalog_tables, catalogue_list)
        
        #create the pipeline config, specifying the batch size
        cls.pipeline_config = os.path.join(cls.workdir,"pipeline_config.txt")
        with open(cls.pipeline_config,"w") as f:
            f.write("SHE_CTE_ObjectIdSplit_batch_size=%d\n"%batch_size)


    @classmethod
    def teardown_class(cls):
        
        #remove the workdir
        shutil.rmtree(cls.workdir)

        return

    def test_objectidsplit(self):
        
        #define program arguments
        parser = defineSpecificProgramOptions()
        args = parser.parse_args([])
        
        args.workdir = self.workdir
        args.data_images=self.data_images
        args.mer_final_catalog_tables = self.mer_final_catalog_tables
        args.object_ids = "object_ids.json"
        args.batch_mer_catalogs = "batch_mer_catalogs.json"
        args.pipeline_config=self.pipeline_config
        
        #run the executable
        mainMethod(args)

        #check that the correct number of batches have been produced (for the mer catalogues and object lists)
        n_batches_expected = math.ceil(num_objects/batch_size)

        catalog_list = read_listfile(os.path.join(self.workdir,"batch_mer_catalogs.json"))
        assert len(catalog_list) == n_batches_expected, (
               "Unexpected number of mer_final_catalog batches created. Got %d expected %d"%(len(catalog_list), n_batches_expected))

        object_list = read_listfile(os.path.join(self.workdir,"object_ids.json"))
        assert len(object_list) == n_batches_expected, (
               "Unexpected number of object_id batches created. Got %d, expected %d"%(len(object_list), n_batches_expected))

        #check the output final catalogs are the correct type, and their tables are in the correct format
        for product_listfile in catalog_list:
            dpds = read_listfile(os.path.join(self.workdir,product_listfile))
            
            #make sure there is only one product in this listfile
            assert len(dpds) == 1, "Expected one data product per MER catalogue listfile, but got %d"%len(dpds)
            
            dpd_xml = dpds[0]

            dpd = read_xml_product(dpd_xml,workdir=self.workdir)

            #ensure the read in data product is of the correct type
            assert type(dpd) == dpdMerFinalCatalog, "Expected output to be dpdMerFinalCatalog, but got %s"%(type(dpd).__name__)
            
            #now read in the table, and ensure it is the correct format
            table_filename = dpd.get_filename()
            qualified_table_filename = os.path.join(self.workdir,table_filename)
            table = Table.read(qualified_table_filename)
            assert is_in_format(table, mfc_tf), "MER final catalog table is not in the correct format"

        # check that the output object id lists are of the correct type
        output_objs = []
        for product in object_list:
            dpd = read_xml_product(product,workdir=self.workdir)

            assert type(dpd) == dpdSheObjectIdList, "Expected output to be dpdSheObjectIdList, but got %s"%(type(dpd).__name__)
            
            #extract the object ids in this product for use below
            ids = dpd.get_id_list()
            output_objs += ids


        #check we have the correct number of objects
        assert len(output_objs) == num_objects, "Expected a total of %d objects, but got %d"%(num_objects, len(output_objs))
        
        #check all the objects are unique
        all_objs = set(output_objs)
        assert len(all_objs) == num_objects, "Not all object_ids are unique"


