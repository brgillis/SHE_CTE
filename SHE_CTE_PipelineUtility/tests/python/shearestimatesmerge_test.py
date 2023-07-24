""" @file shearestimatesmerge_test.py

    Created 17 Oct 2022

    Smoke test for SHE_CTE_ShearEstimatesMerge
"""

__updated__ = "2022-10-17"

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
import os
import json

import numpy as np

from astropy.table import Table

from SHE_PPT.file_io import get_allowed_filename, write_xml_product, read_xml_product
from SHE_PPT.table_utility import init_table

from SHE_PPT.table_formats.she_lensmc_measurements import tf as lmc_tf
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksb_tf
from SHE_PPT.table_formats.she_regauss_measurements import tf as regauss_tf
from SHE_PPT.table_formats.she_momentsml_measurements import tf as momentsml_tf

from SHE_PPT.table_formats.she_lensmc_chains import tf as chains_tf

from SHE_PPT.products.she_measurements import create_she_measurements_product
from SHE_PPT.products.she_lensmc_chains import create_dpd_she_lensmc_chains

from SHE_CTE import __version__ as CTE_VERSION

from SHE_CTE_PipelineUtility.ShearEstimatesMerge import defineSpecificProgramOptions, mainMethod

OBS_IDS = [12345, 6789]
POINTING_IDS = [1,2,3]
TILE_IDS = [4,5] 

# Which methods to make products for
# NOTE: Make some False so we can check that it correctly handles some methods missing
LENSMC = True
KSB = True
REGAUSS = False
MOMENTSML = False
REQUESTED_METHODS = [LENSMC, KSB, REGAUSS, MOMENTSML]
METHODS = ["LensMC", "KSB", "REGAUSS", "MomentsML"]

# NOTE: To add to PPT eventually?
def create_mock_chains_product(workdir, ids=[], include_data=True, obs_ids=OBS_IDS, pointing_ids=POINTING_IDS, tile_ids=TILE_IDS):
    """Creates a mock LensMC chains product from a list of object ids

    Inputs:
      - workdir: the workdir for the product
      - ids: a list of object IDs
      - include_data: bool - actually create a FITS file (T) or not (F)
      - obs_ids: The list of observation ids
      - pointing_ids: the list of pointing ids
      - tile_ids: the list of MER tile ids

    Returns:
      - product_filename: The name of the created product
    """

    datadir = os.path.join(workdir, "data")

    ids = np.asarray(ids, np.int64)

    t_chains = init_table(
        chains_tf,
        init_cols={
            lmc_tf.ID: ids,
        },
    )
    chains_fname = get_allowed_filename("CHAINS", "LENSMC", version=CTE_VERSION, subdir="")
    if include_data:
        qualified_chains_filename = os.path.join(datadir, chains_fname)
        t_chains.write(qualified_chains_filename)

    dpd = create_dpd_she_lensmc_chains()
    dpd.set_data_filename(chains_fname)
    dpd.Data.ObservationIdList = obs_ids
    dpd.Data.PointingIdList = pointing_ids
    dpd.Data.TileIndexList = tile_ids

    product_filename = get_allowed_filename("CHAINS", "PROD", version=CTE_VERSION, subdir="", extension=".xml")

    write_xml_product(dpd, product_filename, workdir=workdir)

    return product_filename


# NOTE: To add to PPT eventually?
def create_mock_measurements_product(workdir, ids=[], lensmc=False, ksb=False, momentsml=False, regauss=False, include_data=True, obs_ids=OBS_IDS, pointing_ids=POINTING_IDS, tile_ids=TILE_IDS):
    """Creates a mock shear measurements product from a list of object ids

    Inputs:
      - workdir: the workdir for the product
      - ids: a list of object IDs
      - lensmc: whether to create LensMC measurements table or not
      - ksb: whether to create KSB measurements table or not
      - momentsml: whether to create MomentsML measurements table or not
      - regauss: whether to create REGAUSS measurements table or not
      - include_data: bool - actually create a FITS file (T) or not (F)
      - obs_ids: the list of observation ids
      - pointing_ids: the list of pointing ids
      - tile_ids: the list of MER tile ids

    Returns:
      - product_filename: The name of the created product
    """

    datadir = os.path.join(workdir, "data")

    ids = np.asarray(ids, np.int64)
    e1 = np.random.random(len(ids))
    e2 = np.random.random(len(ids))

    if lensmc:
        t_lmc = init_table(lmc_tf, init_cols={lmc_tf.ID: ids, lmc_tf.e1: e1, lmc_tf.e2: e2})
        lmc_fname = get_allowed_filename("MEASUREMENTS", "LENSMC", version=CTE_VERSION, subdir="")
        qualified_lmc_filename = os.path.join(datadir, lmc_fname)
        if include_data:
            t_lmc.write(qualified_lmc_filename)

    if ksb:
        t_ksb = init_table(ksb_tf, init_cols={ksb_tf.ID: ids, ksb_tf.e1: e1, ksb_tf.e2: e2})
        ksb_fname = get_allowed_filename("MEASUREMENTS", "KSB", version=CTE_VERSION, subdir="")
        qualified_ksb_filename = os.path.join(datadir, ksb_fname)
        if include_data:
            t_ksb.write(qualified_ksb_filename)

    if regauss:
        t_regauss = init_table(regauss_tf, init_cols={regauss_tf.ID: ids, regauss_tf.e1: e1, regauss_tf.e2: e2})
        regauss_fname = get_allowed_filename("MEASUREMENTS", "REGAUSS", version=CTE_VERSION, subdir="")
        qualified_regauss_filename = os.path.join(datadir, regauss_fname)
        if include_data:
            t_regauss.write(qualified_regauss_filename)

    if momentsml:
        t_momentsml = init_table(
            momentsml_tf, init_cols={momentsml_tf.ID: ids, momentsml_tf.e1: e1, momentsml_tf.e2: e2}
        )
        momentsml_fname = get_allowed_filename("MEASUREMENTS", "MOMENTSML", version=CTE_VERSION, subdir="")
        qualified_momentsml_filename = os.path.join(datadir, momentsml_fname)
        if include_data:
            t_momentsml.write(qualified_momentsml_filename)

    dpd = create_she_measurements_product()
    dpd.Data.ObservationIdList = obs_ids
    dpd.Data.PointingIdList = pointing_ids
    dpd.Data.TileIndexList = tile_ids

    if lensmc:
        dpd.set_LensMC_filename(lmc_fname)

    if ksb:
        dpd.set_KSB_filename(ksb_fname)

    if regauss:
        dpd.set_REGAUSS_filename(regauss_fname)

    if momentsml:
        dpd.set_MomentsML_filename(momentsml_fname)

    product_filename = get_allowed_filename("MEASUREMENTS", "PROD", version=CTE_VERSION, subdir="", extension=".xml")

    write_xml_product(dpd, product_filename, workdir=workdir)

    return product_filename


def validate_shear_measurements_product(measurements_product, workdir, expected_obj_ids):
    q_measurements_prod = os.path.join(workdir, measurements_product)

    assert os.path.exists(q_measurements_prod), "Output measurements product does not exist"

    # check validity of measurements product
    dpd = read_xml_product(q_measurements_prod)

    assert dpd.Data.ObservationIdList == OBS_IDS, "Output product's ObservationID is an unexpected value"
    assert dpd.Data.TileIndexList == TILE_IDS, "Output product's TileList is an unexpected value"
    assert dpd.Data.PointingIdList == POINTING_IDS, "Output product's ObservationID is an unexpected value"

    lmc_fits = dpd.get_LensMC_filename()
    ksb_fits = dpd.get_KSB_filename()
    regauss_fits = dpd.get_REGAUSS_filename()
    momentsml_fits = dpd.get_MomentsML_filename()

    fits_files = [lmc_fits, ksb_fits, regauss_fits, momentsml_fits]

    for method, request, fits in zip(METHODS, REQUESTED_METHODS, fits_files):
        # If this method was requested, make sure its file exists and all the objects are there
        if request:
            assert fits is not None, "%s file is not set in the product but it should be" % method
            t = Table.read(os.path.join(workdir, fits))
            assert set(t["OBJECT_ID"]) == set(expected_obj_ids), "%s object ids do not match those expected" % method
            assert len(t) == len(expected_obj_ids), "Output table has an unexpected number of rows (expected %d, got %d)"%(len(t), len(ids1))
        else:
            assert fits is None, "%s file is set in the product but it should not be" % method


def validate_chains_product(chains_product, workdir, expected_obj_ids):
    q_chains_prod = os.path.join(workdir, chains_product)
    assert os.path.exists(q_chains_prod), "Output chains product does not exist"

    # check validity of chains product
    dpd = read_xml_product(q_chains_prod)

    assert dpd.Data.ObservationIdList == OBS_IDS, "Output product's ObservationID is an unexpected value"
    assert dpd.Data.TileIndexList == TILE_IDS, "Output product's TileList is an unexpected value"
    assert dpd.Data.PointingIdList == POINTING_IDS, "Output product's ObservationID is an unexpected value"
    
    fits = dpd.get_data_filename()

    t = Table.read(os.path.join(workdir, fits))
    assert set(t["OBJECT_ID"]) == set(expected_obj_ids), "%s object ids do not match those expected"



@pytest.fixture
def workdir(tmpdir):
    """Fixture that generates a workdir, and creates the 'data' subdirectory"""
    datadir = os.path.join(tmpdir, "data")
    os.mkdir(datadir)
    return tmpdir


class TestCase:
    """
    Tests for the SHE_CTE_ShearEstimatesMerge executable
    """

    def test_shearestimatesmerge_withchains(self, workdir):
        """Smoketest for the executable - with LensMC chains"""

        # create input data

        ids1 = [i for i in range(10)]
        ids2 = [i for i in range(10, 20)]

        measurements_prod1 = create_mock_measurements_product(
            workdir=workdir, ids=ids1, lensmc=LENSMC, ksb=KSB, regauss=REGAUSS, momentsml=MOMENTSML
        )
        measurements_prod2 = create_mock_measurements_product(
            workdir=workdir, ids=ids2, lensmc=LENSMC, ksb=KSB, regauss=REGAUSS, momentsml=MOMENTSML
        )

        measurements_listfile = "mes.json"

        with open(os.path.join(workdir, measurements_listfile), "w") as f:
            json.dump([measurements_prod1, measurements_prod2], f)

        chains_prod1 = create_mock_chains_product(workdir=workdir, ids=ids1)
        chains_prod2 = create_mock_chains_product(workdir=workdir, ids=ids2)

        chains_listfile = "chains.json"

        with open(os.path.join(workdir, chains_listfile), "w") as f:
            json.dump([chains_prod1, chains_prod2], f)

        # create args for executable

        chains_out = "combined_chains.xml"
        measurements_out = "combined_measurements.xml"

        commandlineargs = [
            "--workdir=%s" % workdir,
            "--shear_estimates_product_listfile=%s" % measurements_listfile,
            "--she_lensmc_chains_listfile=%s" % chains_listfile,
            "--merged_she_measurements=%s" % measurements_out,
            "--merged_she_lensmc_chains=%s" % chains_out,
        ]

        # Run the executable
        parser = defineSpecificProgramOptions()
        args = parser.parse_args(commandlineargs)
        mainMethod(args)

        validate_shear_measurements_product(measurements_out, workdir, ids1+ids2)
        validate_chains_product(chains_out, workdir, ids1+ids2)


    def test_shearestimatesmerge_nochains(self, workdir):
        """Smoketest for the executable - without LensMC chains"""

        # create input data

        ids1 = [i for i in range(10)]
        ids2 = [i for i in range(10, 20)]

        measurements_prod1 = create_mock_measurements_product(
            workdir=workdir, ids=ids1, lensmc=LENSMC, ksb=KSB, regauss=REGAUSS, momentsml=MOMENTSML
        )
        measurements_prod2 = create_mock_measurements_product(
            workdir=workdir, ids=ids2, lensmc=LENSMC, ksb=KSB, regauss=REGAUSS, momentsml=MOMENTSML
        )

        measurements_listfile = "mes.json"

        with open(os.path.join(workdir, measurements_listfile), "w") as f:
            json.dump([measurements_prod1, measurements_prod2], f)


        # create args for executable

        measurements_out = "combined_measurements.xml"

        parser = defineSpecificProgramOptions()

        # first test that the executable correctly throws an error if one of the chains input/output products is specified but not the other

        commandlineargs = [
            "--workdir=%s" % workdir,
            "--shear_estimates_product_listfile=%s" % measurements_listfile,
            "--she_lensmc_chains_listfile=%s" % "dummyfile",
            "--merged_she_measurements=%s" % measurements_out,
        ]
        
        args = parser.parse_args(commandlineargs)

        with pytest.raises(ValueError):
            mainMethod(args)


        commandlineargs = [
            "--workdir=%s" % workdir,
            "--shear_estimates_product_listfile=%s" % measurements_listfile,
            "--merged_she_measurements=%s" % measurements_out,
            "--merged_she_lensmc_chains=%s" % "dummyfile",
        ]

        args = parser.parse_args(commandlineargs)

        with pytest.raises(ValueError):
            mainMethod(args)
        
        # Now run the executable proper

        commandlineargs = [
            "--workdir=%s" % workdir,
            "--shear_estimates_product_listfile=%s" % measurements_listfile,
            "--merged_she_measurements=%s" % measurements_out,
        ]

        args = parser.parse_args(commandlineargs)
        mainMethod(args)

        # validate outputs
        validate_shear_measurements_product(measurements_out, workdir, ids1+ids2)


    def test_shearestimatesmerge_missing_fits(self, workdir):
        """Smoketest for the executable - with products with missing FITS tables"""

        # create input data

        ids1 = [i for i in range(10)]
        ids2 = [i for i in range(10, 20)]

        measurements_prod1 = create_mock_measurements_product(
            workdir=workdir, ids=ids1, lensmc=LENSMC, ksb=KSB, regauss=REGAUSS, momentsml=MOMENTSML
        )
        measurements_prod2 = create_mock_measurements_product(
            workdir=workdir, ids=ids2, lensmc=LENSMC, ksb=KSB, regauss=REGAUSS, momentsml=MOMENTSML, include_data=False
        )

        measurements_listfile = "mes.json"

        with open(os.path.join(workdir, measurements_listfile), "w") as f:
            json.dump([measurements_prod1, measurements_prod2], f)

        chains_prod1 = create_mock_chains_product(workdir=workdir, ids=ids1)
        chains_prod2 = create_mock_chains_product(workdir=workdir, ids=ids2, include_data=False)

        chains_listfile = "chains.json"

        with open(os.path.join(workdir, chains_listfile), "w") as f:
            json.dump([chains_prod1, chains_prod2], f)

        # create args for executable

        chains_out = "combined_chains.xml"
        measurements_out = "combined_measurements.xml"

        commandlineargs = [
            "--workdir=%s" % workdir,
            "--shear_estimates_product_listfile=%s" % measurements_listfile,
            "--she_lensmc_chains_listfile=%s" % chains_listfile,
            "--merged_she_measurements=%s" % measurements_out,
            "--merged_she_lensmc_chains=%s" % chains_out,
        ]

        # Run the executable
        parser = defineSpecificProgramOptions()
        args = parser.parse_args(commandlineargs)
        mainMethod(args)

        # validate outputs
        validate_shear_measurements_product(measurements_out, workdir, ids1)
        validate_chains_product(chains_out, workdir, ids1)




    def test_shearestimatesmerge_missing_products(self, workdir):
        """Smoketest for the executable - with missing measurement products"""

        # create input data

        ids1 = [i for i in range(10)]
        

        measurements_prod1 = create_mock_measurements_product(
            workdir=workdir, ids=ids1, lensmc=LENSMC, ksb=KSB, regauss=REGAUSS, momentsml=MOMENTSML
        )
        measurements_prod2 = "i_am_not_a_measurements_file.xml"

        measurements_listfile = "mes.json"

        with open(os.path.join(workdir, measurements_listfile), "w") as f:
            json.dump([measurements_prod1, measurements_prod2], f)

        chains_prod1 = create_mock_chains_product(workdir=workdir, ids=ids1)
        chains_prod2 = "i_am_not_a_chains_file.xml"

        chains_listfile = "chains.json"

        with open(os.path.join(workdir, chains_listfile), "w") as f:
            json.dump([chains_prod1, chains_prod2], f)

        # create args for executable

        chains_out = "combined_chains.xml"
        measurements_out = "combined_measurements.xml"

        commandlineargs = [
            "--workdir=%s" % workdir,
            "--shear_estimates_product_listfile=%s" % measurements_listfile,
            "--she_lensmc_chains_listfile=%s" % chains_listfile,
            "--merged_she_measurements=%s" % measurements_out,
            "--merged_she_lensmc_chains=%s" % chains_out,
        ]

        # Run the executable
        parser = defineSpecificProgramOptions()
        args = parser.parse_args(commandlineargs)
        mainMethod(args)


        # validate outputs
        validate_shear_measurements_product(measurements_out, workdir, ids1)
        validate_chains_product(chains_out, workdir, ids1)






