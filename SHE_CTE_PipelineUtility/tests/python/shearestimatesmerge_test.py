""" @file objectidsplit_smoketest.py

    Created 17 Oct 2022

    Smoke test for SHE_CTE_ObjectIdSplit
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

# NOTE: To add to PPT eventually?
def create_mock_chains_product(workdir, ids=[]):
    """Creates a mock LensMC chains product from a list of object ids

    Inputs:
      - workdir: the workdir for the product
      - ids: a list of object IDs

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
    qualified_chains_filename = os.path.join(datadir, chains_fname)
    t_chains.write(qualified_chains_filename)

    dpd = create_dpd_she_lensmc_chains()
    dpd.set_data_filename(chains_fname)

    product_filename = get_allowed_filename("CHAINS", "PROD", version=CTE_VERSION, subdir="", extension=".xml")

    write_xml_product(dpd, product_filename, workdir=workdir)

    return product_filename


# NOTE: To add to PPT eventually?
def create_mock_measurements_product(workdir, ids=[], lensmc=False, ksb=False, momentsml=False, regauss=False):
    """Creates a mock shear measurements product from a list of object ids

    Inputs:
      - workdir: the workdir for the product
      - ids: a list of object IDs
      - lensmc: whether to create LensMC measurements table or not
      - ksb: whether to create KSB measurements table or not
      - momentsml: whether to create MomentsML measurements table or not
      - regauss: whether to create REGAUSS measurements table or not

    Returns:
      - product_filename: The name of the created product
    """

    datadir = os.path.join(workdir, "data")

    ids = np.asarray(ids, np.int64)
    g1 = np.random.random(len(ids))
    g2 = np.random.random(len(ids))

    if lensmc:
        t_lmc = init_table(lmc_tf, init_cols={lmc_tf.ID: ids, lmc_tf.g1: g1, lmc_tf.g2: g2})
        lmc_fname = get_allowed_filename("MEASUREMENTS", "LENSMC", version=CTE_VERSION, subdir="")
        qualified_lmc_filename = os.path.join(datadir, lmc_fname)
        t_lmc.write(qualified_lmc_filename)

    if ksb:
        t_ksb = init_table(ksb_tf, init_cols={ksb_tf.ID: ids, ksb_tf.g1: g1, ksb_tf.g2: g2})
        ksb_fname = get_allowed_filename("MEASUREMENTS", "KSB", version=CTE_VERSION, subdir="")
        qualified_ksb_filename = os.path.join(datadir, ksb_fname)
        t_ksb.write(qualified_ksb_filename)

    if regauss:
        t_regauss = init_table(regauss_tf, init_cols={regauss_tf.ID: ids, regauss_tf.g1: g1, regauss_tf.g2: g2})
        regauss_fname = get_allowed_filename("MEASUREMENTS", "REGAUSS", version=CTE_VERSION, subdir="")
        qualified_regauss_filename = os.path.join(datadir, regauss_fname)
        t_regauss.write(qualified_regauss_filename)

    if momentsml:
        t_momentsml = init_table(
            momentsml_tf, init_cols={momentsml_tf.ID: ids, momentsml_tf.g1: g1, momentsml_tf.g2: g2}
        )
        momentsml_fname = get_allowed_filename("MEASUREMENTS", "MOMENTSML", version=CTE_VERSION, subdir="")
        qualified_momentsml_filename = os.path.join(datadir, momentsml_fname)
        t_momentsml.write(qualified_momentsml_filename)

    dpd = create_she_measurements_product()

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

    def test_shearestimatesmerge(self, workdir):
        """Smoketest for the executable"""

        # Which methods to make products for
        # NOTE: Make some False so we can check that it correctly handles some methods missing
        lensmc = True
        ksb = True
        regauss = False
        momentsml = False

        # create input data

        ids1 = [i for i in range(10)]
        ids2 = [i for i in range(10, 20)]

        measurements_prod1 = create_mock_measurements_product(
            workdir=workdir, ids=ids1, lensmc=lensmc, ksb=ksb, regauss=regauss, momentsml=momentsml
        )
        measurements_prod2 = create_mock_measurements_product(
            workdir=workdir, ids=ids2, lensmc=lensmc, ksb=ksb, regauss=regauss, momentsml=momentsml
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

        # Check existence of outputs

        q_measurements_out = os.path.join(workdir, measurements_out)
        q_chains_out = os.path.join(workdir, chains_out)

        assert os.path.exists(q_measurements_out), "Output measurements product does not exist"
        assert os.path.exists(q_chains_out), "Output chains product does not exist"

        # check validity of measurements product
        dpd = read_xml_product(q_measurements_out)

        lmc_fits = dpd.get_LensMC_filename()
        ksb_fits = dpd.get_KSB_filename()
        regauss_fits = dpd.get_REGAUSS_filename()
        momentsml_fits = dpd.get_MomentsML_filename()

        methods = ["LensMC", "KSB", "REGAUSS", "MomentsML"]
        requested = [lensmc, ksb, regauss, momentsml]
        fits_files = [lmc_fits, ksb_fits, regauss_fits, momentsml_fits]

        for method, request, fits in zip(methods, requested, fits_files):
            # If this method was requested, make sure its file exists and all the objects are there
            if request:
                assert fits is not None, "%s file is not set in the product but it should be" % method
                t = Table.read(os.path.join(workdir, fits))
                assert set(t["OBJECT_ID"]) == set(ids1 + ids2), "%s object ids do not match those expected" % method
            else:
                assert fits is None, "%s file is set in the product but it should not be" % method

        # check validity of chains product
        dpd = read_xml_product(q_chains_out)
        fits = dpd.get_data_filename()

        t = Table.read(os.path.join(workdir, fits))
        assert set(t["OBJECT_ID"]) == set(ids1 + ids2), "%s object ids do not match those expected"
