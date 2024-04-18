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
:file: python/SHE_CTE_PipelineUtility/SHE_CTE_PruneBatchExposures.py

:date: 02/03/24
:author: Gordon Gibb

This executable produces for each batch of objects, a listfile of _only_ the VIS exposures that cover each batch (and
optionally listfile of the HDF5 versions of the VIS images that cover each batch)

The output listfile is of the form:
[
    [exposure1.xml, exposure2.xml, ...],
    [exposure1.xml, exposure3.xml, ...],
    ...
]
It contains num_batches lists, each one being the list of exposures that cover that batch.

The pipeline runner can use this output listfile to generate a listfile of exposures for each batch in a parallel split
"""

import argparse
import json
import pathlib
import warnings

from dataclasses import dataclass

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.units import degree, UnitsWarning
from astropy.io.fits.verify import VerifyWarning

import ElementsKernel.Logging as log

from ST_DM_DmUtils import DmUtils
from ST_DM_DPTools import JsonTools

from SHE_PPT.argument_parser import dir_path
from SHE_PPT.she_io.vis_exposures import read_vis_data

logger = log.getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    # pipeline args

    parser.add_argument("--workdir", type=dir_path, default=".", help="Workdir")

    parser.add_argument("--logdir", type=dir_path, help="logdir")

    # input args

    parser.add_argument("--vis_exposure_listfile", type=str, required=True, help="Listfile pointing to VIS exposures")

    parser.add_argument(
        "--segmap_exposure_listfile", type=str, required=True,
        help="Listfile pointing to reprojected segmaps"
    )

    parser.add_argument(
        "--psf_exposure_listfile", type=str, required=True,
        help="Listfile of PSF field parameters for each exposure"
    )

    parser.add_argument("--batch_mer_listfile", type=str, required=True, help="Listfile of batch MER catalogues")

    parser.add_argument("--hdf5_listfile", type=str, help="Listfile of HDF5 files (optional)")

    # output args

    parser.add_argument(
        "--pruned_batch_vis_listfile", type=str, required=True,
        help="Listfile of listfiles pointing to vis exposures in each batch"
    )

    parser.add_argument(
        "--pruned_batch_seg_listfile", type=str, required=True,
        help="Listfile of listfiles pointing to the segmaps in each batch"
    )

    parser.add_argument(
        "--pruned_batch_psf_listfile", type=str, required=True,
        help="Listfile of listfiles pointing to the PSF field params in each batch"
    )

    parser.add_argument(
        "--pruned_batch_hdf5_listfile",
        type=str,
        help="Listfile of listfiles pointing to hdf5 VIS exposures in each batch (optional)",
    )

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger.info("#")
    logger.info("# Entering SHE_CTE_PruneBatchExposures mainMethod()")
    logger.info("#")

    workdir = args.workdir

    # Sanity check - if we have an HDF5 input, we should have an HDF5 output!
    if bool(args.hdf5_listfile) != bool(args.pruned_batch_hdf5_listfile):
        raise ValueError(
            f"Inconsistent arguments. hdf5_listfile={args.hdf5_listfile}, "
            f"pruned_batch_hdf5_listfile={args.pruned_batch_hdf5_listfile}. "
            "These must be either both set or both unset",
            args.hdf5_listfile,
            args.pruned_batch_hdf5_listfile,
        )

    vis_files = JsonTools.read_json(workdir / args.vis_exposure_listfile)
    seg_files = JsonTools.read_json(workdir / args.segmap_exposure_listfile)
    psf_files = JsonTools.read_json(workdir / args.psf_exposure_listfile)

    # To avoid code branches further down, if no HDF5 files are passed in, we just make a list of Nones
    # the same length as the VIS filelist
    if args.hdf5_listfile:
        method = "hdf5"
        hdf5_files = JsonTools.read_json(workdir / args.hdf5_listfile)
    else:
        method = "astropy"
        hdf5_files = [None for _ in vis_files]

    vis_exposures = read_vis_data(
        vis_files, workdir=workdir, hdf5_files=hdf5_files, method=method
    )

    exposure_wcs_lists = [exp.get_wcs_list() for exp in vis_exposures]

    exposure_metadata_list = [
        ExposureMetadata(w, v, s, p, h) for w, v, s, p, h in zip(
            exposure_wcs_lists, vis_files, seg_files, psf_files, hdf5_files
        )
    ]

    batch_mer_filenames = JsonTools.read_json(workdir / args.batch_mer_listfile)

    vis_batches = []
    seg_batches = []
    psf_batches = []
    hdf5_batches = []

    for i, batch_mer_filename in enumerate(batch_mer_filenames):
        mer_dpd = DmUtils.read_product_metadata(workdir / batch_mer_filename)
        mer_fits = mer_dpd.Data.DataStorage.DataContainer.FileName

        # Read batch catalogue in, ignoring warnings from MER catalogues not being in the correct format
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=VerifyWarning)
            warnings.simplefilter("ignore", category=UnitsWarning)
            mer_cat = Table.read(workdir / "data" / mer_fits, memmap=False)

        skycoords = SkyCoord(ra=mer_cat["RIGHT_ASCENSION"], dec=mer_cat["DECLINATION"], unit=degree)

        # get the lists of VIS/HDF5 exposure files that contain the objects for this batch
        vis_batch, seg_batch, psf_batch, hdf5_batch = _get_exposure_lists_with_objects(
            skycoords, exposure_metadata_list
        )

        logger.info("Batch %d: Objects are present in %d exposures", i, len(vis_batch))

        vis_batches.append(vis_batch)
        seg_batches.append(seg_batch)
        psf_batches.append(psf_batch)
        hdf5_batches.append(hdf5_batch)

    logger.info("Creating output listfile %s", workdir / args.pruned_batch_vis_listfile)
    _write_listfile(vis_batches, workdir / args.pruned_batch_vis_listfile)

    logger.info("Creating output listfile %s", workdir / args.pruned_batch_seg_listfile)
    _write_listfile(seg_batches, workdir / args.pruned_batch_seg_listfile)

    logger.info("Creating output listfile %s", workdir / args.pruned_batch_psf_listfile)
    _write_listfile(psf_batches, workdir / args.pruned_batch_psf_listfile)

    if args.pruned_batch_hdf5_listfile:
        logger.info("Creating output listfile %s", workdir / args.pruned_batch_hdf5_listfile)
        _write_listfile(hdf5_batches, workdir / args.pruned_batch_hdf5_listfile)

    logger.info("#")
    logger.info("# Exiting SHE_CTE_PruneBatchExposures mainMethod()")
    logger.info("#")


@dataclass
class ExposureMetadata:
    wcs_list: list["astropy.wcs.WCS"]  # NOQA F821
    vis_filename: str
    seg_filename: str
    psf_filename: str
    hdf5_filename: str = None


def _write_listfile(data, filename):
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


def _exposure_contains_objects(skycoords, exposure_meta):
    for wcs in exposure_meta.wcs_list:
        if wcs.footprint_contains(skycoords).any():
            return True
    return False


def _get_exposure_lists_with_objects(skycoords, exposure_metadata_list):
    """Determines which exposures a set of objects are in

    Inputs:
     - skycoords: a set of object SkyCoords
     - exposure_metadata_list: a list of ExposureMetadata objects
    Returns:
     - batch_vis_list: list of VIS products that contain the objects
     - batch_seg_list: list of Segmap products that contain the objects
     - batch_vis_list: list of PSF field params products that contain the objects
     - batch_hdf5_list: list of HDF5 files that contain the objects
    """

    batch_vis_list = []
    batch_seg_list = []
    batch_psf_list = []
    batch_hdf5_list = []

    for exposure_meta in exposure_metadata_list:

        if _exposure_contains_objects(skycoords, exposure_meta):
            batch_vis_list.append(exposure_meta.vis_filename)
            batch_seg_list.append(exposure_meta.seg_filename)
            batch_psf_list.append(exposure_meta.psf_filename)
            batch_hdf5_list.append(exposure_meta.hdf5_filename)

    return batch_vis_list, batch_seg_list, batch_psf_list, batch_hdf5_list
