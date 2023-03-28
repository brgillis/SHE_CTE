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
:file: python/SHE_CTE_PipelineUtility/ConvertToHDF5.py

:date: 02/11/2022
:author: Gordon Gibb

"""

import argparse
import os

import ElementsKernel.Logging as log

from ST_DM_DmUtils.DmUtils import read_product_metadata

from SHE_PPT.she_io.hdf5 import convert_to_hdf5


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--vis_frame_prod", type=str, help="The VIS product to convert")

    parser.add_argument("--remapped_seg_prod", type=str, help="The Remapped Segmentation Map product")

    parser.add_argument("--output_hdf5", type=str, help="The name of the output hdf5 file")

    parser.add_argument("--workdir", type=str, default=".", help="The working directory")

    parser.add_argument("--chunksize", type=int, default=100, help="The chunk size in pixels (0 means no chunking)")

    parser.add_argument("--logdir", type=str, default=".", help="The log directory")

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger("ConvertToHDF5")

    logger.info("#")
    logger.info("# Entering ConvertToHDF5 mainMethod()")
    logger.info("#")

    # read the input products
    qualified_vis_fname = os.path.join(args.workdir, args.vis_frame_prod)
    if not os.path.exists(qualified_vis_fname):
        raise FileNotFoundError("VIS frame product %s does not exist" % qualified_vis_fname)

    vis_dpd = read_product_metadata(qualified_vis_fname)

    qualified_seg_fname = os.path.join(args.workdir, args.remapped_seg_prod)
    if not os.path.exists(qualified_seg_fname):
        raise FileNotFoundError("Remapped segmentation map product %s does not exist" % qualified_seg_fname)

    seg_dpd = read_product_metadata(qualified_seg_fname)

    # get the filenames of each FITS file we need

    def qualified_fits_fname(filename):
        return os.path.join(args.workdir, "data", filename)

    det_file = qualified_fits_fname(vis_dpd.Data.DataStorage.DataContainer.FileName)
    wgt_file = qualified_fits_fname(vis_dpd.Data.WeightStorage.DataContainer.FileName)
    bkg_file = qualified_fits_fname(vis_dpd.Data.BackgroundStorage.DataContainer.FileName)
    seg_file = qualified_fits_fname(seg_dpd.Data.DataStorage.DataContainer.FileName)

    output_filename = os.path.join(args.workdir, args.output_hdf5)

    chunksize = args.chunksize
    if chunksize > 0:
        chunk = (chunksize, chunksize)
    else:
        chunk = None

    # create the HDF5 file

    convert_to_hdf5(det_file, bkg_file, wgt_file, seg_file, output_filename, chunk=chunk)

    logger.info("Input files successfully converted to HDF5")

    logger.info("#")
    logger.info("# Exiting ConvertToHDF5 mainMethod()")
    logger.info("#")

