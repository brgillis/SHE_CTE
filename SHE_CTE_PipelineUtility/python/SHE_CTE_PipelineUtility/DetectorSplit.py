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
:file: python/SHE_CTE_PipelineUtility/DetectorSplit.py

:date: 2024/01/22
:author: @rrollins

"""

import json
from argparse import ArgumentParser, ArgumentTypeError
from copy import deepcopy
from itertools import islice
from pathlib import Path

from astropy.io import fits
from ElementsKernel.Logging import getLogger
from SHE_PPT.file_io import read_xml_product, save_product_metadata
from ST_DM_FilenameProvider.FilenameProvider import FileNameProvider

from SHE_CTE import SHE_CTE_RELEASE_STRING
from SHE_CTE_PipelineUtility.utils import batched


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    def dir_path(path):
        if not Path(path).is_dir():
            raise ArgumentTypeError(f"{path} is not a valid path")
        return path

    parser = ArgumentParser()
    parser.add_argument('--workdir', type=dir_path, help='Name of the working directory')
    parser.add_argument('--logdir', type=str, help='Log directory')
    parser.add_argument('--vis_calibrated_frame', type=str, help='VIS calibrated frames data product.')
    parser.add_argument('--vis_detector_frames', type=str, help='Listfile of VIS calibrated frames data products.')

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger('SHE_CTE_DetectorSplit')

    logger.info('#')
    logger.info('# Entering SHE_CTE_DetectorSplit mainMethod()')
    logger.info('#')

    # Paths
    workdir = Path(args.workdir)
    datadir = workdir / 'data'

    # Read VIS product
    vis_product = read_xml_product(workdir / args.vis_calibrated_frame)
    observation_id = vis_product.Data.ObservationSequence.ObservationId
    dither_id = vis_product.Data.ObservationSequence.DitherObservation
    exposure_id = vis_product.Data.ObservationSequence.Exposure
    vis_det_list = vis_product.Data.DetectorList.Detector
    vis_det_file = datadir / vis_product.Data.DataStorage.DataContainer.FileName
    vis_bkg_file = datadir / vis_product.Data.BackgroundStorage.DataContainer.FileName
    vis_wgt_file = datadir / vis_product.Data.WeightStorage.DataContainer.FileName

    filename_provider = FileNameProvider()
    vis_products = list()

    with fits.open(vis_det_file) as det_hdus, fits.open(vis_bkg_file) as bkg_hdus, fits.open(vis_wgt_file) as wgt_hdus:
        batched_det_hdus = batched(islice(det_hdus, 1, None), 3)
        ibkg_hdus = islice(bkg_hdus, 1, None)
        iwgt_hdus = islice(wgt_hdus, 1, None)
        for detector, det, bkg, wgt in zip(vis_det_list, batched_det_hdus, ibkg_hdus, iwgt_hdus):

            detector_id = detector.DetectorId
            instance_id = f'{observation_id:06}_{dither_id:02}_{exposure_id}_{detector_id}'

            # VIS DET
            vis_detector_det_filename = filename_provider.get_allowed_filename(
                type_name="DET",
                instance_id=instance_id,
                extension=".fits",
                release=SHE_CTE_RELEASE_STRING,
                processing_function="SHE",
            )
            fits.HDUList([fits.PrimaryHDU(), *det]).writeto(datadir / vis_detector_det_filename)
            logger.info('Wrote VIS detector data file [filename=%s]', vis_detector_det_filename)

            # VIS BKG
            vis_detector_bkg_filename = filename_provider.get_allowed_filename(
                type_name="BKG",
                instance_id=instance_id,
                extension=".fits",
                release=SHE_CTE_RELEASE_STRING,
                processing_function="SHE",
            )
            fits.HDUList([fits.PrimaryHDU(), bkg]).writeto(datadir / vis_detector_bkg_filename)
            logger.info('Wrote VIS detector background file [filename=%s]', vis_detector_bkg_filename)

            # VIS WGT
            vis_detector_wgt_filename = filename_provider.get_allowed_filename(
                type_name="WGT",
                instance_id=instance_id,
                extension=".fits",
                release=SHE_CTE_RELEASE_STRING,
                processing_function="SHE",
            )
            fits.HDUList([fits.PrimaryHDU(), wgt]).writeto(datadir / vis_detector_wgt_filename)
            logger.info('Wrote VIS detector weights file [filename=%s]', vis_detector_wgt_filename)

            # VIS xml
            det_xml = deepcopy(vis_product)
            det_xml.Data.DetectorList.Detector = [detector]
            det_xml.Data.WLQ2Exposures.PRNUCorrectedExposure.DataContainer.FileName = vis_detector_det_filename
            det_xml.Data.DataStorage.DataContainer.FileName = vis_detector_det_filename
            det_xml.Data.BackgroundStorage.DataContainer.FileName = vis_detector_bkg_filename
            det_xml.Data.WeightStorage.DataContainer.FileName = vis_detector_wgt_filename
            det_xml_filename = filename_provider.get_allowed_filename(
                type_name="DET",
                instance_id=instance_id,
                extension=".xml",
                release=SHE_CTE_RELEASE_STRING,
                processing_function="SHE",
            )
            save_product_metadata(det_xml, workdir / det_xml_filename)
            vis_products.append(det_xml_filename)
            logger.info('Wrote VIS detector product [filename=%s]', det_xml_filename)

    # Write listfile of VIS detector products
    with open(workdir / args.vis_detector_frames, 'w') as fs:
        json.dump(vis_products, fs)
    logger.info('Wrote listfile of VIS detector products [filename=%s]', args.vis_detector_frames)

    logger.info('#')
    logger.info('# Exiting SHE_CTE_DetectorSplit mainMethod()')
    logger.info('#')

    return


if __name__ == '__main__':
    parser = defineSpecificProgramOptions()
    args = parser.parse_args()
    mainMethod(args)
