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
:file: python/SHE_CTE_PipelineUtility/DetectorCatalog.py

:date: 2024/01/23
:author: @rrollins

"""

from argparse import ArgumentParser, ArgumentTypeError
from copy import deepcopy
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from ElementsKernel.Logging import getLogger
from SHE_PPT.file_io import read_xml_product, save_product_metadata
from SHE_PPT.products import she_object_id_list

from SHE_CTE_PipelineUtility.object_id_split import skycoords_in_wcs
from SHE_CTE_PipelineUtility.utils import write_fits_table


MER_DETECTOR_CATALOG_COLUMNS = ['OBJECT_ID', 'RIGHT_ASCENSION', 'DECLINATION']
WCS_HDU = 1


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
    parser.add_argument('--vis_detector_frame', type=str, help='VIS detector data product')
    parser.add_argument('--mer_final_catalog', type=str, help='MER final catalog data product')
    parser.add_argument('--mer_detector_catalog', type=str, help='MER detector catalog data product')
    parser.add_argument('--she_object_id_list', type=str, help='SHE Object ID list data product')

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger('SHE_CTE_DetectorCatalog')

    logger.info('#')
    logger.info('# Entering SHE_CTE_DetectorCatalog mainMethod()')
    logger.info('#')

    # Paths
    workdir = Path(args.workdir)
    datadir = workdir / 'data'

    # Read VIS product
    vis_product = read_xml_product(workdir / args.vis_detector_frame)
    if len(vis_product.Data.DetectorList.Detector) != 1:
        raise ValueError('SHE_CTE_DetectorCatalog can only run on preprocessed '
                         'VIS products containing data for a single detector')
    observation_id = vis_product.Data.ObservationSequence.ObservationId
    dither_id = vis_product.Data.ObservationSequence.DitherObservation
    exposure_id = vis_product.Data.ObservationSequence.Exposure
    pointing_id = vis_product.Data.ObservationSequence.PointingId
    detector_id = vis_product.Data.DetectorList.Detector[0].DetectorId
    vis_det_file = datadir / vis_product.Data.DataStorage.DataContainer.FileName

    # VIS detector WCS
    with fits.open(vis_det_file) as vis_hdulist:
        if len(vis_hdulist) != 4:
            raise ValueError('SHE_CTE_DetectorCatalog can only run on preprocessed '
                             'VIS products containing data for a single detector')
        detector_wcs = WCS(vis_hdulist[WCS_HDU].header)

    # Read MER product
    mer_catalog_product = read_xml_product(workdir / args.mer_final_catalog)
    tile_id = mer_catalog_product.Data.TileIndex
    instance_id = f'{tile_id}_{observation_id:06}_{dither_id:02}_{exposure_id}_{detector_id}'
    mer_cat_file = datadir / mer_catalog_product.Data.DataStorage.DataContainer.FileName
    mer_catalog = Table.read(mer_cat_file)
    mer_catalog.keep_columns(MER_DETECTOR_CATALOG_COLUMNS)
    ra = mer_catalog['RIGHT_ASCENSION']
    dec = mer_catalog['DECLINATION']
    sky_coordinates = SkyCoord(ra, dec, unit='degree', frame='icrs')

    # Mask objects outside the VIS detector footprint and write catalog to FITS
    in_detector = skycoords_in_wcs(sky_coordinates, detector_wcs)
    detector_catalog = mer_catalog[in_detector]
    detector_catalog_filename = write_fits_table('DET-CAT', instance_id, detector_catalog, datadir)
    logger.info('Wrote detector catalog FITS file [filename=%s]', detector_catalog_filename)

    # MER batch catalog xml data product
    mer_xml = deepcopy(mer_catalog_product)
    mer_xml.Data.DataStorage.DataContainer.FileName = detector_catalog_filename
    save_product_metadata(mer_xml, workdir / args.mer_detector_catalog)
    logger.info('Wrote MER detector catalog product [filename=%s]', args.mer_detector_catalog)

    # Object ID list xml
    id_list = list(detector_catalog['OBJECT_ID'])
    dpd = she_object_id_list.create_dpd_she_object_id_list(id_list=id_list)
    dpd.Data.BatchIndex = 0
    dpd.Data.PointingIdList = [pointing_id]
    dpd.Data.ObservationIdList = [observation_id]
    dpd.Data.TileList = [tile_id]
    save_product_metadata(dpd, workdir / args.she_object_id_list)
    logger.info('Wrote SHE object ID list product [filename=%s]', args.she_object_id_list)

    logger.info('#')
    logger.info('# Exiting SHE_CTE_DetectorCatalog mainMethod()')
    logger.info('#')

    return


if __name__ == '__main__':
    parser = defineSpecificProgramOptions()
    args = parser.parse_args()
    mainMethod(args)
