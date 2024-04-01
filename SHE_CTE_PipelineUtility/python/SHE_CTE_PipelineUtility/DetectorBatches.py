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
:file: python/SHE_CTE_PipelineUtility/DetectorBatches.py

:date: 2024/04/01
:author: @rrollins

"""

import json
import numpy as np
from argparse import ArgumentParser
from copy import deepcopy
from itertools import islice, repeat
from pathlib import Path

from astropy.table import Table
from ElementsKernel.Logging import getLogger
from SHE_PPT.argument_parser import dir_path
from SHE_PPT.file_io import read_xml_product
from SHE_PPT.products import she_object_id_list
from ST_DM_DmUtils.DqcDmUtils import set_quality_parameters

from SHE_CTE_PipelineUtility.utils import ceiling_division, write_data_product, write_fits_table


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    parser = ArgumentParser()
    parser.add_argument('--workdir', type=dir_path, help='Name of the working directory')
    parser.add_argument('--logdir', type=str, help='Log directory')
    parser.add_argument('--vis_detector_frame', type=str, help='VIS detector frame data product.')
    parser.add_argument('--mer_detector_catalog', type=str, help='MER detector catalog data product')
    parser.add_argument('--vis_detector_batches', type=str, help='Listfile of VIS detector frame products')
    parser.add_argument('--mer_catalog_batches', type=str, help='Listfile of MER batch catalog products')
    parser.add_argument('--she_objectid_lists', type=str, help='Listfile of SHE ObjectIdList products')
    parser.add_argument('--batch_size', type=int, help='Fixed size of each batch')
    parser.add_argument('--max_batches', type=int, help='Limit on the maximum number of batches to produce')

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger('SHE_CTE_DetectorBatches')

    logger.info('#')
    logger.info('# Entering SHE_CTE_DetectorBatches mainMethod()')
    logger.info('#')

    # Paths
    workdir = Path(args.workdir)
    datadir = workdir / 'data'

    # Read VIS product
    vis_product = read_xml_product(workdir / args.vis_detector_frame)
    if len(vis_product.Data.DetectorList.Detector) != 1:
        raise ValueError('')
    observation_id = vis_product.Data.ObservationSequence.ObservationId
    dither_id = vis_product.Data.ObservationSequence.DitherObservation
    exposure_id = vis_product.Data.ObservationSequence.Exposure
    pointing_id = vis_product.Data.ObservationSequence.PointingId
    detector_id = vis_product.Data.DetectorList.Detector[0].DetectorId

    # Read MER product
    mer_catalog_product = read_xml_product(workdir / args.mer_detector_catalog)
    tile_id = mer_catalog_product.Data.TileIndex
    mer_cat_file = datadir / mer_catalog_product.Data.DataStorage.DataContainer.FileName
    mer_catalog = Table.read(mer_cat_file)

    n_objects = len(mer_catalog)
    n_batches = ceiling_division(n_objects, args.batch_size)
    batch_assignment = np.repeat(range(n_batches), args.batch_size)[:n_objects]
    batched_catalog = mer_catalog.group_by(batch_assignment)
    n_batches = min(n_batches, args.max_batches) if args.max_batches else n_batches
    logger.info('Batching objects [n_objects=%s, batch_size=%s, n_batches=%s]', n_objects, args.batch_size, n_batches)

    vis_products = list(repeat(args.vis_detector_frame, n_batches))
    mer_products = list()
    she_products = list()

    for batch_id, batch_catalog in enumerate(islice(batched_catalog.groups, n_batches)):

        # Unique instance ID for batch
        instance_id = f'{tile_id}_{observation_id:06}_{dither_id:02}_{exposure_id}_{detector_id}_{batch_id}'

        # Write batch catalog to FITS file
        batch_catalog_filename = write_fits_table('BATCH-CAT', instance_id, batch_catalog, datadir)
        logger.info('Wrote detector catalog FITS file [filename=%s]', batch_catalog_filename)

        # MER batch catalog xml data product
        mer_xml = deepcopy(mer_catalog_product)
        mer_xml.Data.DataStorage.DataContainer.FileName = batch_catalog_filename
        object_count = {"value": len(batch_catalog), "flagged": False}
        set_quality_parameters(mer_xml.Data.QualityParams, {'ObjectCount': object_count})
        mer_product_filename = write_data_product('BATCH-CAT', instance_id, mer_xml, workdir)
        mer_products.append(mer_product_filename)
        logger.info('Wrote MER detector catalog product [filename=%s]', mer_product_filename)

        # Object ID list xml
        id_list = list(batch_catalog['OBJECT_ID'])
        dpd = she_object_id_list.create_dpd_she_object_id_list(id_list=id_list)
        dpd.Data.BatchIndex = batch_id
        dpd.Data.PointingIdList = [pointing_id]
        dpd.Data.ObservationIdList = [observation_id]
        dpd.Data.TileList = [tile_id]
        id_list_filename = write_data_product('OBJ-ID-LIST', instance_id, dpd, workdir)
        she_products.append(id_list_filename)
        logger.info('Wrote SHE object ID list product [filename=%s]', id_list_filename)

    # Write listfile of VIS detector products
    with open(workdir / args.vis_detector_batches, 'w') as fs:
        json.dump(vis_products, fs)
    logger.info('Wrote listfile of VIS detector products [filename=%s]', args.vis_detector_batches)

    # Write listfile of VIS detector products
    with open(workdir / args.mer_catalog_batches, 'w') as fs:
        json.dump(mer_products, fs)
    logger.info('Wrote listfile of MER batch catalog products [filename=%s]', args.mer_catalog_batches)

    # Write listfile of VIS detector products
    with open(workdir / args.she_objectid_lists, 'w') as fs:
        json.dump(she_products, fs)
    logger.info('Wrote listfile of SheObjectIdList products [filename=%s]', args.she_objectid_lists)

    logger.info('#')
    logger.info('# Exiting SHE_CTE_DetectorBatches mainMethod()')
    logger.info('#')

    return


if __name__ == '__main__':
    parser = defineSpecificProgramOptions()
    args = parser.parse_args()
    mainMethod(args)
