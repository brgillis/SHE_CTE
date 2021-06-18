""" @file object_id_split.py

    Created 14 Mar 2019

    Functions to handle split over object IDs in a MER catalog
"""

__updated__ = "2021-05-27"

# Copyright (C) 2012-2020 Euclid Science Ground Segment
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

from copy import deepcopy
import math
import os

from SHE_PPT.file_io import (write_listfile,
                             write_xml_product,
                             get_allowed_filename)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import read_analysis_config, AnalysisConfigKeys
from SHE_PPT.products import she_object_id_list, mer_final_catalog
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.table_formats.mer_final_catalog import initialise_mer_final_catalog, tf as mfc_tf

import SHE_CTE
import numpy as np

default_batch_size = 400
default_max_batches = 0

default_sub_batch_size = 20
default_max_sub_batches = 0


# Pipeline config keys for the outer loop (sub_batch=False)
defaults_dict_batch = {AnalysisConfigKeys.OID_BATCH_SIZE.value: default_batch_size,
                       AnalysisConfigKeys.OID_MAX_BATCHES.value: default_max_batches,
                       AnalysisConfigKeys.OID_IDS.value: None,
                       AnalysisConfigKeys.PIP_PROFILE.value: "False"}

# Pipeline config keys for the inner loop (sub_batch=True)
defaults_dict_sub_batch = {AnalysisConfigKeys.SOID_BATCH_SIZE.value: default_sub_batch_size,
                           AnalysisConfigKeys.SOID_MAX_BATCHES.value: default_max_sub_batches,
                           AnalysisConfigKeys.PIP_PROFILE.value: "False"}

logger = getLogger(__name__)


def read_oid_split_pipeline_config(pipeline_config,
                                   workdir,
                                   sub_batch=False):

    if not sub_batch:
        defaults_dict = defaults_dict_batch
        batch_size_key = AnalysisConfigKeys.OID_BATCH_SIZE.value
        max_batches_key = AnalysisConfigKeys.OID_MAX_BATCHES.value
        ids_key = AnalysisConfigKeys.OID_IDS.value

    else:
        defaults_dict = defaults_dict_sub_batch
        batch_size_key = AnalysisConfigKeys.SOID_BATCH_SIZE.value
        max_batches_key = AnalysisConfigKeys.SOID_MAX_BATCHES.value
        ids_key = None

    if type(pipeline_config) is str:
        pipeline_config = read_analysis_config(config_filename=pipeline_config,
                                            workdir=workdir,
                                            cline_args=None,
                                            defaults=defaults_dict)
    elif type(pipeline_config) is dict:
        pass
    else:
        raise TypeError("read_oid_split_pipeline_config: pipeline_config is an unexpected type - %s"%type(pipeline_config))


    # Convert values into the proper type and log values

    batch_size = int(pipeline_config[batch_size_key])
    logger.info(f"Using batch size of: {batch_size}")

    max_batches = int(pipeline_config[max_batches_key])
    logger.info(f"Using max batches of: {max_batches}")

    if ids_key is not None:
        ids_to_use = pipeline_config[ids_key]
    else:
        ids_to_use = None

    # Check for the IDs key
    if ids_to_use is None:
        logger.info("Using all IDs")
    else:

        ids_to_use_str = pipeline_config[ids_key].split()

        if len(ids_to_use_str) == 0:
            ids_to_use = None
            logger.info("Using all IDs")
        else:
            ids_to_use = list(map(int, ids_to_use_str))
            logger.info(f"Using limited selection of IDs: {ids_to_use}")

    return ids_to_use, max_batches, batch_size


def get_tile_list(data_stack):

    tile_list = []
    for detections_catalogue_product in data_stack.detections_catalogue_products:
        tile_index = detections_catalogue_product.Data.TileIndex
        tile_product_id = detections_catalogue_product.Header.ProductId
        tile_list.append((tile_index, tile_product_id))

    return tile_list


def get_pointing_list_and_observation_id(data_stack):
    pointing_list = []
    observation_id = None
    for exposure_product in data_stack.exposure_products:
        pointing_list.append(exposure_product.Data.ObservationSequence.PointingId)
        if observation_id is None:
            observation_id = exposure_product.Data.ObservationSequence.ObservationId
        else:
            new_observation_id = exposure_product.Data.ObservationSequence.ObservationId
            if not observation_id == new_observation_id:
                # Check it's consistent and warn if not
                logger.warning(f"Inconsistent observation IDs: {observation_id} and {new_observation_id}.")

    return pointing_list, observation_id


def get_ids_array(ids_to_use, max_batches, batch_size, data_stack, sub_batch=False):

    if ids_to_use is not None:

        all_ids_array = np.array(ids_to_use)

    else:

        # Get the number of IDs desired - 0 indicates all available
        num_ids_desired = max_batches * batch_size

        all_ids = set(data_stack.detections_catalogue[mfc_tf.ID].data)

        logger.info("Finished reading in IDs from mer_final_catalog.")

        # Prune IDs that aren't in the images
        good_ids = []
        num_good_ids = 0

        for gal_id in all_ids:

            # If this is a sub-batch, checking has already been done so we can trust all IDs are good
            if sub_batch:
                good_ids.append(gal_id)
                num_good_ids += 1
                continue

            # Get a stack of the galaxy images
            gal_stamp_stack = data_stack.extract_galaxy_stack(gal_id=gal_id,
                                                              width=1,)

            # Do we have any data for this object?
            if not gal_stamp_stack.is_empty():
                good_ids.append(gal_id)
                num_good_ids += 1
                if num_ids_desired > 0 and num_good_ids >= num_ids_desired:
                    break

        all_ids_array = np.array(good_ids)

    return all_ids_array


def read_oid_input_data(data_images,
                        mer_final_catalog_tables,
                        workdir,
                        ids_to_use,
                        max_batches,
                        batch_size,
                        sub_batch=False):

    # Read in the data images
    logger.info("Reading data images...")

    data_stack = SHEFrameStack.read(exposure_listfile_filename=data_images,
                                    detections_listfile_filename=mer_final_catalog_tables,
                                    workdir=workdir,
                                    save_products=True,
                                    memmap=True,
                                    mode='denywrite',
                                    load_images=False)

    first_mer_final_catalog_product = data_stack.detections_catalogue_products[0]

    # Get the tile list, pointing list, and observation id from the input data products

    tile_list = get_tile_list(data_stack)

    pointing_list, observation_id = get_pointing_list_and_observation_id(data_stack)

    all_ids_array = get_ids_array(ids_to_use, max_batches, batch_size, data_stack, sub_batch=sub_batch)

    num_ids = len(all_ids_array)

    # If batch size is zero, use all IDs in one batch
    if batch_size == 0:
        batch_size = num_ids

    num_batches = int(math.ceil(num_ids / batch_size))

    if max_batches > 0:
        limited_num_batches = np.min((num_batches, max_batches))
    else:
        limited_num_batches = num_batches

    logger.info(f"Splitting IDs into {limited_num_batches} batches of size {batch_size}.")

    id_split_indices = np.linspace(batch_size, (num_batches - 1) * batch_size,
                                   num_batches - 1, endpoint=True, dtype=int)

    id_arrays = np.split(all_ids_array, id_split_indices)

    # Perform some quick sanity checks
    assert len(id_arrays) == num_batches
    assert len(id_arrays[0]) == batch_size or num_batches == 1
    assert len(id_arrays[-1]) <= batch_size

    return (limited_num_batches,
            id_arrays,
            pointing_list,
            observation_id,
            tile_list,
            data_stack,
            first_mer_final_catalog_product)


def write_oid_batch(workdir,
                    id_arrays,
                    pointing_list,
                    observation_id,
                    tile_list,
                    data_stack,
                    first_mer_final_catalog_product,
                    i):

    # For the filenames, we want to set it up in a subfolder so we don't get too many files
    subfolder_number = i % 256
    subfolder_name = f"data/s{subfolder_number}"

    qualified_subfolder_name = os.path.join(workdir, subfolder_name)

    if not os.path.exists(qualified_subfolder_name):
        # Can we create it?
        try:
            os.mkdir(qualified_subfolder_name)
        except Exception as e:
            logger.error(f"Directory {qualified_subfolder_name} does not exist and cannot be created.")
            raise e

    logger.debug(f"Writing ID list #{i} to product.")

    # Get a filename for this batch and store it in the list
    batch_id_list_product_filename = get_allowed_filename(type_name="OBJ-ID-LIST",
                                                          instance_id=str(i),
                                                          extension=".xml",
                                                          version=SHE_CTE.__version__,
                                                          subdir=subfolder_name,
                                                          processing_function="SHE")

    # Create the product and fill out metadata
    batch_id_list_product = she_object_id_list.create_dpd_she_object_id_list(id_list=list(id_arrays[i]))

    batch_id_list_product.Data.BatchIndex = i
    batch_id_list_product.Data.PointingIdList = pointing_list
    batch_id_list_product.Data.ObservationId = observation_id

    num_tiles = len(tile_list)

    try:
        base_tile_object = deepcopy(batch_id_list_product.Data.TileList[0])
        batch_id_list_product.Data.TileList = [base_tile_object] * num_tiles
        tile_list_binding = batch_id_list_product.Data.TileList

        for tile_i, (tile_index, tile_product_id) in enumerate(tile_list):
            tile_list_binding[tile_i].TileIndex = tile_index
            tile_list_binding[tile_i].TileProductId = tile_product_id
    except TypeError as e:
        if not "object does not support indexing" in str(e):
            raise
        logger.warning("Cannot list all tiles in data product; will only list first tile.")
        tile_list_binding = batch_id_list_product.Data.TileList
        tile_list_binding.TileIndex = tile_list[0][0]
        tile_list_binding.TileProductId = tile_list[0][1]

    # Save the product
    write_xml_product(batch_id_list_product, batch_id_list_product_filename, workdir=workdir)

    logger.debug(f"Successfully wrote ID list #{i} to product: {batch_id_list_product_filename}")

    # Create and write out the batch catalog

    # Get a filename for the batch catalog
    batch_mer_catalog_filename = get_allowed_filename(type_name="BATCH-MER-CAT",
                                                      instance_id=str(i),
                                                      extension=".fits",
                                                      version=SHE_CTE.__version__,
                                                      subdir=subfolder_name,
                                                      processing_function="SHE")

    # Init the catalog and copy over metadata
    batch_mer_catalog = initialise_mer_final_catalog(optional_columns=data_stack.detections_catalogue.colnames)
    for key in data_stack.detections_catalogue.meta:
        batch_mer_catalog.meta[key] = data_stack.detections_catalogue.meta[key]

    # Add a row to the catalog for each ID, taking the row from the combined catalogue
    for object_id in id_arrays[i]:
        object_row = data_stack.detections_catalogue.loc[object_id]
        batch_mer_catalog.add_row(object_row)

    # Write out the catalog
    batch_mer_catalog.write(os.path.join(workdir, batch_mer_catalog_filename))

    logger.debug(f"Successfully wrote batch MER catalog #{i} to: {batch_mer_catalog_filename}")

    # Create and write out the batch MER product
    logger.debug(f"Writing batch MER catalog product #{i}.")

    # Get a filename for this product and store it in the list
    batch_mer_catalog_product_filename = get_allowed_filename(type_name="P-BATCH-MER-CAT",
                                                              instance_id=str(i),
                                                              extension=".xml",
                                                              version=SHE_CTE.__version__,
                                                              subdir=subfolder_name,
                                                              processing_function="SHE")

    # Create the product and copy metadata from the first catalogue
    batch_mer_catalog_product = mer_final_catalog.create_dpd_mer_final_catalog()

    for attr in ['CatalogDescription', 'CutoutsCatalogStorage', 'DataStorage', 'ObservationIdList', 'ProcessingMode',
                 'ProcessingSteps', 'QualityParams', 'SpatialCoverage', 'SpectralCoverage', 'TileIndex']:
        setattr(batch_mer_catalog_product.Data, attr, getattr(first_mer_final_catalog_product.Data, attr))

    # Overwrite the data filename with the new catalog
    batch_mer_catalog_product.set_data_filename(batch_mer_catalog_filename)

    # Save the product
    write_xml_product(batch_mer_catalog_product, batch_mer_catalog_product_filename, workdir=workdir)

    logger.debug(f"Successfully wrote batch MER catalog product #{i} to: {batch_mer_catalog_product_filename}")

    return (batch_id_list_product_filename,
            batch_mer_catalog_product_filename)


def object_id_split_from_args(args,
                              sub_batch=False):
    """ Core function for implementing a split by object ID
    """

    logger.debug('# Entering object_id_split_from_args(args)')

    (ids_to_use,
     max_batches,
     batch_size) = read_oid_split_pipeline_config(pipeline_config=args.pipeline_config,
                                                  workdir=args.workdir,
                                                  sub_batch=sub_batch)

    (limited_num_batches,
     id_arrays,
     pointing_list,
     observation_id,
     tile_list,
     data_stack,
     first_mer_final_catalog_product) = read_oid_input_data(data_images=args.data_images,
                                                            mer_final_catalog_tables=args.mer_final_catalog_tables,
                                                            workdir=args.workdir,
                                                            ids_to_use=ids_to_use,
                                                            max_batches=max_batches,
                                                            batch_size=batch_size,
                                                            sub_batch=sub_batch)

    # Start outputting the ID lists for each batch and creating trimmed catalogs

    # Keep a list of all product filenames
    id_list_product_filename_list = []
    batch_mer_catalog_product_filename_list = []

    logger.info("Writing ID lists and batch catalogs into products.")

    for i in range(limited_num_batches):

        (batch_id_list_product_filename,
         batch_mer_catalog_product_filename) = write_oid_batch(workdir=args.workdir,
                                                               id_arrays=id_arrays,
                                                               pointing_list=pointing_list,
                                                               observation_id=observation_id,
                                                               tile_list=tile_list,
                                                               data_stack=data_stack,
                                                               first_mer_final_catalog_product=first_mer_final_catalog_product,
                                                               i=i)

        id_list_product_filename_list.append(batch_id_list_product_filename)
        batch_mer_catalog_product_filename_list.append(batch_mer_catalog_product_filename)

    # Output the listfiles
    write_listfile(os.path.join(args.workdir, args.object_ids), id_list_product_filename_list)
    logger.info(f"Finished writing listfile of object ID list products to {args.object_ids}")

    write_listfile(os.path.join(args.workdir, args.batch_mer_catalogs), batch_mer_catalog_product_filename_list)
    logger.info(f"Finished writing listfile of batch MER catalog products to {args.batch_mer_catalogs}")

    logger.debug('# Exiting object_id_split_from_args normally')

    return
