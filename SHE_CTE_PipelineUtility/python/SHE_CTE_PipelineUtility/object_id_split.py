""" @file object_id_split.py

    Created 14 Mar 2019

    Functions to handle split over object IDs in a MER catalog
"""

__updated__ = "2021-08-19"

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

import math
import os
from copy import deepcopy
from typing import List, Optional, Sequence, Set
import time

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.units import deg as degrees
from astropy.table import Column

import SHE_CTE
from SHE_PPT.file_io import (get_allowed_filename, write_listfile, write_xml_product)
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import AnalysisConfigKeys
from SHE_PPT.products import mer_final_catalog, she_object_id_list
from SHE_PPT.she_frame_stack import SHEFrameStack
from SHE_PPT.she_image_stack import SHEImageStack
from SHE_PPT.table_formats.mer_final_catalog import initialise_mer_final_catalog, tf as mfc_tf
from SHE_PPT.coordinates import reproject_to_equator, euclidean_metric, DTOR, RTOD
from SHE_PPT.clustering import identify_all_blends, partition_into_batches

from SHE_PPT.table_utility import is_in_format

logger = getLogger(__name__)


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


def get_ids_array(ids_to_use: Optional[Sequence[int]],
                  max_batches: int,
                  batch_size: int,
                  data_stack: SHEFrameStack,
                  sub_batch: bool = False):
    """Gets the list of all valid object IDs from the catalogue that are in the observation, and the indices of these objects in the catalogue
       if ids_to_use is not None, get the valid object IDs and indices from this list rather than the whole catalogue"""

    s_all_ids: Set[int]
    ids_supplied: bool
    num_ids_desired: int
    
    #get the list of object ids we want, and their indices in the catalogue
    if ids_to_use is not None:

        ids_supplied = True

        # Set the number of IDs desired - 0 indicates all available
        num_ids_desired = 0

        #get the unique objects in the catalogue and their indices
        catalogue_ids, unique_catalogue_indices = np.unique(data_stack.detections_catalogue[mfc_tf.ID].data,return_index=True)

        #get the unique ids_to_use
        ids_to_use = np.unique(np.asarray(ids_to_use))

        #get the indices of the catalogue object that are in ids_to_use (e.g. we don't want ids that are not in the catalogue)
        ids_to_use_indices = np.where(np.isin(catalogue_ids,ids_to_use))
        
        #now get the indices of these points in the original catalogue, and their ids
        indices_in_catalogue = unique_catalogue_indices[ids_to_use_indices]
        s_all_ids = catalogue_ids[ids_to_use_indices]

        
        # #select only the ids_to_use that are in the catalogue, removing duplicates
        # catalogue_ids = data_stack.detections_catalogue[mfc_tf.ID].data
        # valid_ids = np.unique(catalogue_ids[np.where(np.isin(catalogue_ids,ids_to_use))])

        # # s_all_ids = set(np.array(ids_to_use))
        # s_all_ids, unique_indices = np.unique(np.array(ids_to_use),return_index=True)

        logger.info("Using input list of IDs.")

    else:

        ids_supplied = False

        # Get the number of IDs desired - 0 indicates all available
        num_ids_desired = max_batches * batch_size

        #s_all_ids = set(data_stack.detections_catalogue[mfc_tf.ID].data)
        #get the unique object ids in the catalogue, and their indices
        s_all_ids, indices_in_catalogue = np.unique(data_stack.detections_catalogue[mfc_tf.ID].data,return_index=True)

        logger.info("Finished reading in IDs from mer_final_catalog.")

   
    #get the objects in the observation
    
    if sub_batch:
        # If this is a sub-batch, checking has already been done so we can trust all IDs are good
        all_ids_indices = indices
        all_ids_array = s_all_ids
    else:
        t0 = time.time()
        logger.info("Extracting objecs in the observation from the catalogue")

        #first get the sky coordinates of all the objects in the candidate list, convert to SkyCoords
        ras = data_stack.detections_catalogue[mfc_tf.gal_x_world].data[indices_in_catalogue]
        decs = data_stack.detections_catalogue[mfc_tf.gal_y_world].data[indices_in_catalogue]
        object_coords = SkyCoord(ras,decs,unit=degrees)
        
        #get the indices of all the objects found in the observation (these are relative to the indices in s_all_ids)
        found_object_indices = data_stack.get_objects_in_observation(object_coords)
        
        #get the (catalogue) indices and ids for these objects
        all_ids_indices = indices_in_catalogue[found_object_indices]
        all_ids_array = s_all_ids[found_object_indices]

        t1 = time.time()
        logger.info(f"Extracted {len(all_ids_array)} objects in observation from catalogue in {t1-t0}s")

    
    #if we don't want all the IDs, trim the list down to the specified size
    if num_ids_desired > 0:
        all_ids_array = all_ids_array[0:num_ids_desired]
        all_ids_indices = all_ids_indices[0:num_ids_desired]


    return all_ids_array, all_ids_indices


def read_oid_input_data(data_images,
                        mer_final_catalog_tables,
                        workdir,
                        ids_to_use,
                        max_batches,
                        batch_size,
                        sub_batch = False):
    # Read in the data images
    logger.info("Reading data images...")

    data_stack = SHEFrameStack.read(exposure_listfile_filename = data_images,
                                    detections_listfile_filename = mer_final_catalog_tables,
                                    workdir = workdir,
                                    save_products = True,
                                    memmap = (not sub_batch),
                                    mode = 'denywrite',
                                    load_images = (not sub_batch))

    first_mer_final_catalog_product = data_stack.detections_catalogue_products[0]

    # Get the tile list, pointing list, and observation id from the input data products

    tile_list = get_tile_list(data_stack)

    pointing_list, observation_id = get_pointing_list_and_observation_id(data_stack)


    #get the list of object IDs (and their indices in the catalogue) for objects in the observation
    all_ids_array, all_ids_indices = get_ids_array(ids_to_use = ids_to_use,
                                                    max_batches = max_batches,
                                                    batch_size = batch_size,
                                                    data_stack = data_stack,
                                                    sub_batch = sub_batch)

    num_ids = len(all_ids_array)

    # If batch size is zero, use all IDs in one batch
    if batch_size == 0:
        batch_size = num_ids

    num_batches = int(math.ceil(num_ids / batch_size))

    if max_batches > 0:
        limited_num_batches = np.min((num_batches, max_batches))
    else:
        limited_num_batches = num_batches

    
    #get ras and decs of the objects
    ras = data_stack.detections_catalogue[mfc_tf.gal_x_world].data[all_ids_indices]
    decs = data_stack.detections_catalogue[mfc_tf.gal_y_world].data[all_ids_indices]
    
    #  change coordinate system to one centred about the celestial equator so that over the half degree or so
    #  range we can approximate ra/dec as cartesian
    ras_eq, decs_eq = reproject_to_equator(ras,decs)

    logger.info("Identifying blended objects...")

    #  first we want to identify blends. In order to do this, we want to further improve the coordinates such
    #  that they're as close to cartesian as possible. To do this, remember that an small distance in
    #  spherical coordinates is delta l**2 = (delta dec)**2 + (cos(dec) * delta lon)**2. As we're only interested
    #  in accurate small distances (for the blend) we can use a cartesian coordinate system where
    #  x = cos(dec)*ra
    #  y = dec
    
    # calculate the new x and y coordinates, converting to arcseconds
    x = np.cos(decs_eq*DTOR)*ras_eq * 3600.
    y = decs_eq * 3600.
    
    #get the blends (separation <= 1"), and the updated positions of all objects for batching (all blends have the same position)
    x_blend, y_blend, blend_ids = identify_all_blends(x, y, sep=1, metric=euclidean_metric)

    logger.info(f"Splitting objects into batches of mean size {batch_size}")
    
    #spatially batch objects into batches of mean size batch_size
    clusters, batch_ids, ns = partition_into_batches(x_blend, y_blend, batchsize = batch_size)

    
    #now construct the lists of IDs and indices for each batch

    id_arrays = []
    index_arrays = []
    blend_arrays=[]

    unique_batch_ids = set(batch_ids)

    n=1
    for batch in unique_batch_ids:
        if n > limited_num_batches:
            #if we've reached the maximum number of batches required, we are done
            break

        inds = np.where(batch_ids == batch)

        ids = all_ids_array[inds]
        indices = all_ids_indices[inds]
        blends = blend_ids[inds]

        id_arrays.append(ids)
        index_arrays.append(indices)
        blend_arrays.append(blends)

        n+=1

    #make sure we do have num_batches batches
    assert len(id_arrays) == num_batches


    return (limited_num_batches,
            id_arrays,
            blend_arrays,
            index_arrays,
            pointing_list,
            observation_id,
            tile_list,
            data_stack,
            first_mer_final_catalog_product)


def write_oid_batch(workdir,
                    id_arrays,
                    blend_arrays,
                    index_arrays,
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
    batch_id_list_product_filename = get_allowed_filename(type_name = "OBJ-ID-LIST",
                                                          instance_id = f"{os.getpid()}-{i}",
                                                          extension = ".xml",
                                                          version = SHE_CTE.__version__,
                                                          subdir = subfolder_name,
                                                          processing_function = "SHE")

    # Create the product and fill out metadata
    batch_id_list_product = she_object_id_list.create_dpd_she_object_id_list(id_list = list(id_arrays[i]))

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
        if "object does not support indexing" not in str(e):
            raise
        logger.warning("Cannot list all tiles in data product; will only list first tile.")
        tile_list_binding = batch_id_list_product.Data.TileList
        tile_list_binding.TileIndex = tile_list[0][0]
        tile_list_binding.TileProductId = tile_list[0][1]

    # Save the product
    write_xml_product(batch_id_list_product, batch_id_list_product_filename, workdir = workdir)

    logger.debug(f"Successfully wrote ID list #{i} to product: {batch_id_list_product_filename}")

    # Create and write out the batch catalog

    # Get a filename for the batch catalog
    batch_mer_catalog_filename = get_allowed_filename(type_name = "BATCH-MER-CAT",
                                                      instance_id = f"{os.getpid()}-{i}",
                                                      extension = ".fits",
                                                      version = SHE_CTE.__version__,
                                                      subdir = subfolder_name,
                                                      processing_function = "SHE")
    
    #get the data from the detections_catalogue for this batch, creating a new table
    # (The filled() method creates a copy with this data only)
    batch_table = data_stack.detections_catalogue[index_arrays[i]].filled()
    
    #IF the group id is not already in the table, add it
    if mfc_tf.GROUP_ID not in batch_table.columns:
        col = Column(name = mfc_tf.GROUP_ID, data = blend_arrays[i])
        batch_table.add_column(col)

    
    # Init the catalog and copy over metadata and the batch's data
    #batch_mer_catalog = initialise_mer_final_catalog(optional_columns = data_stack.detections_catalogue.colnames, init_cols = batch_table.columns)
    batch_mer_catalog = batch_table

    for key in data_stack.detections_catalogue.meta:
        batch_mer_catalog.meta[key] = data_stack.detections_catalogue.meta[key]
    
    # This could be included, but it ads 0.5s time per table... 150s alltogether on runtime
    assert(is_in_format(batch_mer_catalog, mfc_tf, verbose=True))

    # Write out the catalog
    batch_mer_catalog.write(os.path.join(workdir, batch_mer_catalog_filename))

    logger.debug(f"Successfully wrote batch MER catalog #{i} to: {batch_mer_catalog_filename}")

    # Create and write out the batch MER product
    logger.debug(f"Writing batch MER catalog product #{i}.")

    # Get a filename for this product and store it in the list
    batch_mer_catalog_product_filename = get_allowed_filename(type_name = "P-BATCH-MER-CAT",
                                                              instance_id = f"{os.getpid()}-{i}",
                                                              extension = ".xml",
                                                              version = SHE_CTE.__version__,
                                                              subdir = subfolder_name,
                                                              processing_function = "SHE")

    # Create the product and copy metadata from the first catalogue
    batch_mer_catalog_product = mer_final_catalog.create_dpd_mer_final_catalog()

    for attr in ['CatalogDescription', 'CutoutsCatalogStorage', 'DataStorage', 'ObservationIdList', 'ProcessingMode',
                 'ProcessingSteps', 'QualityParams', 'SpatialCoverage', 'SpectralCoverage', 'TileIndex']:
        setattr(batch_mer_catalog_product.Data, attr, getattr(first_mer_final_catalog_product.Data, attr))

    # Overwrite the data filename with the new catalog
    batch_mer_catalog_product.set_data_filename(batch_mer_catalog_filename)

    # Save the product
    write_xml_product(batch_mer_catalog_product, batch_mer_catalog_product_filename, workdir = workdir)

    logger.debug(f"Successfully wrote batch MER catalog product #{i} to: {batch_mer_catalog_product_filename}")

    # Write a listfile for this product for a consistent interface
    batch_mer_catalog_listfile_filename = get_allowed_filename(type_name = "L-BATCH-MER-CAT",
                                                               instance_id = f"{os.getpid()}-{i}",
                                                               extension = ".json",
                                                               version = SHE_CTE.__version__,
                                                               subdir = subfolder_name,
                                                               processing_function = "SHE")
    write_listfile(os.path.join(workdir, batch_mer_catalog_listfile_filename), [batch_mer_catalog_product_filename])

    return (batch_id_list_product_filename,
            batch_mer_catalog_listfile_filename)


def object_id_split_from_args(args,
                              sub_batch = False):
    """ Core function for implementing a split by object ID
    """

    logger.debug('# Entering object_id_split_from_args(args)')

    if not sub_batch:
        ids_to_use = args.pipeline_config[AnalysisConfigKeys.OID_IDS]
        max_batches = args.pipeline_config[AnalysisConfigKeys.OID_MAX_BATCHES]
        batch_size = args.pipeline_config[AnalysisConfigKeys.OID_BATCH_SIZE]
    else:
        ids_to_use = args.pipeline_config[AnalysisConfigKeys.SOID_IDS]
        max_batches = args.pipeline_config[AnalysisConfigKeys.SOID_MAX_BATCHES]
        batch_size = args.pipeline_config[AnalysisConfigKeys.SOID_BATCH_SIZE]

    (limited_num_batches,
     id_arrays,
     blend_arrays,
     index_arrays,
     pointing_list,
     observation_id,
     tile_list,
     data_stack,
     first_mer_final_catalog_product) = read_oid_input_data(data_images = args.data_images,
                                                            mer_final_catalog_tables = args.mer_final_catalog_tables,
                                                            workdir = args.workdir,
                                                            ids_to_use = ids_to_use,
                                                            max_batches = max_batches,
                                                            batch_size = batch_size,
                                                            sub_batch = sub_batch)

    # Start outputting the ID lists for each batch and creating trimmed catalogs

    # Keep a list of all product filenames
    id_list_product_filename_list = []
    batch_mer_catalog_listfile_filename_list = []

    logger.info("Writing ID lists and batch catalogs into products.")

    for i in range(limited_num_batches):
        t0 = time.time()

        (batch_id_list_product_filename,
         batch_mer_catalog_listfile_filename) = write_oid_batch(
            workdir = args.workdir,
            id_arrays = id_arrays,
            blend_arrays = blend_arrays,
            index_arrays = index_arrays,
            pointing_list = pointing_list,
            observation_id = observation_id,
            tile_list = tile_list,
            data_stack = data_stack,
            first_mer_final_catalog_product = first_mer_final_catalog_product,
            i = i)

        id_list_product_filename_list.append(batch_id_list_product_filename)
        batch_mer_catalog_listfile_filename_list.append(batch_mer_catalog_listfile_filename)

        t1=time.time()
        logger.info(f"Written batch {i+1} of {limited_num_batches} in {t1-t0} seconds")

    # Output the listfiles
    write_listfile(os.path.join(args.workdir, args.object_ids), id_list_product_filename_list)
    logger.info(f"Finished writing listfile of object ID list products to {args.object_ids}")

    write_listfile(os.path.join(args.workdir, args.batch_mer_catalogs), batch_mer_catalog_listfile_filename_list)
    logger.info(f"Finished writing listfile of batch MER catalog products to {args.batch_mer_catalogs}")

    logger.debug('# Exiting object_id_split_from_args normally')
