""" @file object_id_split.py

    Created 14 Mar 2019

    Functions to handle split over object IDs in a MER catalog
"""

__updated__ = "2022-09-21"

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

import json
import os
import warnings
import math
from copy import deepcopy
import re

import numpy as np

from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.units import degree
from astropy.coordinates import SkyCoord

from EL_PythonUtils.utilities import hash_any

from ST_DM_FilenameProvider.FilenameProvider import FileNameProvider
from ST_DM_DmUtils.DmUtils import save_product_metadata, read_product_metadata, get_product_name, get_header_from_wcs_binding

from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import AnalysisConfigKeys
from SHE_PPT.coordinates import reproject_to_equator, euclidean_metric, DTOR

from SHE_PPT.table_formats.mer_final_catalog import tf as mer_tf
from SHE_PPT.table_utility import is_in_format

from SHE_PPT.products import she_object_id_list

from SHE_CTE import SHE_CTE_RELEASE_STRING

from .clustering import spatially_batch, find_groups

logger = getLogger(__name__)


def read_vis_frames(vis_frame_listfile, workdir):
    """
    Reads in a number of VIS frames, returning the list of WCSs for each detector, 
    the pointing IDs and the Observation ID

    Inputs:
      - vis_frame_listfile: The name of the listfile pointing to VIS products
      - workdir: The path to the workdir

    Outputs:
      - wcs_list: list of astropy.WCS objects
      - pointing_ids: list of pointing_ids for the exposures
      - observation_id: observation_id of the exposure
    """
    
    qualified_listfile_filename = os.path.join(workdir, vis_frame_listfile)
    if not os.path.exists(qualified_listfile_filename):
        logger.error("VIS frame listfile %s does not exist", qualified_listfile_filename)
        raise FileNotFoundError("VIS frame listfile %s does not exist"%qualified_listfile_filename)

    logger.debug("Opening VIS listfile %s",qualified_listfile_filename)

    with open(qualified_listfile_filename) as f:
        product_list = json.load(f)

    wcs_list = []
    pointing_ids = []
    observation_ids = []

    for product_filename in product_list:
        qualified_product_filename = os.path.join(workdir, product_filename)
        if not os.path.exists(qualified_product_filename):
            logger.error("VIS frame product %s does not exist", qualified_product_filename)
            raise FileNotFoundError("VIS frame product %s does not exist"%qualified_product_filename)
        
        logger.debug("Reading VIS product %s",qualified_product_filename)
        dpd = read_product_metadata(qualified_product_filename)
        product_name = get_product_name(dpd)
        if product_name == "DpdVisCalibratedFrame":
            detectors = dpd.Data.DetectorList.Detector
            quadrants = False
        elif product_name == "DpdVisCalibratedQuadFrame":
            # NOTE: Yes! There is a typo in the Data product definition - Qadrant used here is deliberate!!
            detectors = dpd.Data.QuadrantList.Qadrant
            quadrants = True
        else:
            error_msg = f"{qualified_product_filename} is type {product_name}. Expected DpdVisCalibrated[Quad]Frame"
            logger.warning(error_msg)
            raise ValueError(error_msg)

        pointing_ids.append(dpd.Data.ObservationSequence.PointingId)
        observation_ids.append(dpd.Data.ObservationSequence.ObservationId)

        try:
            # Try to get the WCS list from the xml metadata directly (to avoid reading from FITS)

            exp_wcs_list = []

            nx, ny = dpd.Data.AxisLengths

            for detector in detectors:
                if quadrants:
                    detector_id = detector.QuadrantId
                else:
                    detector_id = detector.DetectorId

                # detector_id should be "x-y[.z - if quadrants]" where x and y range from 0-6 and z is E,F,G or H
                # In mock SHE products, this is a longer string with random contents
                if not re.match('^[1-6]-[1-6](.[EFGH])?$', detector_id):
                    raise ValueError("Detector ID is invalid: %s" % detector_id)

                wcs_hdr = get_header_from_wcs_binding(detector.WCS)
                wcs = WCS(wcs_hdr)
                wcs._naxis = (nx, ny)

                exp_wcs_list.append(wcs)


        except Exception as e:

            # Failsafe - read the WCS from the FITS
            logger.warning("Unable to extract WCS information from the data product with error %s -  will extract from FITS instead", e)


            fits_filename = dpd.Data.DataStorage.DataContainer.FileName
            qualified_fits_filename = os.path.join(workdir,"data",fits_filename)

            if not os.path.exists(qualified_fits_filename):
                logger.error("VIS FITS file %s does not exist", qualified_fits_filename)
                raise FileNotFoundError("VIS FITS file %s does not exist"%qualified_fits_filename)
            
            logger.debug("Reading VIS DET FITS %s",qualified_fits_filename)
            with fits.open(qualified_fits_filename) as hdul:
                # The VIS fits file contains n_detectors x 3 hdus, plus a PrimaryHDU
                # Some SHE mock files (made by SHE_GST) do not contain the PrimaryHDU
                n_hdus = len(hdul)
                offset = n_hdus%3

                det_hdus = hdul[offset::3]

                exp_wcs_list = [ WCS(hdu.header) for hdu in det_hdus]
            
            logger.debug("Extracted %d WCS(s) from %s", len(exp_wcs_list), qualified_fits_filename)
        
        wcs_list += exp_wcs_list

    observation_ids = list(set(observation_ids))
    
    logger.info("ObservationIds = %s", observation_ids)
    logger.info("PointingIds = %s", pointing_ids)
    logger.info("Extracted %d WCS(s) from %d VIS exposure(s)", len(wcs_list), len(product_list))

    return wcs_list, pointing_ids, observation_ids


def read_mer_catalogs(mer_catalog, workdir):
    """
    Reads in the MER final catalog(ue)s from the listfile, returning the concatenated catalog.

    Inputs:
      - mer_catalog: Listfile pointing to MER final catalog data products, or a single MER Final Catalog product
      - workdir: The path to the workdir
    
    Returns:
      - mer_final_catalog: an astropy.Table MER final catalog
      - tile_ids: list of the MER tile ids
      - mfc_dpd: the last dpdMerFinalCatalog data product object read in (used as a template to create output products)
    """
    
    ext = os.path.splitext(mer_catalog)[-1]

    if ext == ".json":
        # Get the list of products from the listfile
        qualified_listfile_name = os.path.join(workdir,mer_catalog)
        if not os.path.exists(qualified_listfile_name):
            logger.error("MER final catalog listfile %s does not exist", qualified_listfile_name)
            raise FileNotFoundError("MER final catalog listfile %s does not exist"%qualified_listfile_name)

        logger.debug("Openeing MER listfile %s", qualified_listfile_name)

        with open(qualified_listfile_name) as f:
            product_list = json.load(f)

    elif ext == ".xml":
        # We have a single product, turn this into a list
        product_list = [mer_catalog]
    
    else:
        # Unknown file type
        raise ValueError("Unknown file extension of mer_final_catalog file: %s"%ext)
    
    mer_catalogues = []
    tile_ids = []

    for product_filename in product_list:
        qualified_product_filename = os.path.join(workdir, product_filename)
        if not os.path.exists(qualified_product_filename):
            logger.error("MER final catalog product %s does not exist", qualified_product_filename)
            raise FileNotFoundError("MER final catalog product %s does not exist"%qualified_product_filename)
        
        logger.debug("Reading MER final catalog product %s", qualified_product_filename)

        dpd = read_product_metadata(qualified_product_filename)
        if get_product_name(dpd) != "DpdMerFinalCatalog":
            logger.error("%s is the wrong product type. Expected DpdMerFinalCatalog - got %s",qualified_product_filename,get_product_name(dpd))
            raise ValueError("%s is the wrong product type. Expected DpdMerFinalCatalog - got %s"%(qualified_product_filename,get_product_name(dpd)))

        tile_ids.append(dpd.Data.TileIndex)

        fits_filename = dpd.Data.DataStorage.DataContainer.FileName

        qualified_fits_filename = os.path.join(workdir,"data",fits_filename)

        if not os.path.exists(qualified_fits_filename):
            logger.error("MER FITS file %s does not exist", qualified_fits_filename)
            raise FileNotFoundError("MER FITS file %s does not exist"%qualified_fits_filename)

        logger.debug("Reading MER Final Catalog table %s", qualified_fits_filename)
        
        # Catch verify warnings caused by the MER tables not complying to the FITS standard
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=VerifyWarning)
            mfc = Table.read(qualified_fits_filename)

        mer_catalogues.append(mfc)

    mer_final_catalog_table = vstack(mer_catalogues, metadata_conflicts="silent")

    logger.info("Extracted %d objects from %d input catalogues",len(mer_final_catalog_table), len(product_list))
    logger.info("Tile IDs = %s", tile_ids)

    return mer_final_catalog_table, tile_ids, dpd

def skycoords_in_wcs(skycoords, wcs):
    """
    Determines which of a set of skycoords are contained within a WCS.

    Whilst one could just use SkyCoord.contained_by or WCS.footprint_contains, there is a bug 
    (fixed in Astropy 5.0.5) with WCS and TPV distortions (which VIS images use) whereby if one
    point fails to converge, all subsequent points also fail. These bad points tend to be very
    far from the detector, so manually excluding points far from the detector before applying
    SkyCoord.contained_by mitigates this error. (This also speeds excution up as we no longer
    need to process as many points!)

    (See https://github.com/astropy/astropy/issues/13750 for details)
    
    Inputs:
      - skycoords: An astropy.coordinates.SkyCoord object, containing one or more object coordinates
      - wcs: An astropy.wcs.WCS object for a detector

    Returns:
      - contained_by: A numpy.ndarray(dtype=Bool) containing True if the object is in the detector, False otherwise
    """
    
    from astropy import __version__ as astropy_version
    from packaging import version

    if version.parse(astropy_version) < version.parse("5.0.5"):
        logger.warning("Astropy version %s is lower than 5.0.5. There is a known WCS bug which may cause incorrect results!", astropy_version)

    ny, nx = wcs.array_shape

    # First determine the coordinates of the centre of the detector, and its corners
    centre_sk = wcs.pixel_to_world(nx/2, ny/2)
    corners = wcs.calc_footprint()
    corners_sk = SkyCoord(corners[:,0], corners[:,1], unit=degree)

    # now get the maximum distance between the centre and the corners (x1.05 buffer)
    max_sep = max(centre_sk.separation(corners_sk)) * 1.05
    
    # We consider points within max_sep from the centre of the detector as candidate 
    # points to determine if they are in the detector's FOV
    candidate_indices = np.where(skycoords.separation(centre_sk) < max_sep)

    # Now detemrine which candidates are within the detector 
    # (can alternatively use WCS.footprint_contains(SkyCoord))
    candidates_contained_by = skycoords[candidate_indices].contained_by(wcs)

    # Create and populate the output array
    contained_by = np.full(len(skycoords), False, dtype=bool)

    contained_by[candidate_indices] = candidates_contained_by

    return contained_by




def prune_objects_not_in_FOV(mer_table, wcs_list, ids_to_use):
    """
    Returns a new mer table with objects not in the FOV (or in ids_to_use) removed

    Inputs:
      - mer_table: A MER Final Catalog table
      - wcs_list: A list of WCS objects for each detector
      - ids_to_use: A list of objects we wish to use. All others will be pruned. If None, use all objects

    Returns:
      - pruned_table: The pruned MER Final Catalog table
    """

    ids = mer_table[mer_tf.ID]

    # prune objects not in ids_to_use
    if ids_to_use:
        inds = [ i for i, id_ in enumerate(ids) if id_ in ids_to_use]
        if len(inds) == 0:
            logger.error("None of the ids_to_use are in the catalogue")
            raise ValueError("None of the ids_to_use are in the catalogue")
        mer_table = mer_table[inds]

    # then check to see if there are any duplicates, remove them if so
    uniques, indices = np.unique(ids, return_index=True)
    if len(uniques) != len(ids):
        logger.warning("%d duplicate objects were detected. Pruning the duplicates", len(ids)-len(uniques))
        mer_table = mer_table[indices]
    
    # construct the SkyCoord object used to check against the WCSs
    ras = mer_table[mer_tf.gal_x_world]
    decs = mer_table[mer_tf.gal_y_world]
    skycoords = SkyCoord(ras, decs, unit=degree)

    n_objs = len(mer_table)
    all_present = np.full(n_objs, False, dtype=bool)
    
    # find the objects in each detector
    for wcs in wcs_list:
        present_in_detector = skycoords_in_wcs(skycoords, wcs)
        if not np.any(present_in_detector):
            logger.warning("No points were found in detector!")
        logger.debug("Found %d objects in detector",np.sum(present_in_detector))
        all_present |= present_in_detector

    pruned_table = mer_table[all_present]

    if len(pruned_table) == 0:
        logger.error("No objects found in the Observation.")
        raise ValueError("No objects found in the Observation.")

    logger.info("Identified %d objects in the FOV", len(pruned_table))

    return pruned_table



def group_and_batch_objects(mer_final_catalog, batch_size, max_batches, grouping_radius):
    """
    Groups objects then splits the objects into batches of mean size batch_size, with a maximum
    number of batches of max_batches, such that all objects belonging to the same group are put
    into the same batch.

    Inputs:
      - mer_final_catalog: The mer_final_catalog
      - batch_size: the mean size of the batches
      - max_batches: the maximum number of batches
      - grouping_radius: the maximum separation between objects belonging to the same group (in arcseconds)

    Outputs:
      - ids_array: list of lists of object_ids for the objects in each batch
      - index_arrays: list of lists of indices for the objects in each batch
      - group_arrays: list of lists of group_ids for the objects in each batch
    """
    

    # First we want to identify blended/grouped objects

    #get ras and decs of the objects
    ras = mer_final_catalog[mer_tf.gal_x_world]
    decs = mer_final_catalog[mer_tf.gal_y_world]

    obj_ids = mer_final_catalog[mer_tf.ID]
    
    #  change coordinate system to one centred about the celestial equator so that over the half degree or so
    #  range we can approximate ra/dec as cartesian
    ras_eq, decs_eq = reproject_to_equator(ras,decs)

    #  Now we want to further improve the coordinates such that they're as close to cartesian as possible. 
    #  To do this, remember that a small distance in spherical coordinates is 
    #  delta l**2 = (delta dec)**2 + (cos(dec) * delta lon)**2. As we're only interested
    #  in accurate small distances (for the blend) we can use a cartesian coordinate system where
    #  x = cos(dec)*ra
    #  y = dec
    
    # calculate the new x and y coordinates, converting to arcseconds
    x = np.cos(decs_eq*DTOR)*ras_eq * 3600.
    y = decs_eq * 3600.

    # get the groups (separation <= grouping_radius), and the updated positions of all objects for batching (all objects
    # belonging to the same group have the same position to ensure they are partitioned into the same batches below)
    x_group, y_group, group_ids = find_groups(x, y, sep=grouping_radius)

    # batch objects with a batch size approximately equal (but generally slightly smaller than) batch_size. 
    # Returns a list of lists of array indices
    batches = spatially_batch(x_group, y_group, batch_size)
    
    # select max_batches batches
    # if max_batches = 0, do not limit the number of batches
    if max_batches:
        if max_batches < 0:
            raise ValueError(f"Maximum number of batches must be positive not {max_batches}! ")
        limited_num_batches = np.min((len(batches), max_batches))
        logger.info("Limiting number of output batches to %d", limited_num_batches)
        batches = batches[:limited_num_batches]
    
    # construct lists of object ids, indices and group_ids for each batch
    id_arrays = [obj_ids[inds] for inds in batches]
    index_arrays = batches
    group_arrays = [group_ids[inds] for inds in batches]

    return id_arrays, index_arrays, group_arrays


def write_id_list_product(ids, pointing_list, observation_id_list, tile_list, workdir, i, batch_uuid):
    """
    Writes a dpdSheObjectIdList product, returning its filename

    Inputs:
      - ids: the list of ids to write to the product
      - pointing_list: The list of pointing ids
      - observation_id: the observation id
      - tile_list: the list of tile_ids
      - workdir: the workdir
      - i: the instance id of this file

    Returns:
      - partially_qualified_id_list_filename: The output filename relative to the workdir
    """

    # Create the product and fill out metadata
    dpd = she_object_id_list.create_dpd_she_object_id_list(id_list = list(ids))

    dpd.Data.BatchIndex = i
    dpd.Data.PointingIdList = pointing_list
    dpd.Data.ObservationIdList = observation_id_list
    dpd.Data.TileList = tile_list

    # Get a filename for the product
    id_list_filename = FileNameProvider().get_allowed_filename(
        type_name="P-OBJ-ID-LIST",
        instance_id=batch_uuid,
        extension=".xml",
        release=SHE_CTE_RELEASE_STRING,
        processing_function="SHE",
    )

    fully_qualified_id_list_filename = os.path.join(workdir, id_list_filename)

    # Save the product
    save_product_metadata(dpd, fully_qualified_id_list_filename)

    return id_list_filename



def write_mer_product(input_mer_cat, inds, groups, mfc_dpd, workdir, batch_uuid):
    """
    Writes a listfile pointing to a single mer catalog product for a batch

    Inputs:
      - input_mer_cat - the catalogue we will be taking rows from
      - inds: The row indices to create the batched catalogue from
      - groups: the group ids of the objects in the batch
      - mfc_dpd: a dpdMerFinalCatalog object (from the inputs) that we use as a template for creating the new catalogue
      - workdir: the working directory
      - i: the instance id of this catalogue

    Returns:
      - partially_qualified_listfile_filename: the output filename relative to workdir

    """

    # get the data from the table for this batch, creating a new table
    # (The filled() method creates a copy with this data only)
    batch_table = input_mer_cat[inds].filled()
    
    #IF the group id is not already in the table, add it
    if mer_tf.GROUP_ID not in batch_table.columns:
        batch_table.add_column(groups, name = mer_tf.GROUP_ID)

    # Get a filename for the FITS table
    table_filename = FileNameProvider().get_allowed_filename(
        type_name="T-BATCH-MER-CAT",
        instance_id=batch_uuid,
        extension=".fits",
        release=SHE_CTE_RELEASE_STRING,
        processing_function="SHE",
    )

    fully_qualified_table_filename = os.path.join(workdir, "data", table_filename)

    batch_table.write(fully_qualified_table_filename)

    # set up the product
    dpd = deepcopy(mfc_dpd)

    dpd.Data.DataStorage.DataContainer.FileName = table_filename


    # Get a filename for the product
    mer_filename = FileNameProvider().get_allowed_filename(
        type_name="P-BATCH-MER-CAT",
        instance_id=batch_uuid,
        extension=".xml",
        release=SHE_CTE_RELEASE_STRING,
        processing_function="SHE",
    )

    fully_qualified_mer_filename = os.path.join(workdir, mer_filename)

    save_product_metadata(dpd, fully_qualified_mer_filename)

    return mer_filename

    

def object_id_split_from_args(args,
                              sub_batch = False):
    """ Core function for implementing a split by object ID
    """

    logger.debug('# Entering object_id_split_from_args(args)')

    workdir = args.workdir
    data_images = args.data_images
    mer_final_catalog_tables = args.mer_final_catalog_tables
    output_object_ids = args.object_ids
    output_batch_mer_catalogs = args.batch_mer_catalogs

    grouping_radius = args.grouping_radius
    ids_to_use = args.ids_to_use
    max_batches = args.max_batches
    batch_size = args.batch_size

    # Read the data we need from the VIS products
    wcs_list, pointing_ids, observation_id = read_vis_frames(data_images, workdir)
    
    # Read the mer final catalogues in, stacking them to one large catalogue
    mer_final_catalog, tile_ids, mfc_dpd = read_mer_catalogs(mer_final_catalog_tables, workdir)

    # Remove objects not detected in VIS
    if not args.include_vis_non_detections:
        n_before = len(mer_final_catalog)
        vis_det = mer_final_catalog["VIS_DET"] == 1
        mer_final_catalog = mer_final_catalog[vis_det]
        n_after = len(mer_final_catalog)
        logger.info("Removed %s non-VIS-detected objects from catalogue", n_before-n_after)

    # TODO: When moving to a Tile-Based pipeline, check that obs_ids from VIS exposures matches obs_ids in MER catalogue dpd
    
    # If we're sub batching, this has already been done in the batching stage
    if not sub_batch:
        # Remove objects outwith the FOV
        mer_final_catalog = prune_objects_not_in_FOV(mer_final_catalog, wcs_list, ids_to_use)
    
    id_arrays, index_arrays, group_arrays = group_and_batch_objects(mer_final_catalog, batch_size, max_batches, grouping_radius)
    
    # Loop over all the batches, writing the products

    id_prods = []
    mer_prods = []
    i = 0
    n_batches = len(id_arrays)
    for ids, inds, groups in zip(id_arrays, index_arrays, group_arrays):
        
        # generate a uuid for this batch (from the list of object_ids in the batch) to
        # use in the filenames of the object_id_list and mer_final_catalog products/files/
        # This greatly reduces the chances of filename conflitcs (which can happen!)
        batch_uuid = hash_any(ids, format="hex", max_length=10)

        id_prod = write_id_list_product(ids, pointing_ids, observation_id, tile_ids, workdir, i, batch_uuid)

        mer_prod = write_mer_product(mer_final_catalog, inds, groups, mfc_dpd, workdir, batch_uuid)

        id_prods.append(id_prod)
        mer_prods.append(mer_prod)

        logger.info("Writen batch %d of %d to file",i+1, n_batches)

        i+=1

    # Finally write the output listfiles

    qualified_mer_listfile = os.path.join(workdir,output_batch_mer_catalogs)
    with open(qualified_mer_listfile,"w") as f:
        json.dump(mer_prods, f)
    logger.info("Written output batch_mer_catalogs listfile to %s", qualified_mer_listfile)

    qualified_obj_ids_listfile = os.path.join(workdir,output_object_ids)
    with open(qualified_obj_ids_listfile,"w") as f:
        json.dump(id_prods, f)
    logger.info("Written output object_ids listfile to %s", qualified_obj_ids_listfile)
