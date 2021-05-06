""" @file match_to_tu.py

    Created 10 May 2019

    Code to implement matching of shear estimates catalogs to SIM's TU galaxy and star catalogs.
"""

__updated__ = "2021-05-06"

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

import os

from SHE_PPT import file_io
from SHE_PPT import products
from SHE_PPT.file_io import read_listfile
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.she_bfd_moments import tf as bfdm_tf
from SHE_PPT.table_formats.she_ksb_measurements import tf as ksbm_tf
from SHE_PPT.table_formats.she_lensmc_measurements import tf as lmcm_tf
from SHE_PPT.table_formats.she_momentsml_measurements import tf as mmlm_tf
from SHE_PPT.table_formats.she_regauss_measurements import tf as regm_tf
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import table_to_hdu
from astropy.table import Table, Column, join, vstack

import SHE_CTE
import numpy as np

logger = getLogger(__name__)

methods = ("BFD", "KSB", "LensMC", "MomentsML", "REGAUSS")

star_index_colname = "STAR_INDEX"
gal_index_colname = "GAL_INDEX"

max_coverage = 1.0  # deg

shear_estimation_method_table_formats = {"KSB": ksbm_tf,
                                         "REGAUSS": regm_tf,
                                         "MomentsML": mmlm_tf,
                                         "LensMC": lmcm_tf,
                                         "BFD": bfdm_tf}

galcat_g1_colname = "GAMMA1"
galcat_g2_colname = "GAMMA2"
galcat_kappa_colname = "KAPPA"
galcat_bulge_angle_colname = "DISK_ANGLE"
galcat_bulge_axis_ratio_colname = "DISK_AXIS_RATIO"
galcat_disk_angle_colname = "DISK_ANGLE"
galcat_disk_axis_ratio_colname = "DISK_AXIS_RATIO"


def select_true_universe_sources(catalog_filenames, ra_range, dec_range, path):
    """ Loads all the True Universe catalog files and selects only those
    sources that fall inside the specified (RA, Dec) region.

    """
    # Loop over the True Universe catalog files and select the relevant sources
    merged_catalog = None

    logger.info("Reading in overlapping sources.")

    for filename in catalog_filenames:

        qualified_filename = file_io.find_file(filename, path=path)

        logger.debug("Reading overlapping sources from " + qualified_filename + ".")

        # Load the catalog table

        try:
            catalog = Table.read(qualified_filename, format="fits")
        except OSError as e:
            logger.error(filename + " is corrupt or missing")
            raise
        except Exception as e:
            logger.error("Error reading in catalog " + filename)
            raise

        # Get the (RA, Dec) columns
        ra = catalog["RA_MAG"] if "RA_MAG" in catalog.colnames else catalog["RA"]
        dec = catalog["DEC_MAG"] if "DEC_MAG" in catalog.colnames else catalog["DEC"]

        # Check which sources fall inside the given (RA, Dec) ranges
        cond_ra = np.logical_and(ra > ra_range[0], ra < ra_range[1])
        cond_dec = np.logical_and(dec > dec_range[0], dec < dec_range[1])
        cond = np.logical_and(cond_ra, cond_dec)

        if np.any(cond):
            # Add the selected sources to the merged catalog
            if merged_catalog is None:
                merged_catalog = catalog[cond]
            else:
                merged_catalog = vstack([merged_catalog, catalog[cond]])

    if merged_catalog is None:
        logger.warning("No TU sources found in region: \n" +
                       f"R.A.: {ra_range}\n" +
                       f"Dec.: {dec_range}")
        merged_catalog = catalog[cond]

    return merged_catalog


def match_to_tu_from_args(args):
    """ @TODO Fill in docstring
    """

    # Get the lists of star and galaxy filenames

    search_path = args.workdir + ":" + args.sim_path

    if args.tu_star_catalog_list is not None:
        qualified_tu_star_catalog_list = file_io.find_file(args.tu_star_catalog_list,
                                                           path=search_path)
        tu_star_catalog_product_filenames = read_listfile(qualified_tu_star_catalog_list)
    elif args.tu_star_catalog is not None:
        tu_star_catalog_product_filenames = [args.tu_star_catalog]
    else:
        raise ValueError("No star catalogs provided to match with.")

    if args.tu_galaxy_catalog_list is not None:
        qualified_tu_galaxy_catalog_list = file_io.find_file(args.tu_galaxy_catalog_list,
                                                             path=search_path)
        tu_galaxy_catalog_product_filenames = read_listfile(qualified_tu_galaxy_catalog_list)
    elif args.tu_galaxy_catalog is not None:
        tu_galaxy_catalog_product_filenames = [args.tu_galaxy_catalog]
    else:
        raise ValueError("No galaxy catalogs provided to match with.")

    # Read in the data products for SIM's TU galaxy and star catalogs, and get the filenames of the fits files from
    # them

    star_catalog_filenames = []
    for tu_star_catalog_product_filename in tu_star_catalog_product_filenames:
        qualified_star_catalog_product_filename = file_io.find_file(tu_star_catalog_product_filename,
                                                                    path=search_path)
        logger.info("Reading in True Universe star catalog product from " + qualified_star_catalog_product_filename)
        star_catalog_product = file_io.read_xml_product(qualified_star_catalog_product_filename)
        star_catalog_filenames.append(star_catalog_product.get_data_filename())

    galaxy_catalog_filenames = []
    for tu_galaxy_catalog_product_filename in tu_galaxy_catalog_product_filenames:
        qualified_galaxy_catalog_product_filename = file_io.find_file(tu_galaxy_catalog_product_filename,
                                                                      path=search_path)
        logger.info("Reading in True Universe galaxy catalog product from " + qualified_galaxy_catalog_product_filename)
        galaxy_catalog_product = file_io.read_xml_product(qualified_galaxy_catalog_product_filename)
        galaxy_catalog_filenames.append(galaxy_catalog_product.get_data_filename())

    # Read in the shear estimates data product, and get the filenames of the tables for each method from it.
    qualified_shear_estimates_product_filename = file_io.find_file(args.she_measurements_product,
                                                                   path=args.workdir)
    logger.info("Reading in Shear Estimates product from " + qualified_shear_estimates_product_filename)
    shear_estimates_product = file_io.read_xml_product(qualified_shear_estimates_product_filename)

    shear_tables = {}

    for method in methods:
        fn = shear_estimates_product.get_method_filename(method)
        if fn is None or fn == "None" or fn == "data/None":
            shear_tables[method] = None
            logger.warning("No filename for method " + method + ".")
        else:
            qualified_filename = os.path.join(args.workdir, fn)
            logger.debug("Reading in shear estimates table from " + qualified_filename)
            shear_tables[method] = Table.read(os.path.join(args.workdir, fn))

    # Determine the ra/dec range covered by the shear estimates file

    good_ra_range = np.array((1e99, -1e99))
    good_dec_range = np.array((1e99, -1e99))
    full_ra_range = np.array((1e99, -1e99))
    full_dec_range = np.array((1e99, -1e99))

    for method in methods:

        # BFD has a bug in determining the range, so skip it here unless it's the only method
        if method == "BFD" and len(methods) > 1:
            continue

        shear_table = shear_tables[method]
        if shear_table is None:
            continue

        sem_tf = shear_estimation_method_table_formats[method]

        ra_col = shear_table[sem_tf.ra]
        dec_col = shear_table[sem_tf.dec]

        flags_col = shear_table[sem_tf.fit_flags]

        good_ra_data = ra_col[flags_col == 0]
        good_dec_data = dec_col[flags_col == 0]

        # Check if the range in this method's table sets a new min/max for ra and dec

        if len(good_ra_data) > 0:
            good_ra_range[0] = np.min((good_ra_range[0], np.min(good_ra_data)))
            good_ra_range[1] = np.max((good_ra_range[1], np.max(good_ra_data)))
            good_dec_range[0] = np.min((good_dec_range[0], np.min(good_dec_data)))
            good_dec_range[1] = np.max((good_dec_range[1], np.max(good_dec_data)))

        full_ra_range[0] = np.min((full_ra_range[0], np.min(ra_col)))
        full_ra_range[1] = np.max((full_ra_range[1], np.max(ra_col)))
        full_dec_range[0] = np.min((full_dec_range[0], np.min(dec_col)))
        full_dec_range[1] = np.max((full_dec_range[1], np.max(dec_col)))

    if good_ra_range[1] < good_ra_range[0] or good_dec_range[1] < good_dec_range[0]:
        logger.warning("Invalid range or no valid data for any method.")
        if not (full_ra_range[1] > full_ra_range[0] and full_dec_range[1] > full_dec_range[0]):
            raise ValueError("No valid position data for any method.")
        else:
            ra_range = full_ra_range
            dec_range = full_dec_range
    else:
        ra_range = good_ra_range
        dec_range = good_dec_range

    # Pad the range by the threshold amount
    ra_range[0] -= args.match_threshold
    ra_range[1] += args.match_threshold
    dec_range[0] -= args.match_threshold
    dec_range[1] += args.match_threshold

    logger.info("Object range is: ")
    logger.info("  RA : " + str(ra_range[0]) + " to " + str(ra_range[1]))
    logger.info("  DEC: " + str(dec_range[0]) + " to " + str(dec_range[1]))

    ra_limits = np.linspace(ra_range[0], ra_range[1], num=int(
        (ra_range[1] - ra_range[0]) / max_coverage) + 2, endpoint=True)
    dec_limits = np.linspace(dec_range[0], dec_range[1], num=int(
        (dec_range[1] - dec_range[0]) / max_coverage) + 2, endpoint=True)

    star_matched_tables = {}
    gal_matched_tables = {}

    for method in methods:
        star_matched_tables[method] = []
        gal_matched_tables[method] = []

    for ra_i in range(len(ra_limits) - 1):
        for dec_i in range(len(dec_limits) - 1):
            local_ra_range = np.array((ra_limits[ra_i], ra_limits[ra_i + 1]))
            local_dec_range = np.array((dec_limits[dec_i], dec_limits[dec_i + 1]))

            logger.info("Processing batch (" + str(ra_i + 1) + "," + str(dec_i + 1) + ") of (" +
                        str(len(ra_limits) - 1) + "," + str(len(dec_limits) - 1) + ")")
            logger.info("Processing ra range: " + str(local_ra_range))
            logger.info("      and dec range: " + str(local_dec_range))

            # Read in the star and galaxy catalogs from the overlapping area
            overlapping_star_catalog = select_true_universe_sources(catalog_filenames=star_catalog_filenames,
                                                                    ra_range=local_ra_range,
                                                                    dec_range=local_dec_range,
                                                                    path=search_path)

            logger.info("Found " + str(len(overlapping_star_catalog)) + " stars in overlapping region.")

            overlapping_galaxy_catalog = select_true_universe_sources(catalog_filenames=galaxy_catalog_filenames,
                                                                      ra_range=local_ra_range,
                                                                      dec_range=local_dec_range,
                                                                      path=search_path)

            logger.info("Found " + str(len(overlapping_galaxy_catalog)) + " galaxies in overlapping region.")

            # Remove unused columns in the star table

            overlapping_star_catalog.remove_columns(['DIST', 'TU_MAG_H_2MASS', 'SED_TEMPLATE',
                                                     'AV', 'TU_FNU_VIS', 'TU_FNU_Y_NISP', 'TU_FNU_J_NISP',
                                                     'TU_FNU_H_NISP', 'TU_FNU_G_DECAM', 'TU_FNU_R_DECAM',
                                                     'TU_FNU_I_DECAM', 'TU_FNU_Z_DECAM', 'TU_FNU_U_MEGACAM',
                                                     'TU_FNU_R_MEGACAM', 'TU_FNU_G_JPCAM', 'TU_FNU_I_PANSTARRS',
                                                     'TU_FNU_Z_PANSTARRS', 'TU_FNU_Z_HSC', 'TU_FNU_G_GAIA',
                                                     'TU_FNU_BP_GAIA', 'TU_FNU_RP_GAIA', 'TU_FNU_U_LSST',
                                                     'TU_FNU_G_LSST', 'TU_FNU_R_LSST', 'TU_FNU_I_LSST',
                                                     'TU_FNU_Z_LSST', 'TU_FNU_Y_LSST', 'TU_FNU_U_KIDS',
                                                     'TU_FNU_G_KIDS', 'TU_FNU_R_KIDS', 'TU_FNU_I_KIDS',
                                                     'TU_FNU_J_2MASS', 'TU_FNU_H_2MASS', 'TU_FNU_KS_2MASS',
                                                     'GAIA_RA', 'GAIA_RA_ERROR', 'GAIA_DEC', 'GAIA_DEC_ERROR',
                                                     'GAIA_PHOT_G_MEAN_FLUX', 'GAIA_PHOT_G_MEAN_FLUX_ERROR',
                                                     'GAIA_PHOT_BP_MEAN_FLUX', 'GAIA_PHOT_BP_MEAN_FLUX_ERROR',
                                                     'GAIA_PHOT_RP_MEAN_FLUX', 'GAIA_PHOT_RP_MEAN_FLUX_ERROR'])

            # Remove unused columns in the galaxy table
            overlapping_galaxy_catalog.remove_columns(['SOURCE_ID', 'HALO_ID', 'RA', 'DEC',
                                                       'EXT_LAW', 'EBV', 'HALPHA_LOGFLAM_EXT_MAG',
                                                       'HBETA_LOGFLAM_EXT_MAG', 'O2_LOGFLAM_EXT_MAG',
                                                       'O3_LOGFLAM_EXT_MAG', 'N2_LOGFLAM_EXT_MAG',
                                                       'S2_LOGFLAM_EXT_MAG', 'AV', 'TU_FNU_VIS_MAG',
                                                       'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG',
                                                       'TU_FNU_G_DECAM_MAG', 'TU_FNU_R_DECAM_MAG',
                                                       'TU_FNU_I_DECAM_MAG', 'TU_FNU_Z_DECAM_MAG',
                                                       'TU_FNU_U_MEGACAM_MAG', 'TU_FNU_R_MEGACAM_MAG',
                                                       'TU_FNU_G_JPCAM_MAG', 'TU_FNU_I_PANSTARRS_MAG',
                                                       'TU_FNU_Z_PANSTARRS_MAG', 'TU_FNU_Z_HSC_MAG',
                                                       'TU_FNU_G_GAIA_MAG', 'TU_FNU_BP_GAIA_MAG', 'TU_FNU_RP_GAIA_MAG',
                                                       'TU_FNU_U_LSST_MAG', 'TU_FNU_G_LSST_MAG', 'TU_FNU_R_LSST_MAG',
                                                       'TU_FNU_I_LSST_MAG', 'TU_FNU_Z_LSST_MAG', 'TU_FNU_Y_LSST_MAG',
                                                       'TU_FNU_U_KIDS_MAG', 'TU_FNU_G_KIDS_MAG', 'TU_FNU_R_KIDS_MAG',
                                                       'TU_FNU_I_KIDS_MAG', 'TU_FNU_J_2MASS_MAG',
                                                       'TU_FNU_H_2MASS_MAG', 'TU_FNU_KS_2MASS_MAG'])

            # Set up star and galaxy tables for matching

            ra_star = overlapping_star_catalog["RA"]
            dec_star = overlapping_star_catalog["DEC"]
            sky_coord_star = SkyCoord(ra=ra_star, dec=dec_star)

            overlapping_star_catalog.add_column(Column(np.arange(len(ra_star)), name=star_index_colname))

            ra_gal = overlapping_galaxy_catalog["RA_MAG"]
            dec_gal = overlapping_galaxy_catalog["DEC_MAG"]
            sky_coord_gal = SkyCoord(ra=ra_gal, dec=dec_gal)

            overlapping_galaxy_catalog.add_column(Column(np.arange(len(ra_gal)), name=gal_index_colname))

            # Perform match to SIM's tables for each method

            for method in methods:

                unpruned_shear_table = shear_tables[method]
                if unpruned_shear_table is None:
                    logger.info("No catalog provided for method " + method + ".")
                    continue

                sem_tf = shear_estimation_method_table_formats[method]

                # Prune any rows with NaN for R.A. or Dec. from the shear table
                good_rows = ~np.logical_or(np.isnan(unpruned_shear_table[sem_tf.ra]),
                                           np.isnan(unpruned_shear_table[sem_tf.dec]))
                shear_table = unpruned_shear_table[good_rows]

                if len(shear_table) == 0:
                    logger.info("No valid rows in catalog for method " + method + ".")
                    continue

                ra_se = shear_table[sem_tf.ra]
                dec_se = shear_table[sem_tf.dec]
                sky_coord_se = SkyCoord(ra=ra_se * units.degree, dec=dec_se * units.degree)

                if len(sky_coord_star) > 0:

                    # Match to star table
                    best_star_id, best_star_distance, _ = sky_coord_se.match_to_catalog_sky(sky_coord_star)

                    # Perform the reverse match as well, and only use a symmetric best-match table
                    best_obj_id_from_star, _best_distance_from_star, _ = sky_coord_star.match_to_catalog_sky(
                        sky_coord_se)

                else:
                    best_star_id, best_star_distance = [], []
                    best_obj_id_from_star, _best_distance_from_star = [], []

                if len(sky_coord_gal) > 0:

                    # Match to galaxy table
                    best_gal_id, best_gal_distance, _ = sky_coord_se.match_to_catalog_sky(sky_coord_gal)

                    # Perform the reverse match as well, and only use a symmetric best-match table
                    best_obj_id_from_gal, _best_distance_from_gal, _ = sky_coord_gal.match_to_catalog_sky(sky_coord_se)

                else:
                    best_gal_id, best_gal_distance = [], []
                    best_obj_id_from_gal, _best_distance_from_gal = [], []

                # Check that the overall best distance is less than the threshold
                best_distance = np.where(best_gal_distance <= best_star_distance,
                                         best_gal_distance, best_star_distance)

                # Mask rows where the match isn't close enough, or to the other type of object, with -99

                in_ra_range = np.logical_and(
                    shear_table[sem_tf.ra] >= local_ra_range[0], shear_table[sem_tf.ra] < local_ra_range[1])
                in_dec_range = np.logical_and(
                    shear_table[sem_tf.dec] >= local_dec_range[0], shear_table[sem_tf.dec] < local_dec_range[1])

                in_range = np.logical_and(in_ra_range, in_dec_range)

                # Mask out with -99 if the distance is outside the threshold or it better matches to the
                # other type of object

                if len(sky_coord_star) > 0:
                    best_star_id = np.where(in_range, np.where(best_distance < args.match_threshold,
                                                               np.where(best_star_distance <
                                                                        best_gal_distance, best_star_id, -99),
                                                               -99), -99)
                    symmetric_star_match = best_obj_id_from_star[best_star_id] == np.arange(len(best_star_id))

                    # Mask out with -99 if we don't have a symmetric match
                    best_star_id = np.where(symmetric_star_match, best_star_id, -99)
                else:
                    best_star_id = np.zeros(len(in_range), dtype=int)

                if len(sky_coord_gal) > 0:
                    best_gal_id = np.where(in_range, np.where(best_distance < args.match_threshold,
                                                              np.where(best_gal_distance <=
                                                                       best_star_distance, best_gal_id, -99),
                                                              -99), -99)
                    symmetric_gal_match = best_obj_id_from_gal[best_gal_id] == np.arange(len(best_gal_id))

                    # Mask out with -99 if we don't have a symmetric match
                    best_gal_id = np.where(symmetric_gal_match, best_gal_id, -99)
                else:
                    best_gal_id = np.zeros(len(in_range), dtype=int)

                # Add columns to the shear estimates table so we can match to it
                if star_index_colname in shear_table.colnames:
                    if len(best_star_id) > 0:
                        shear_table[star_index_colname] = best_star_id
                else:
                    shear_table.add_column(Column(best_star_id, name=star_index_colname))
                if gal_index_colname in shear_table.colnames:
                    if len(best_gal_id) > 0:
                        shear_table[gal_index_colname] = best_gal_id
                else:
                    shear_table.add_column(Column(best_gal_id, name=gal_index_colname))

                # Match to the star and galaxy tables

                if len(sky_coord_star) > 0:
                    star_matched_table = join(shear_table, overlapping_star_catalog, keys=star_index_colname)
                    logger.info("Matched " + str(len(star_matched_table)) + " objects to stars.")
                else:
                    star_matched_table = shear_table[False * np.ones(len(shear_table), dtype=bool)]

                if len(sky_coord_gal) > 0:
                    gal_matched_table = join(shear_table, overlapping_galaxy_catalog, keys=gal_index_colname)
                    logger.info("Matched " + str(len(gal_matched_table)) + " objects to galaxies.")
                else:
                    gal_matched_table = shear_table[False * np.ones(len(shear_table), dtype=bool)]

                # Remove matched rows from the shear table
                matched_rows = np.logical_or(best_star_id > 0, best_gal_id > 0)
                matched_indices = (matched_rows * np.arange(len(matched_rows)))[matched_rows]
                shear_table.remove_rows(matched_indices)

                # Remove extra columns we no longer need
                star_matched_table.remove_columns([star_index_colname, gal_index_colname])
                gal_matched_table.remove_columns([star_index_colname, gal_index_colname])

                # Add these tables to the dictionaries of tables
                star_matched_tables[method].append(star_matched_table)
                gal_matched_tables[method].append(gal_matched_table)

                # Skip to next batch if we didn't match any galaxies
                if len(gal_matched_table) == 0:
                    continue

                # Add extra useful columns to the galaxy-matched table for analysis

                # Details about estimated shear

                if not method == "BFD":

                    gal_matched_table.add_column(
                        Column(np.arctan2(gal_matched_table[sem_tf.g2].data, gal_matched_table[sem_tf.g1].data) * 90 / np.pi,
                               name="Beta_Est_Shear"))

                    g_mag = np.sqrt(gal_matched_table[sem_tf.g1].data ** 2 + gal_matched_table[sem_tf.g2].data ** 2)
                    gal_matched_table.add_column(Column(g_mag, name="Mag_Est_Shear"))

                # Details about the input shear

                g1_in = gal_matched_table[galcat_g1_colname] / (1 - gal_matched_table[galcat_kappa_colname])
                g2_in = gal_matched_table[galcat_g2_colname] / (1 - gal_matched_table[galcat_kappa_colname])

                gal_matched_table.add_column(
                    Column(np.arctan2(g2_in, g1_in) * 90 / np.pi, name="Beta_Input_Shear"))

                gal_matched_table.add_column(
                    Column(np.sqrt(g1_in ** 2 + g2_in ** 2), name="Mag_Input_Shear"))

                # Details about the input bulge shape

                bulge_angle = gal_matched_table[galcat_bulge_angle_colname] + 90
                regularized_bulge_angle = np.where(bulge_angle < -90, bulge_angle + 180,
                                                   np.where(bulge_angle > 90, bulge_angle - 180, bulge_angle))
                gal_matched_table.add_column(Column(regularized_bulge_angle,
                                                    name="Beta_Input_Bulge_Unsheared_Shape"))

                # Details about the input disk shape

                disk_angle = gal_matched_table[galcat_disk_angle_colname] + 90
                regularized_disk_angle = np.where(disk_angle < -90, disk_angle + 180,
                                                  np.where(disk_angle > 90, disk_angle - 180, disk_angle))
                gal_matched_table.add_column(Column(regularized_disk_angle,
                                                    name="Beta_Input_Disk_Unsheared_Shape"))

    # Create output data product
    matched_catalog_product = products.she_measurements.create_dpd_she_measurements()
    for method in methods:

        if len(gal_matched_tables[method]) == 0:
            matched_catalog_product.set_method_filename(method, "None")
            continue

        gal_matched_table = vstack(gal_matched_tables[method])
        if len(gal_matched_table) == 0:
            logger.warn(f"No measurements with method {method} were matched to galaxies.")
        star_matched_table = vstack(star_matched_tables[method])
        if len(star_matched_table) == 0:
            logger.warn(f"No measurements with method {method} were matched to stars.")
        unmatched_table = shear_tables[method]

        method_filename = file_io.get_allowed_filename("SHEAR-SIM-MATCHED-CAT",
                                                       instance_id=method.upper() + "-" + str(os.getpid()),
                                                       extension=".fits", version=SHE_CTE.__version__, subdir="data",)
        matched_catalog_product.set_method_filename(method, method_filename)

        # Turn each table into an HDU and add it to an HDU list
        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU())

        # Add the galaxy table first, since it's more relevant
        hdulist.append(table_to_hdu(gal_matched_table))
        hdulist.append(table_to_hdu(star_matched_table))
        hdulist.append(table_to_hdu(unmatched_table))

        # Write out the HDU list to a file
        logger.info("Writing output matched catalogs for method " + method + " to " + os.path.join(args.workdir,
                                                                                                   method_filename))
        hdulist.writeto(os.path.join(args.workdir, method_filename), overwrite=True)

    # Write the data product
    logger.info("Writing output matched catalog data product to " + os.path.join(args.workdir,
                                                                                 args.matched_catalog))
    file_io.write_xml_product(matched_catalog_product, args.matched_catalog,
                              workdir=args.workdir)

    return
