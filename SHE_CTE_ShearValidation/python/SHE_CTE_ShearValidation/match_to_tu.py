""" @file match_to_tu.py

    Created 10 May 2019

    Code to implement matching of shear estimates catalogs to SIM's TU galaxy and star catalogs.
"""

__updated__ = "2020-07-02"

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
from SHE_PPT.logging import getLogger
from SHE_PPT.table_utility import table_to_hdu
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column, join, vstack, unique

import SHE_CTE
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
import numpy as np

logger = getLogger(__name__)

methods = ("BFD", "KSB", "LensMC", "MomentsML", "REGAUSS")

star_index_colname = "STAR_INDEX"
gal_index_colname = "GAL_INDEX"

max_coverage = 1.0  # deg


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

    return merged_catalog


def match_to_tu_from_args(args):
    """ @TODO Fill in docstring
    """

    # Read in the data products for SIM's TU galaxy and star catalogs, and get the filenames of the fits files from
    # them
    qualified_star_catalog_product_filename = file_io.find_file(args.tu_star_catalog,
                                                                path=args.workdir + ":" + args.sim_path)
    logger.info("Reading in True Universe star catalog product from " + qualified_star_catalog_product_filename)
    star_catalog_product = file_io.read_xml_product(qualified_star_catalog_product_filename)
    star_catalog_filenames = star_catalog_product.get_data_filenames()

    qualified_galaxy_catalog_product_filename = file_io.find_file(args.tu_galaxy_catalog,
                                                                  path=args.workdir + ":" + args.sim_path)
    logger.info("Reading in True Universe galaxy catalog product from " + qualified_galaxy_catalog_product_filename)
    galaxy_catalog_product = file_io.read_xml_product(qualified_galaxy_catalog_product_filename)
    galaxy_catalog_filenames = galaxy_catalog_product.get_data_filenames()

    # Read in the shear estimates data product, and get the filenames of the tables for each method from it.
    qualified_shear_estimates_product_filename = file_io.find_file(args.she_measurements_product,
                                                                   path=args.workdir)
    logger.info("Reading in Shear Estimates product from " + qualified_shear_estimates_product_filename)
    shear_estimates_product = file_io.read_xml_product(qualified_shear_estimates_product_filename)

    shear_tables = {}

    for method in methods:
        fn = shear_estimates_product.get_method_filename(method)
        if fn is None or fn == "None":
            shear_tables[method] = None
            logger.warning("No filename for method " + method + ".")
        else:
            qualified_filename = os.path.join(args.workdir, fn)
            logger.debug("Reading in shear estimates table from " + qualified_filename)
            shear_tables[method] = Table.read(os.path.join(args.workdir, fn))

    # Determine the ra/dec range covered by the shear estimates file

    ra_range = np.array((1e99, -1e99))
    dec_range = np.array((1e99, -1e99))

    for method in methods:

        # BFD has a bug in determining the range, so skip it here unless it's the only method
        if method == "BFD" and len(methods) > 1:
            continue

        shear_table = shear_tables[method]
        if shear_table is None:
            continue

        ra_col = shear_table[sm_tf.x_world]
        dec_col = shear_table[sm_tf.y_world]

        flags_col = shear_table[sm_tf.flags]

        good_ra_data = ra_col[flags_col == 0]
        good_dec_data = dec_col[flags_col == 0]

        if len(good_ra_data) == 0:
            continue

        # Check if the range in this method's table sets a new min/max for ra and dec
        ra_range[0] = np.min((ra_range[0], np.min(good_ra_data)))
        ra_range[1] = np.max((ra_range[1], np.max(good_ra_data)))

        dec_range[0] = np.min((dec_range[0], np.min(good_dec_data)))
        dec_range[1] = np.max((dec_range[1], np.max(good_dec_data)))

    if ra_range[1] < ra_range[0] or dec_range[1] < dec_range[0]:
        raise ValueError("Invalid range")

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
                                                                    path=args.sim_path)

            logger.info("Found " + str(len(overlapping_star_catalog)) + " stars in overlapping region.")

            overlapping_galaxy_catalog = select_true_universe_sources(catalog_filenames=galaxy_catalog_filenames,
                                                                      ra_range=local_ra_range,
                                                                      dec_range=local_dec_range,
                                                                      path=args.sim_path)

            logger.info("Found " + str(len(overlapping_galaxy_catalog)) + " galaxies in overlapping region.")

            # Remove unused columns in the star table

            overlapping_star_catalog.remove_columns(['DIST', 'TU_MAG_H_2MASS', 'SED_TEMPLATE', 'AV', 'TU_FNU_VIS',
                                                     'TU_FNU_Y_NISP', 'TU_FNU_J_NISP', 'TU_FNU_H_NISP', 'TU_FNU_G_DECAM',
                                                     'TU_FNU_R_DECAM', 'TU_FNU_I_DECAM', 'TU_FNU_Z_DECAM',
                                                     'TU_FNU_U_MEGACAM', 'TU_FNU_R_MEGACAM', 'TU_FNU_G_JPCAM',
                                                     'TU_FNU_I_PANSTARRS', 'TU_FNU_Z_PANSTARRS', 'TU_FNU_Z_HSC',
                                                     'TU_FNU_G_GAIA', 'TU_FNU_BP_GAIA', 'TU_FNU_RP_GAIA', 'TU_FNU_U_LSST',
                                                     'TU_FNU_G_LSST', 'TU_FNU_R_LSST', 'TU_FNU_I_LSST', 'TU_FNU_Z_LSST',
                                                     'TU_FNU_Y_LSST', 'TU_FNU_U_KIDS', 'TU_FNU_G_KIDS', 'TU_FNU_R_KIDS',
                                                     'TU_FNU_I_KIDS', 'TU_FNU_J_2MASS', 'TU_FNU_H_2MASS',
                                                     'TU_FNU_KS_2MASS', ])

            # Remove unused columns in the galaxy table

            overlapping_galaxy_catalog.remove_columns(['SOURCE_ID', 'RA', 'DEC', 'TU_MAG_R01_SDSS_ABS', 'TU_MAG_R01_SDSS',
                                                       'EXT_LAW', 'EBV', 'HALPHA_LOGFLAM_EXT', 'HBETA_LOGFLAM_EXT',
                                                       'O2_LOGFLAM_EXT', 'O3_LOGFLAM_EXT', 'N2_LOGFLAM_EXT',
                                                       'S2_LOGFLAM_EXT', 'AV', 'TU_FNU_VIS', 'TU_FNU_Y_NISP',
                                                       'TU_FNU_J_NISP', 'TU_FNU_H_NISP', 'TU_FNU_G_DECAM',
                                                       'TU_FNU_R_DECAM', 'TU_FNU_I_DECAM', 'TU_FNU_Z_DECAM',
                                                       'TU_FNU_U_MEGACAM', 'TU_FNU_R_MEGACAM', 'TU_FNU_G_JPCAM',
                                                       'TU_FNU_I_PANSTARRS', 'TU_FNU_Z_PANSTARRS', 'TU_FNU_Z_HSC',
                                                       'TU_FNU_G_GAIA', 'TU_FNU_BP_GAIA', 'TU_FNU_RP_GAIA',
                                                       'TU_FNU_U_LSST', 'TU_FNU_G_LSST', 'TU_FNU_R_LSST',
                                                       'TU_FNU_I_LSST', 'TU_FNU_Z_LSST', 'TU_FNU_Y_LSST',
                                                       'TU_FNU_U_KIDS', 'TU_FNU_G_KIDS', 'TU_FNU_R_KIDS',
                                                       'TU_FNU_I_KIDS', 'TU_FNU_J_2MASS', 'TU_FNU_H_2MASS',
                                                       'TU_FNU_KS_2MASS',
                                                       ])
            # Set up star and galaxy tables for matching

            ra_star = overlapping_star_catalog["RA"]
            dec_star = overlapping_star_catalog["DEC"]
            sky_coord_star = SkyCoord(ra=ra_star * units.degree, dec=dec_star * units.degree)

            overlapping_star_catalog.add_column(Column(np.arange(len(ra_star)), name=star_index_colname))

            ra_gal = overlapping_galaxy_catalog["RA_MAG"]
            dec_gal = overlapping_galaxy_catalog["DEC_MAG"]
            sky_coord_gal = SkyCoord(ra=ra_gal * units.degree, dec=dec_gal * units.degree)

            overlapping_galaxy_catalog.add_column(Column(np.arange(len(ra_gal)), name=gal_index_colname))

            # Perform match to SIM's tables for each method

            for method in methods:

                shear_table = shear_tables[method]
                if shear_table is None or len(shear_table) == 0:
                    continue

                logger.info("Performing match for method " + method + ".")

                ra_se = shear_table[sm_tf.x_world]
                dec_se = shear_table[sm_tf.y_world]
                sky_coord_se = SkyCoord(ra=ra_se * units.degree, dec=dec_se * units.degree)

                # Match to both star and galaxy tables, and determine which is best match
                best_star_id, best_star_distance, _ = sky_coord_se.match_to_catalog_sky(sky_coord_star)
                best_gal_id, best_gal_distance, _ = sky_coord_se.match_to_catalog_sky(sky_coord_gal)

                # Perform the reverse match as well, and only use a symmetric best-match table
                best_obj_id_from_star, _best_distance_from_star, _ = sky_coord_star.match_to_catalog_sky(sky_coord_se)
                best_obj_id_from_gal, _best_distance_from_gal, _ = sky_coord_gal.match_to_catalog_sky(sky_coord_se)

                # Check that the overall best distance is less than the threshold
                best_distance = np.where(best_gal_distance <= best_star_distance,
                                         best_gal_distance, best_star_distance)

                # Mask rows where the match isn't close enough, or to the other type of object, with -99

                in_ra_range = np.logical_and(
                    shear_table[sm_tf.x_world] >= local_ra_range[0], shear_table[sm_tf.x_world] < local_ra_range[1])
                in_dec_range = np.logical_and(
                    shear_table[sm_tf.y_world] >= local_dec_range[0], shear_table[sm_tf.y_world] < local_dec_range[1])

                in_range = np.logical_and(in_ra_range, in_dec_range)

                # Mask out with -99 if the distance is outside the threshold or it better matches to the
                # other type of object
                best_star_id = np.where(in_range, np.where(best_distance < args.match_threshold,
                                                           np.where(best_star_distance <
                                                                    best_gal_distance, best_star_id, -99),
                                                           -99), -99)
                best_gal_id = np.where(in_range, np.where(best_distance < args.match_threshold,
                                                          np.where(best_gal_distance <=
                                                                   best_star_distance, best_gal_id, -99),
                                                          -99), -99)

                # Check for symmetric matches
                symmetric_star_match = best_obj_id_from_star[best_star_id] == np.arange(len(best_star_id))
                symmetric_gal_match = best_obj_id_from_gal[best_gal_id] == np.arange(len(best_gal_id))

                # Note - this will give junk results for those that are already flagged -99, but since we
                # do an or-wise update in the next step, it isn't an issue here. Be warned though if
                # do anything else with these new symmetric_*_match arrays

                # Mask out with -99 if we don't have a symmetric match
                best_star_id = np.where(symmetric_star_match, best_star_id, -99)
                best_gal_id = np.where(symmetric_gal_match, best_gal_id, -99)

                # Add columns to the shear estimates table so we can match to it
                if star_index_colname in shear_table.colnames:
                    shear_table[star_index_colname] = best_star_id
                else:
                    shear_table.add_column(Column(best_star_id, name=star_index_colname))
                if gal_index_colname in shear_table.colnames:
                    shear_table[gal_index_colname] = best_gal_id
                else:
                    shear_table.add_column(Column(best_gal_id, name=gal_index_colname))

                # Match to the star and galaxy tables

                star_matched_table = join(shear_table, overlapping_star_catalog, keys=star_index_colname)
                logger.info("Matched " + str(len(star_matched_table)) + " objects to stars.")

                gal_matched_table = join(shear_table, overlapping_galaxy_catalog, keys=gal_index_colname)
                logger.info("Matched " + str(len(gal_matched_table)) + " objects to galaxies.")

                # Remove matched rows from the shear table
                matched_rows = np.logical_or(best_star_id > 0, best_gal_id > 0)
                matched_indices = (matched_rows * np.arange(len(matched_rows)))[matched_rows]
                shear_table.remove_rows(matched_indices)

                # Remove extra columns we no longer need
                star_matched_table.remove_columns([star_index_colname, gal_index_colname])
                gal_matched_table.remove_columns([star_index_colname, gal_index_colname])

                # Add extra useful columns to the galaxy-matched table for analysis

                # Details about estimated shear

                if not method == "BFD":

                    gal_matched_table.add_column(
                        Column(np.arctan2(gal_matched_table["G2"].data, gal_matched_table["G1"].data) * 90 / np.pi,
                               name="Beta_Est_Shear"))

                    g_mag = np.sqrt(gal_matched_table["G1"].data ** 2 + gal_matched_table["G2"].data ** 2)
                    gal_matched_table.add_column(Column(g_mag, name="Mag_Est_Shear"))

                    gal_matched_table.add_column(Column((1 - g_mag) / (1 + g_mag), name="Axis_Ratio_Est_Shear"))

                # Details about the input shear

                g1_in = gal_matched_table["GAMMA1"]
                g2_in = gal_matched_table["GAMMA2"]

                gal_matched_table.add_column(
                    Column(np.arctan2(g2_in, g1_in) * 90 / np.pi, name="Beta_Input_Shear"))

                gal_matched_table.add_column(
                    Column(np.sqrt(g1_in ** 2 + g2_in ** 2) * 90 / np.pi, name="Mag_Input_Shear"))

                # Details about the input bulge shape

                bulge_angle = gal_matched_table["DISK_ANGLE"] + 90
                regularized_bulge_angle = np.where(bulge_angle < -90, bulge_angle + 180,
                                                   np.where(bulge_angle > 90, bulge_angle - 180, bulge_angle))
                gal_matched_table.add_column(Column(regularized_bulge_angle,
                                                    name="Beta_Input_Bulge_Unsheared_Shape"))

                bulge_axis_ratio = gal_matched_table["DISK_AXIS_RATIO"]
                bulge_g_mag = (1 - bulge_axis_ratio) / (1 + bulge_axis_ratio)
                gal_matched_table.add_column(Column(bulge_g_mag, name="Mag_Input_Bulge_Unsheared_Shape"))

                gal_matched_table.add_column(Column(bulge_g_mag * np.cos(bulge_angle * np.pi / 90),
                                                    name="G1_Input_Bulge_Unsheared_Shape"))
                gal_matched_table.add_column(Column(bulge_g_mag * np.sin(bulge_angle * np.pi / 90),
                                                    name="G2_Input_Bulge_Unsheared_Shape"))

                # Details about the input disk shape

                disk_angle = gal_matched_table["DISK_ANGLE"] + 90
                regularized_disk_angle = np.where(disk_angle < -90, disk_angle + 180,
                                                  np.where(disk_angle > 90, disk_angle - 180, disk_angle))
                gal_matched_table.add_column(Column(regularized_disk_angle,
                                                    name="Beta_Input_Disk_Unsheared_Shape"))

                disk_axis_ratio = gal_matched_table["DISK_AXIS_RATIO"]
                disk_g_mag = (1 - disk_axis_ratio) / (1 + disk_axis_ratio)
                gal_matched_table.add_column(Column(disk_g_mag, name="Mag_Input_Disk_Unsheared_Shape"))

                gal_matched_table.add_column(Column(disk_g_mag * np.cos(disk_angle * np.pi / 90),
                                                    name="G1_Input_Disk_Unsheared_Shape"))
                gal_matched_table.add_column(Column(disk_g_mag * np.sin(disk_angle * np.pi / 90),
                                                    name="G2_Input_Disk_Unsheared_Shape"))

                # Add these tables to the dictionaries of tables
                star_matched_tables[method].append(star_matched_table)
                gal_matched_tables[method].append(gal_matched_table)

    # Create output data product
    matched_catalog_product = products.she_measurements.create_dpd_she_measurements()
    for method in methods:

        if len(gal_matched_tables[method]) == 0:
            matched_catalog_product.set_method_filename(method, "None")
            continue

        gal_matched_table = vstack(gal_matched_tables[method])
        star_matched_table = vstack(star_matched_tables[method])
        unmatched_table = shear_tables[method]
        unmatched_table.remove_columns([star_index_colname, gal_index_colname])

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
        hdulist.writeto(os.path.join(args.workdir, method_filename), clobber=True)

    # Write the data product
    logger.info("Writing output matched catalog data product to " + os.path.join(args.workdir,
                                                                                 args.matched_catalog))
    file_io.write_xml_product(matched_catalog_product, args.matched_catalog,
                              workdir=args.workdir)

    return
