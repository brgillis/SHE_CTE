""" @file match_to_tu.py

    Created 10 May 2019

    Code to implement matching of shear estimates catalogs to SIM's TU galaxy and star catalogs.
"""

__updated__ = "2019-05-10"

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

import SHE_CTE
from SHE_PPT import file_io
from SHE_PPT import products
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.shear_estimates import tf as setf
from SHE_PPT.table_utility import table_to_hdu
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column, join, vstack
import numpy as np


logger = getLogger(__name__)

methods = ("BFD", "KSB", "LensMC", "MomentsML", "REGAUSS")

star_index_colname = "STAR_INDEX"
gal_index_colname = "GAL_INDEX"


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
        catalog = Table.read(qualified_filename, format="fits")

        # Get the (RA, Dec) columns
        ra = catalog["ra_mag"] if "ra_mag" in catalog.colnames else catalog["RA"]
        dec = catalog["dec_mag"] if "dec_mag" in catalog.colnames else catalog["DEC"]

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
    qualified_shear_estimates_product_filename = file_io.find_file(args.shear_estimates_product,
                                                                   path=args.workdir)
    logger.info("Reading in Shear Estimates product from " + qualified_shear_estimates_product_filename)
    shear_estimates_product = file_io.read_xml_product(qualified_shear_estimates_product_filename)

    shear_tables = {}

    for method in methods:
        fn = shear_estimates_product.get_method_filename(method)
        if fn is None or fn == "None":
            shear_tables[method] = None
            logger.warn("No filename for method " + method + ".")
        else:
            qualified_filename = os.path.join(args.workdir, fn)
            logger.debug("Reading in shear estimates table from " + qualified_filename)
            shear_tables[method] = Table.read(os.path.join(args.workdir, fn))

    # Determine the ra/dec range covered by the shear estimates file

    ra_range = np.array((1e99, -1e99))
    dec_range = np.array((1e99, -1e99))

    for method in methods:

        shear_table = shear_tables[method]
        if shear_table is None:
            continue

        ra_col = shear_table[setf.x_world]
        dec_col = shear_table[setf.y_world]

        # Check if the range in this method's table sets a new min/max for ra and dec
        ra_range[0] = np.min((ra_range[0], np.min(ra_col.data)))
        ra_range[1] = np.max((ra_range[1], np.max(ra_col.data)))

        dec_range[0] = np.min((dec_range[0], np.min(dec_col.data)))
        dec_range[1] = np.max((dec_range[1], np.max(dec_col.data)))

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

    # Read in the star and galaxy catalogs from the overlapping area
    overlapping_star_catalog = select_true_universe_sources(catalog_filenames=star_catalog_filenames,
                                                            ra_range=ra_range,
                                                            dec_range=dec_range,
                                                            path=args.sim_path)

    logger.info("Found " + str(len(overlapping_star_catalog)) + " stars in overlapping region.")

    overlapping_galaxy_catalog = select_true_universe_sources(catalog_filenames=galaxy_catalog_filenames,
                                                              ra_range=ra_range,
                                                              dec_range=dec_range,
                                                              path=args.sim_path)

    logger.info("Found " + str(len(overlapping_galaxy_catalog)) + " galaxies in overlapping region.")

    # Set up star and galaxy tables for matching

    ra_star = overlapping_star_catalog["RA"]
    dec_star = overlapping_star_catalog["DEC"]
    sky_coord_star = SkyCoord(ra=ra_star * units.degree, dec=dec_star * units.degree)

    overlapping_star_catalog.add_column(Column(np.arange(len(ra_star)), name=star_index_colname))

    ra_gal = overlapping_galaxy_catalog["ra_mag"]
    dec_gal = overlapping_galaxy_catalog["dec_mag"]
    sky_coord_gal = SkyCoord(ra=ra_gal * units.degree, dec=dec_gal * units.degree)

    overlapping_galaxy_catalog.add_column(Column(np.arange(len(ra_gal)), name=gal_index_colname))

    # Perform match to SIM's tables for each method

    star_matched_tables = {}
    gal_matched_tables = {}

    for method in methods:

        shear_table = shear_tables[method]
        if shear_table is None:
            star_matched_tables[method] = None
            gal_matched_tables[method] = None
            continue

        logger.info("Performing match for method " + method + ".")

        ra_se = shear_table[setf.x_world]
        dec_se = shear_table[setf.y_world]
        sky_coord_se = SkyCoord(ra=ra_se * units.degree, dec=dec_se * units.degree)

        # Match to both star and galaxy tables, and determine which is best match
        best_star_id, best_star_distance, _ = sky_coord_se.match_to_catalog_sky(sky_coord_star)
        best_gal_id, best_gal_distance, _ = sky_coord_se.match_to_catalog_sky(sky_coord_gal)

        assert(len(best_star_id) == len(best_gal_id))

        # Check that the overall best distance is less than the threshold
        best_distance = np.where(best_gal_distance <= best_star_distance, best_gal_distance, best_star_distance)

        # Mask rows where the match isn't close enough, or to the other type of object, with -99
        best_star_id = np.where(best_distance < args.match_threshold,
                                np.where(best_star_distance < best_gal_distance, best_star_id, -99),
                                -99)
        best_gal_id = np.where(best_distance < args.match_threshold,
                               np.where(best_gal_distance <= best_star_distance, best_gal_id, -99),
                               -99)

        # Add columns to the shear estimates table so we can match to it
        shear_table.add_column(Column(best_star_id, name=star_index_colname))
        shear_table.add_column(Column(best_gal_id, name=gal_index_colname))

        # Match to the star and galaxy tables
        star_matched_tables[method] = join(shear_table, overlapping_star_catalog, keys=star_index_colname)
        gal_matched_tables[method] = join(shear_table, overlapping_galaxy_catalog, keys=gal_index_colname)

    # Create output data product
    matched_catalog_product = products.shear_estimates.create_shear_estimates_product()
    for method in methods:

        if gal_matched_tables[method] is None:
            matched_catalog_product.set_method_filename(method, "None")
            continue

        method_filename = file_io.get_allowed_filename("SHEAR-SIM-MATCHED-CAT", instance_id=str(os.getpid()),
                                                       extension=".fits", version=SHE_CTE.__version__, subdir="data",)
        matched_catalog_product.set_method_filename(method, method_filename)

        # Turn each table into an HDU and add it to an HDU list
        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU())

        # Add the galaxy table first, since it's more relevant
        hdulist.append(table_to_hdu(gal_matched_tables[method]))
        hdulist.append(table_to_hdu(star_matched_tables[method]))

        # Write out the HDU list to a file
        logger.debug("Writing output matched catalogs for method " + method + " to " + os.path.join(args.workdir,
                                                                                                    method_filename))
        hdulist.write(os.path.join(args.workdir, method_filename), clobber=True)

    # Write the data product
    logger.info("Writing output matched catalog data product to " + os.path.join(args.workdir,
                                                                                 args.matched_catalog))
    file_io.write_xml_product(product=matched_catalog_product, filename=os.path.join(args.workdir,
                                                                                     args.matched_catalog))

    return
