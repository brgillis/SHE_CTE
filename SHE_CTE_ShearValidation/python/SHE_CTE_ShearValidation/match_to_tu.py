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
from astropy.table import Table
import numpy as np

logger = getLogger(__name__)

methods = ("BFD", "KSB", "LensMC", "MomentsML", "REGAUSS")


def select_true_universe_sources(catalog_file_names, ra_range, dec_range, path):
    """ Loads all the True Universe catalog files and selects only those
    sources that fall inside the specified (RA, Dec) region.

    """
    # Loop over the True Universe catalog files and select the relevant sources
    merged_catalog = None

    for file_name in catalog_file_names:
        # Load the catalog table
        catalog = Table.read(file_io.find_file(file_name, path=path), format="fits")

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
    qualified_galaxy_catalog_product_filename = file_io.find_file(args.tu_galaxy_catalog,
                                                                  path=args.workdir + ":" + args.sim_path)

    star_catalog_product = file_io.read_xml_product(qualified_star_catalog_product_filename)
    galaxy_catalog_product = file_io.read_xml_product(qualified_galaxy_catalog_product_filename)

    star_catalog_filenames = star_catalog_product.get_data_filenames()
    galaxy_catalog_filenames = galaxy_catalog_product.get_data_filenames()

    # Read in the shear estimates data product, and get the filenames of the tables for each method from it.
    shear_estimates_product = file_io.read_xml_product(file_io.find_file(args.shear_estimates_product,
                                                                         path=args.workdir))
    shear_tables = {}

    for method in methods:
        fn = shear_estimates_product.get_method_filename(method)
        if fn is None:
            shear_tables[method] = None
            logger.warn("No filename for method " + method + ".")
        else:
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

        ra_range[0] = np.min((ra_range[0], np.min(ra_col.data)))
        ra_range[1] = np.max((ra_range[1], np.max(ra_col.data)))

        dec_range[0] = np.min((dec_range[0], np.min(dec_col.data)))
        dec_range[1] = np.max((dec_range[1], np.max(dec_col.data)))

    if ra_range[1] < ra_range[0] or dec_range[1] < dec_range[0]:
        raise ValueError("Invalid range")

    # @TODO Perform match to SIM's tables

    # Create output data product
    matched_catalog_product = products.shear_estimates.create_shear_estimates_product()
    for method in methods:
        method_filename = file_io.get_allowed_filename("SHEAR-SIM-MATCHED-CAT", instance_id=str(os.getpid()),
                                                       extension=".fits", version=SHE_CTE.__version__, subdir="data",)

        logger.debug("Writing output matched catalog for method " + method + " to " + os.path.join(args.workdir,
                                                                                                   method_filename))
        # @TODO output matched table to args.workdir/method_filename

    # Write the data product
    logger.info("Writing output matched catalog data product to " + os.path.join(args.workdir,
                                                                                 args.matched_catalog))
    file_io.write_xml_product(product=matched_catalog_product, filename=os.path.join(args.workdir,
                                                                                     args.matched_catalog))

    return
