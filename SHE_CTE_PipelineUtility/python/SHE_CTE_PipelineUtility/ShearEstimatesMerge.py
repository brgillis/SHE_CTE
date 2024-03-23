""" @file ShearEstimatesMerge.py

    Created 14 Mar 2019

    Merge point executable for the split of shear estimation over object ID, merging into a single output product
    per Field of View.
"""

__updated__ = "2021-08-18"

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

import argparse
import multiprocessing as mp
import os
import itertools

from astropy import table

from ST_DM_DmUtils import DmUtils

from SHE_PPT import products
from SHE_PPT.argument_parser import dir_path
from SHE_PPT.constants.shear_estimation_methods import ShearEstimationMethods
from SHE_PPT.file_io import get_allowed_filename, read_listfile, read_xml_product, write_xml_product
from SHE_PPT.logging import getLogger
from SHE_PPT.table_formats.she_lensmc_chains import tf as lmcc_tf
from SHE_PPT.table_formats.she_measurements import tf as sm_tf
from SHE_PPT.table_utility import is_in_format

import SHE_CTE

METHOD_CATALOG_NAMES = {
    ShearEstimationMethods.LENSMC: "She-Lensmc-Shear",
    ShearEstimationMethods.KSB: "She-Ksb-Shear",
    ShearEstimationMethods.REGAUSS: "She-Regauss-Shear",
    ShearEstimationMethods.MOMENTSML: "She-Momentsml-Shear"
}

METHOD_CATALOG_PATHS = {
    ShearEstimationMethods.LENSMC: "Data.LensMcShearMeasurements.DataStorage.DataContainer.FileName",
    ShearEstimationMethods.KSB: "Data.KsbShearMeasurements.DataStorage.DataContainer.FileName",
    ShearEstimationMethods.REGAUSS: "Data.Regauss.DataStorage.DataContainer.FileName",
    ShearEstimationMethods.MOMENTSML: "Data.MomentsMl.DataStorage.DataContainer.FileName"
}

logger = getLogger(__name__)


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    # Pipeline args
    parser.add_argument("--workdir", type=dir_path, default=".", help="Workdir")
    parser.add_argument("--logdir", type=dir_path, default=".", help="Logdir")

    # Input args
    parser.add_argument(
        '--shear_estimates_product_listfile', type=str,
        required=True, help="Listfile of shear measurements"
    )
    parser.add_argument('--she_lensmc_chains_listfile', type=str, help="Listfile of LensMC chains (optional)")

    # Output args
    parser.add_argument(
        '--merged_she_measurements', type=str, required=True,
        help="Filename of output shear measurements product"
    )
    parser.add_argument(
        '--merged_she_lensmc_chains', type=str,
        help='XML data product to contain LensMC chains data (optional)'
    )

    # Options
    parser.add_argument('--num_procs', type=int, help='Number of parallel processes to use')

    return parser


def read_lensmc_chains_tables(she_lensmc_chains_table_product_filename, workdir):
    try:

        logger.debug("Loading chains from file: " + she_lensmc_chains_table_product_filename)

        # Read in the product and get the filename of the table

        she_lensmc_chains_product = read_xml_product(
            os.path.join(workdir, she_lensmc_chains_table_product_filename))

        if not isinstance(she_lensmc_chains_product, products.she_lensmc_chains.dpdSheLensMcChains):
            raise TypeError(f"Chains product is of invalid type: {type(she_lensmc_chains_product)}")

        she_lensmc_chains_table_filename = she_lensmc_chains_product.get_filename()

        if she_lensmc_chains_table_filename is None or she_lensmc_chains_table_filename == "None":

            logger.debug("No chains avaialble from file: " + she_lensmc_chains_table_product_filename)

        else:

            she_lensmc_chains_table = table.Table.read(
                os.path.join(workdir, she_lensmc_chains_table_filename)
            )

            if not is_in_format(she_lensmc_chains_table, lmcc_tf, verbose=True, ignore_metadata=True):
                raise TypeError("Input chains table is of invalid format.")

    except Exception as e:

        logger.warning("Cannot read chains table from %s: %s", she_lensmc_chains_table_product_filename, str(e))
        return

    logger.debug("Finished loading chains from file: " + she_lensmc_chains_table_product_filename)

    return she_lensmc_chains_table


def read_method_estimates_tables(she_measurements_table_product_filename, workdir):

    try:

        logger.debug("Loading shear estimates from file: " + she_measurements_table_product_filename)

        # Read in the product and get the filename of the table

        she_measurements_table_product = read_xml_product(
            os.path.join(workdir, she_measurements_table_product_filename))

        if not isinstance(she_measurements_table_product, products.she_measurements.dpdSheMeasurements):
            raise TypeError(f"Shear product is of invalid type: {type(she_measurements_table_product)}")

        observation_ids = she_measurements_table_product.Data.ObservationIdList
        observation_time = she_measurements_table_product.Data.ObservationDateTime
        pointing_id_list = she_measurements_table_product.Data.PointingIdList
        tile_list = she_measurements_table_product.Data.TileIndexList
        spatial_coverage = she_measurements_table_product.Data.SpatialCoverage

    except Exception as e:

        logger.warning("Error reading %s: %s", she_measurements_table_product_filename, str(e))
        return None, None, None, None, None, None

    # Loop over methods and read in the table

    she_measurements_tables = {}

    for method in ShearEstimationMethods:

        try:

            she_measurements_method_table_filename = she_measurements_table_product.get_method_filename(method)

            if she_measurements_method_table_filename is None or she_measurements_method_table_filename == "None":
                logger.debug(f"No shear estimates available for method {method}")
                continue

            she_measurements_method_table = table.Table.read(os.path.join(
                workdir, she_measurements_method_table_filename))

            if not is_in_format(
                she_measurements_method_table, sm_tf, verbose=False, ignore_metadata=True, strict=False
            ):
                raise TypeError(f"Input shear estimates table for method {method} is of invalid format.")

            # Append the table to the list of tables
            she_measurements_tables[method] = she_measurements_method_table

        except Exception as e:
            logger.warning("Error occurred accessing %s measurements: %s", method, str(e))

        logger.debug("Finished loading shear estimates from file: " + she_measurements_table_product_filename)

    return (she_measurements_tables, observation_ids, observation_time, pointing_id_list, tile_list, spatial_coverage)


def she_measurements_merge_from_args(args):
    """ Core function for implementing a merge of shear estimates tables
    """

    logger.debug('# Entering she_measurements_merge_from_args(args)')

    # Sanity check arguments (either input and output chains are both None, or both not None)
    if bool(args.she_lensmc_chains_listfile) != bool(args.merged_she_lensmc_chains):
        raise ValueError(
            "Inconsistent input arguments. Input chains listfile is %s but output chains product is %s" % (
                args.she_lensmc_chains_listfile, args.merged_she_lensmc_chains
            )
        )

    # Determine how many processes we'll use
    # (If this is not set, assume one process)
    num_procs = args.num_procs if args.num_procs else 1

    # Keep a list of all shear estimates tables for each method
    she_measurements_tables = dict.fromkeys(ShearEstimationMethods)
    for method in ShearEstimationMethods:
        # Start with an empty list of the tables
        she_measurements_tables[method] = []

    logger.info("Loading shear estimates from files listed in: " + args.shear_estimates_product_listfile)

    measurements_product_filenames = read_listfile(
        os.path.join(args.workdir, args.shear_estimates_product_listfile))

    # Read the measurements tables

    input_tuples = [
        (she_measurements_table_product_filename, args.workdir) for
        she_measurements_table_product_filename in measurements_product_filenames
    ]

    with mp.Pool(processes=num_procs) as pool:
        she_measurements_tables_and_metadata = pool.starmap(read_method_estimates_tables, input_tuples)

    (
        full_l_she_measurements_tables,
        full_l_observation_ids,
        full_l_observation_times,
        full_l_pointing_id_lists,
        full_l_tile_lists,
        full_l_spatial_coverage
    ) = zip(*she_measurements_tables_and_metadata)

    l_she_measurements_tables = [x for x in full_l_she_measurements_tables if x is not None]
    l_observation_ids = [x for x in full_l_observation_ids if x is not None]
    l_observation_times = [x for x in full_l_observation_times if x is not None]
    l_pointing_id_lists = [x for x in full_l_pointing_id_lists if x is not None]
    l_tile_lists = [x for x in full_l_tile_lists if x is not None]
    l_spatial_coverage = [x for x in full_l_spatial_coverage if x is not None]

    # Get list of unique observations, pointings, etc for this run
    observation_ids = list(set(itertools.chain(*l_observation_ids)))
    pointing_id_list = list(set(itertools.chain(*l_pointing_id_lists)))
    tile_list = list(set(itertools.chain(*l_tile_lists)))
    # NOTE: observation time makes no sense when there are more than one observation being used in a tile
    observation_time = l_observation_times[0]
    # All spatial coverages should be the same, so use the first from the list
    spatial_coverage = l_spatial_coverage[0]

    # Sort the tables into the expected format
    for method in ShearEstimationMethods:

        for i in range(len(l_she_measurements_tables)):
            if method not in l_she_measurements_tables[i]:
                continue
            t = l_she_measurements_tables[i][method]
            if t is None or len(t) == 0:
                continue
            she_measurements_tables[method].append(t)

    logger.info("Finished loading shear estimates from files listed in: " + args.shear_estimates_product_listfile)

    logger.info("Combining shear estimates tables.")

    # Combine the shear estimates tables for each method and output them

    combined_she_measurements_tables = dict.fromkeys(ShearEstimationMethods)

    # Create the output products
    combined_she_measurements_product = products.she_measurements.create_she_measurements_product()

    # Set the metadata for the measurements product
    combined_she_measurements_product.Data.ObservationIdList = observation_ids
    combined_she_measurements_product.Data.ObservationDateTime = observation_time
    combined_she_measurements_product.Data.PointingIdList = pointing_id_list
    combined_she_measurements_product.Data.TileIndexList = tile_list
    combined_she_measurements_product.Data.SpatialCoverage = spatial_coverage
    combined_she_measurements_product.Data.CatalogDescription = []
    combined_she_measurements_product.Data.SpectralCoverage = None
    combined_she_measurements_product.Data.PatchIdList = None
    combined_she_measurements_product.Data.CalblockIdList = None
    combined_she_measurements_product.Data.CalblockVariantList = None
    combined_she_measurements_product.Data.IslandLabelList = None
    combined_she_measurements_product.Data.PlaceholderData = None
    combined_she_measurements_product.Data.NumberExposures = len(pointing_id_list)
    combined_she_measurements_product.Data.BatchIndex = None

    # Create (and write to disk) the output tables for each method (if available)
    for method in ShearEstimationMethods:

        # Skip if no data for this method
        if len(she_measurements_tables[method]) == 0:
            combined_she_measurements_product.set_method_filename(method, None)
            continue

        catalog_description = DmUtils.create_catalog_description(
            METHOD_CATALOG_NAMES[method], METHOD_CATALOG_PATHS[method]
        )
        combined_she_measurements_product.Data.CatalogDescription.append(catalog_description)

        # Combine the tables
        combined_she_measurements_tables[method] = table.vstack(
            she_measurements_tables[method],
            metadata_conflicts="silent"
        )

        # Get a filename for the table
        combined_she_measurements_table_filename = get_allowed_filename(
            type_name="SHEAR-EST-" + method.name,
            instance_id='MERGED',
            extension=".fits.gz",
            version=SHE_CTE.__version__,
            subdir="data",
            processing_function="SHE"
        )
        combined_she_measurements_product.set_method_filename(method, combined_she_measurements_table_filename)

        logger.info(
            "Combined shear estimates for method %s output to %s",
            method.value, combined_she_measurements_table_filename
        )

        # Output the combined table
        combined_she_measurements_tables[method].write(
            os.path.join(args.workdir, combined_she_measurements_table_filename),
            format="fits"
        )

    # Save the measurement product
    write_xml_product(combined_she_measurements_product, args.merged_she_measurements, args.workdir)
    logger.info("Combined shear estimates product output to: " + args.merged_she_measurements)

    # Process LensMC Chains (if they are provided)
    if args.she_lensmc_chains_listfile is not None:
        logger.info("Loading chains from files listed in: " + args.she_lensmc_chains_listfile)

        chains_product_filenames = read_listfile(
            os.path.join(args.workdir, args.she_lensmc_chains_listfile))

        # Read the LensMC chains tables
        input_tuples = [
            (she_lensmc_chains_table_product_filename, args.workdir) for
            she_lensmc_chains_table_product_filename in chains_product_filenames
        ]

        with mp.Pool(processes=num_procs) as pool:
            she_lensmc_chains_tables = pool.starmap(read_lensmc_chains_tables, input_tuples)

        logger.info("Finished loading chains from files listed in: " + args.she_lensmc_chains_listfile)

        she_lensmc_chains_tables = [t for t in she_lensmc_chains_tables if t is not None]

        # Combine the chains tables
        combined_she_lensmc_chains_tables = table.vstack(
            she_lensmc_chains_tables,
            metadata_conflicts="silent"
        )

        # Get a filename for the table
        combined_she_lensmc_chains_table_filename = get_allowed_filename(
            type_name="SHEAR-CHAIN",
            instance_id='MERGED',
            extension=".fits.gz",
            version=SHE_CTE.__version__,
            subdir="data",
            processing_function="SHE"
        )

        # Create the output product
        combined_she_lensmc_chains_product = products.she_lensmc_chains.create_lensmc_chains_product()

        # Set the metadata for the chains product
        combined_she_lensmc_chains_product.Data.ObservationIdList = observation_ids
        combined_she_lensmc_chains_product.Data.ObservationDateTime = observation_time
        combined_she_lensmc_chains_product.Data.PointingIdList = pointing_id_list
        combined_she_lensmc_chains_product.Data.TileIndexList = tile_list
        combined_she_lensmc_chains_product.Data.SpatialCoverage = spatial_coverage
        combined_she_lensmc_chains_product.Data.SpectralCoverage = None
        combined_she_lensmc_chains_product.Data.PatchIdList = None
        combined_she_lensmc_chains_product.Data.CalblockIdList = None
        combined_she_lensmc_chains_product.Data.CalblockVariantList = None
        combined_she_lensmc_chains_product.Data.IslandLabelList = None
        combined_she_lensmc_chains_product.Data.PlaceholderData = None
        combined_she_lensmc_chains_product.Data.NumberExposures = len(pointing_id_list)
        combined_she_lensmc_chains_product.Data.BatchIndex = None
        combined_she_lensmc_chains_product.Data.ShearEstimationMethod = "LensMC"

        catalog_description = DmUtils.create_catalog_description(
            "She-Lensmc-Chains", "Data.DataStorage.DataContainer.FileName"
        )
        combined_she_lensmc_chains_product.Data.CatalogDescription = [catalog_description]

        combined_she_lensmc_chains_product.set_filename(combined_she_lensmc_chains_table_filename)

        # Output the combined chains table
        combined_she_lensmc_chains_tables.write(
            os.path.join(args.workdir, combined_she_lensmc_chains_table_filename),
            format="fits"
        )

        write_xml_product(combined_she_lensmc_chains_product, args.merged_she_lensmc_chains, args.workdir)
        logger.info("Combined chains product output to: " + args.merged_she_lensmc_chains)
    else:
        logger.info("Not merging LensMC chains as none were provided")

    logger.debug('# Exiting she_measurements_merge_from_args normally')


def mainMethod(args):
    """
    @brief
        The "main" method for this program, execute a pipeline.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger.info('#')
    logger.info('# Entering SHE_CTE_ShearEstimatesMerge mainMethod()')
    logger.info('#')

    she_measurements_merge_from_args(args)

    logger.info('#')
    logger.info('# Exiting SHE_CTE_ShearEstimatesMerge mainMethod()')
    logger.info('#')


def main():
    """
    @brief
        Alternate entry point for non-Elements execution.
    """

    parser = defineSpecificProgramOptions()

    args = parser.parse_args()

    mainMethod(args)


if __name__ == "__main__":
    main()
