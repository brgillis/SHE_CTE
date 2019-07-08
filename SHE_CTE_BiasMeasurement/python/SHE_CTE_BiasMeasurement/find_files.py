""" @file find_files.py

    Created 20 Aug 2018

    Functions to help find archived shear statistics files
"""

__updated__ = "2019-07-08"

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
from SHE_PPT.logging import getLogger
import multiprocessing as mp

shear_statistics_filename = "shear_bias_statistics.xml"


def recursive_find_files(base_dir=".", required_filename=shear_statistics_filename, endpoint=True,
                         number_threads=1):
    """Search recursively through the provided directory to find any files matching the desired filename.
    """

    logger = getLogger(__name__)

    logger.debug("Entering recursive_find_files with base_dir = " + base_dir)

    # Silently remove trailing slash from base_dir if present
    if base_dir[-1] == "/":
        base_dir = base_dir[0:-1]

    # Start an empty list of all files we've found
    files_found = []

    files_and_dirs = os.listdir(base_dir)

    # Loop through this directory and add any we find, including recursively
    dirs = []
    for file_or_dir in files_and_dirs:
        qualified_name = os.path.join(base_dir, file_or_dir)

        # If it's a file, does it match the required pattern?
        if os.path.isfile(qualified_name):
            if file_or_dir == required_filename:
                # It does, so add it to the list
                files_found.append(qualified_name)
        # If it's a directory, search through it for files
        if os.path.isdir(qualified_name):
            dirs.append(qualified_name)

    # Look through directories now. If number_threads is 1, do this simply, otherwise use a pool
    if number_threads == 1:
        for qualified_name in dirs:
            more_files_found = recursive_find_files(qualified_name, required_filename, endpoint=False)
            files_found += more_files_found
    else:
        pool = mp.Pool(processes=number_threads)
        pool_results = [pool.apply(recursive_find_files, args=(
            qualified_name, required_filename, False, 1)) for qualified_name in dirs]
        for more_files_found in pool_results:
            files_found += more_files_found

    # If we're at the endpoint, remove the base dir from all files in the list
    if endpoint:
        logger.debug("At endpoint of recursive_find_files, so cleaning up base_dir from filenames (" +
                     str(files_found) + ")")
        for i in range(len(files_found)):
            files_found[i] = files_found[i].replace(base_dir + "/", "")

    logger.debug("Exiting recursive_find_files with files found: " + str(files_found))

    # Return the list of files found
    return files_found
