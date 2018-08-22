""" @file find_files.py

    Created 20 Aug 2018

    Functions to help find archived shear statistics files
"""

__updated__ = "2018-08-22"

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

shear_statistics_filename = "shear_bias_statistics.xml"


def recursive_find_files(base_dir=".", required_filename=shear_statistics_filename, endpoint=True):
    """Search recursively through the provided directory to find any files matching the desired filename.
    """

    # Silently remove trailing slash from base_dir if present
    if base_dir[-1] == "/":
        base_dir = base_dir[0:-1]

    # Start an empty list of all files we've found
    files_found = []

    files_and_dirs = os.listdir(base_dir)

    # Loop through this directory and add any we find, including recursively
    for file_or_dir in files_and_dirs:
        qualified_name = os.path.join(base_dir, file_or_dir)

        # If it's a file, does it match the required pattern?
        if os.path.isfile(qualified_name):
            if file_or_dir == required_filename:
                # It does, so add it to the list
                files_found.append(qualified_name)
        # If it's a directory, search through it for files
        if os.path.isdir(qualified_name):
            more_files_found = recursive_find_files(qualified_name, required_filename, endpoint=False)
            files_found += more_files_found

    # End "for file_or_dir in files_and_dirs:"

    # If we're at the endpoint, remove the base dir from all files in the list
    if endpoint:
        for i in range(len(files_found)):
            files_found[i] = files_found[i].replace(base_dir + "/", "")

    # Return the list of files found
    return files_found
