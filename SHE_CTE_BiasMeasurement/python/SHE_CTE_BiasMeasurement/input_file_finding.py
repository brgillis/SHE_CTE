""" @file input_file_finding.py

    Created 7 Apr 2017

    Functions for locating input files for bias measurement.

    ---------------------------------------------------------------------

    Copyright (C) 2017 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os

from icebrgpy.logging import getLogger

from SHE_CTE_BiasMeasurement import magic_values as mv
from SHE_CTE_ShearEstimation import magic_values as est_mv

def get_input_files(root_dir, required_input_pattern=None, depth=0, input_files=None):
    """
    @brief
        Get input files for bias measurement found in a given directory.
    
    @param root_dir <str>
    @param required_input_pattern <str>
    @param depth <int> Subfolder depth to search.
    
    @return None
    """
    
    logger = getLogger(mv.logger_name)
    logger.debug("Entering get_input_files in directory " + root_dir + ".")
    
    if input_files is None:
        input_files = []
    
    # Start in the base path
    files_and_dirs_in_path = os.listdir(root_dir)

    for file_or_dir_name in files_and_dirs_in_path:

        joined_name = os.path.join(root_dir, file_or_dir_name)

        # Check if this is a file

        if os.path.isfile(joined_name):
            
            # Check for the required input pattern if necessary
            if required_input_pattern is not None:
                if required_input_pattern not in file_or_dir_name:
                    continue

            # Check if this name corresponds to a valid measurements file
            if est_mv.output_tail not in file_or_dir_name:
                continue
            
            test_details_file_name = None
            
            # Check if a corresponding details file exists
            test_details_file_name_1 = joined_name.replace(est_mv.output_tail,"_details.fits")
            if os.path.isfile(test_details_file_name_1):
                test_details_file_name = test_details_file_name_1
            else:
                # Try using the required_input_pattern as well 
                test_details_file_name_2 = joined_name.replace(required_input_pattern+est_mv.output_tail,"_details.fits")
                if os.path.isfile(test_details_file_name_2):
                    test_details_file_name = test_details_file_name_2
                else:
                    test_details_file_name_3 = joined_name.replace("_"+required_input_pattern+est_mv.output_tail,"_details.fits")
                    if os.path.isfile(test_details_file_name_3):
                        test_details_file_name = test_details_file_name_3
                    else:
                        logger.warning("Shear measurements file " + joined_name + " has no corresponding details file. " + 
                                       "(Expected " + str(test_details_file_name_1) + ", " +
                                       str(test_details_file_name_2) + ", or " +
                                       str(test_details_file_name_3) + ".)")

            if test_details_file_name is not None:
                input_files.append((joined_name,test_details_file_name))
        else:

            # This is a directory, so let's go into it and look for files there if we
            # don't exceed depth
            if depth > 0:
                get_input_files(joined_name, required_input_pattern=required_input_pattern,
                                depth=depth-1, input_files=input_files)
            
    logger.debug("Exiting get_input_files from directory " + root_dir + ".")
    
    return input_files
