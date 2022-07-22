""" @file utils.py

    Created 25 Aug 2021
    

    Contains utilities for the pipeline
"""

__updated__ = "2021-08-25"

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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA

import json
import os
import logging

from SHE_PPT.logging import getLogger

logger=getLogger(__name__)



def get_qualified_filename(workdir,fname):
    """Checks if fname is an absolute path. If so, returns it, otherwise, assumes it is relative to the workdir and returns the path joined to the workdir"""
    if os.path.isabs(fname):
        return fname
    else:
        return os.path.join(workdir,fname)


def check_input(workdir,fname, name):
    """Checks for the existence of an input"""

    qualified_fname = get_qualified_filename(workdir,fname)
    exists = os.path.exists(qualified_fname)
    if exists:
        logger.info("Input %s exists at %s",name,qualified_fname)
    else:
        logger.warn("Input %s does not exist at %s",name,qualified_fname)
        raise FileNotFoundError("%s not found"%qualified_fname)


def create_dummy_output(workdir,fname, name):
    """Creates a dummy output file for a dry run"""

    qualified_fname = get_qualified_filename(workdir, fname)

    with open(qualified_fname,"w") as f:
        f.write("TEST FILE")
    
    logger.info("Written dummy file for output %s at %s",name,qualified_fname)


def create_dummy_output_listfile(workdir, fname, name, count):
    """Creates a dummy listfile for a dry run wiht the correct number of items in it"""

    qualified_fname = get_qualified_filename(workdir, fname)

    qualified_dummyfile = get_qualified_filename(workdir,"dummy_listfile_entry")
    with open(qualified_dummyfile,"w") as f:
        f.write("TEST FILE")

    dummyfilelist= [qualified_dummyfile]*count

    with open(qualified_fname,"w") as f:
        json.dump(dummyfilelist,f)
    
    logger.info("Written dummy listfile for output %s at %s",name,qualified_fname)


def read_config(filename):
    """Reads in the config file for the pipeline, returning a dictionary"""

    #default values
    defaults = {
        "HDF5": 1,
        "chunked": 1,
        "maxbatches": 4,
        "batchsize": 20,
        "memmap": 1,
        "spatial_batching": 0,
        "dry_run": 1,
        "mean_compute_time": 10,
        "compression": 0
    }

    if not os.path.exists(filename):
        raise FileNotFoundError("Config file %s not found"%filename)

    with open(filename,"r") as f:
        lines = f.readlines()

    #remove whitespace and comments
    strippedlines = []
    for line in lines:
        strippedline = line.strip()
        if strippedline != "" and strippedline[0] != "#":
            strippedlines.append(strippedline)
    

    #create config dictionary
    config = {}
    
    #populate it
    for line in strippedlines:
        key, val = line.split("=")
        key = key.strip()
        val = int(val.strip())

        config[key] = val
    

    #make fill in any missing keys with the defaults
    for key in defaults.keys():
        if key not in config:
            logger.warning("Key %s not in %s. Setting to default value of %d",key,filename,defaults[key])
            config[key] = defaults[key]
    
    #check for unexpected keys and warn about their presence
    for key in config.keys():
        if key not in defaults:
            logger.warning("Unexpected key %s in %s. This will be ignored",key,filename)
    

    return config
        
        
