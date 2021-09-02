""" @file combine_split_fits_listfile.py

    Created 25 Aug 2021
    

    Combines files in the listfile from SplitFits into one .
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


import os
import json

from SHE_PPT.logging import getLogger

from .utils import read_config, get_qualified_filename, check_input, create_dummy_output



logger = getLogger(__name__)



def combine_split_fits_listfile_from_args(args):

    config=read_config(get_qualified_filename(args.workdir,args.pipeline_config))
    dry_run = bool(config["dry_run"])

    if dry_run:
        #check all the input files exist, create dummy output files, then exit
        logger.info("Dry run:")

        check_input(args.workdir,args.input_listfile,"input_listfile")
        create_dummy_output(args.workdir,args.output_json,"output_json")

        return
    
    workdir = args.workdir
    listfile = args.input_listfile
    
    #read the input listfile
    with open(get_qualified_filename(workdir,listfile),"r") as f:
        input_files = json.load(f)

    output_dict = {}
    
    #combine the dictionaries within each file into one master dict
    for input_file in input_files:
        if type(input_file) is list:
            input_file=input_file[0]
        with open(get_qualified_filename(workdir,input_file),"r") as f:
            ccd_dict = json.load(f)
        
        #merge the dictionary from the input file into the final_dict
        output_dict.update(ccd_dict)
    
    #then write this to file
    outfile=get_qualified_filename(workdir,args.output_json)
    logger.info("Writing combined file to %s",outfile)
    with open(outfile,"w") as f:
        json.dump(output_dict,f,indent=4)

    
    






