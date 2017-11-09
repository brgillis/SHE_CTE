""" @file prepare_configs.py

    Created 25 Aug 2017

    Contains primary functions for preparing config files.
"""

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

from SHE_PPT.io import get_allowed_filename, write_listfile
from SHE_PPT.utility import hash_any

def prepare_configs_from_args(args):
    """
        @brief Primary function for preparing configuration files.
        
        @param args Parsed arguments
    """
    
    image_models, full_options = get_image_models_and_options(args)
    
    options_hash = hash_any(frozenset(full_options.items))
    
    model_seed_start = args.model_seed_start
    if args.vary_noise:
        noise_seed_start = args.noise_seed_start
    else:
        noise_seed_start = 0
        
    filenames = []
    
    # Loop through and write a config for each image
    for i in range(len(image_models)):
        
        if args.vary_noise:
            model_seed = model_seed_start
            noise_seed = noise_seed_start + i
        else:
            model_seed = model_seed_start + i
            noise_seed = noise_seed_start
        
        instance_id = args.tag + "-" + str(hash_any((options_hash,model_seed,noise_seed),8))
        
        filename = get_allowed_filename("CTE-CFG",instance_id,".cfg",)
        
        # Write out this config file
        write_config( filename, full_options, model_seed, noise_seed )
        
        filenames.append(filename)
        
    # Save a listfile of the filenames
    write_listfile(args.output_listfile,filenames)
    
    return
        