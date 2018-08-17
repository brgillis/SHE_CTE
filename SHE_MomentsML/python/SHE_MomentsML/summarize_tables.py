#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

"""
File: python/SHE_MomentsML/summarize_tables.py

Created on: 11/07/17
Author: user
"""

import logging

import astropy
import numpy as np

from SHE_PPT.table_formats.detections import tf as dtf
from SHE_PPT.table_formats.shear_estimates import initialise_shear_estimates_table, tf as setf
from SHE_PPT.table_utility import is_in_format

from . import utils_table


logger = logging.getLogger(__name__)


def join(table_list):
    """

    """
    logger.info("Joining {} tables...".format(len(table_list)))

    #joined_table = table_list[0]  # Implement this...
    #print(utils_table.get_info(joined_table))
    # logger.info(utils_table.get_info(joined_table))

    #joined_table["G1"] = joined_table["pre_s1"]
    #joined_table["G2"] = joined_table["pre_s2"]
    #joined_table["G1_ERR"] = 0.5
    #joined_table["G2_ERR"] = 0.5

    #joined_table.keep_columns(["SOURCE_ID", "G1", "G2", "G1_ERR", "G2_ERR"])
    #shear_estimates_table = initialise_shear_estimates_table(
    #    detector=int(joined_table.meta[dtf.m.extname].split(".")[0]))
    
    shear_estimates_table = initialise_shear_estimates_table()

    # logger.info(utils_table.get_info(shear_estimates_table))
    
    #for table in table_list:
    #    print(utils_table.get_info(table))
    #    print(table["pre_s1", "pre_s2", "adamom_snr", "adamom_flux", "adamom_x", "adamom_y"])
    
    #exit()
    for row_index in range(len(table_list[0])):
        
        object_id = table_list[0]["ObjectId"][row_index]
        
        #print("Object {}:".format(object_id))
        #for table in table_list:
        #    print(table["pre_s1", "pre_s2", "adamom_snr", "adamom_flux"][row_index])
        
        s1_mean = np.mean([table["pre_s1"][row_index] for table in table_list])
        s2_mean = np.mean([table["pre_s2"][row_index] for table in table_list])
        s1w_mean = np.mean([table["pre_s1w"][row_index] for table in table_list])
        s2w_mean = np.mean([table["pre_s2w"][row_index] for table in table_list])
      
        sigma_shape = 0.25
        s1_err = sigma_shape / np.sqrt(s1w_mean)
        s2_err = sigma_shape / np.sqrt(s2w_mean)
       
        snr = np.mean([table["adamom_snr"][row_index] for table in table_list])
                       
        shear_estimates_table.add_row({setf.ID: object_id,
                                       setf.g1: s1_mean,
                                       setf.g2: s2_mean,
                                       setf.g1_err: s1_err,
                                       setf.g2_err: s2_err,
                                       setf.snr: snr
                                       })

    # logger.info(utils_table.get_info(shear_estimates_table))
    # logger.info(str(shear_estimates_table))

    if not is_in_format(shear_estimates_table, setf):
        raise ValueError("Output table not in correct format!")

    logger.info("Returning joined table with {} rows...".format(len(shear_estimates_table)))
    return shear_estimates_table
