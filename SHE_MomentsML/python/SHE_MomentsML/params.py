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
File: python/SHE_MomentsML/params.py

Created on: 09/07/17
Author: Malte Tewes
"""

import os
import configparser
from configparser import SafeConfigParser


import SHE_CalibMLCore.com
from SHE_MomentsML import utils_io




import logging
logger = logging.getLogger(__name__)


class MomentsMLParams(object):
    """Class grouping all parameters required for a MomentsML run
    
    A (maybe temporary) key feature is that this object can be read from a single file, produced at or after the 
    training stage.
    
    Attributes
    ----------
    meas : Config for the feature measurement
    ml : General config for the machine learning
    tenbilac_conflist : list of (dataconf, tenbilac) tuples
    
    
    """
    
    
    def __init__(self, name="Empty"):
        
        self.name = name
        self.meas = None
        self.ml = None
        self.tenbilac_conflist = []
    
    def __str__(self):
        return "MomentsMLParams({}) with {} tenbilacs".format(self.name, len(self.tenbilac_conflist))
    
    
    @classmethod
    def read_dir(cls, dirpath):
        """Assembles the object from a human-readable directory structure containing the config files and tenbilac
        pickles.
        """
        logger.info("Reading MomentsMLParams from '{}'...".format(dirpath))
        name = os.path.split(dirpath)[-1]
        
        new_params = MomentsMLParams(name)
        new_params.meas = read_cfg_file(os.path.join(dirpath, "meas.cfg"))
        new_params.ml = read_cfg_file(os.path.join(dirpath, "ml.cfg"))
        
        # The Sections of the ml config point to the tenbilacs to use.
        tenbilac_steps = [item for item in new_params.ml.sections() if item.startswith("tenbilac_")]
        # We loop through those, and load one config and one tenbilac per section:
        for section in tenbilac_steps:
            dataconf = read_cfg_file(os.path.join(dirpath, new_params.ml.get(section, "data_config_file")))
            tenbilac = SHE_CalibMLCore.com.Tenbilac(
                os.path.join(dirpath, new_params.ml.get(section, "ml_dir"))
                )
            tenbilac.load()
            new_params.tenbilac_conflist.append([dataconf, tenbilac])
        
        return new_params
    
    def write_pickle(self, filepath):
        """Writes the object to disk"""
        utils_io.write_pickle(self, filepath)
    
     
    @classmethod
    def read_pickle(cls, filepath):
        """Reads-in a MomentsMLParams from disk"""
        return utils_io.read_pickle(filepath)
    
    @classmethod
    def make_pickle_from_dir(cls, dirpath, filepath):
        new_params = cls.read_dir(dirpath)
        new_params.write_pickle(filepath)
    
    
def read_cfg_file(filepath):
    """Reads-in a config file and returns a SafeConfigParser
    
    """
    logger.info("Reading '{}'...".format(filepath))
    
    config = SafeConfigParser(allow_no_value=True)
    config.read(filepath)
    try:
        version = config.get("setup", "version")
    except configparser.Error:
        version = None
    sections = config.sections()
    n_sections = len(sections)
    logger.debug("Version is '{}', {} sections: {}".format(version, n_sections, sections))
    
    return config


    
    
    
    