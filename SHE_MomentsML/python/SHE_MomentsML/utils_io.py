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
File: python/SHE_MomentsML/utils_io.py

Created on: 09/07/17
Author: Malte Tewes
"""

import gzip
import logging
import os
import pickle


logger = logging.getLogger(__name__)


def write_pickle(obj, filepath, protocol=-1):
    """Write a python object into a pickle file

    Parameters
    ----------
    obj : object
        A python object to be written into the pickle file
    filepath : str
        Filepath to write to.
        If filepath ends with .gz, uses gzip to compress the pickle.
    protocol : int
        Protocol to use. Leave this to -1 : it will use the latest binary protocol of pickle.
    """
    if os.path.splitext(filepath)[1] == ".gz":
        pkl_file = gzip.open(filepath, 'wb')
    else:
        pkl_file = open(filepath, 'wb')

    pickle.dump(obj, pkl_file, protocol)
    pkl_file.close()
    logger.info("Wrote '{}'".format(filepath))


def read_pickle(filepath):
    """Reads a pickle file and returns the object it contains.

    Parameters
    ----------
    filepath : str
        Filepath to read from.
        If filepath ends with .gz, uses gzip to uncompress the pickle.

    """
    if os.path.splitext(filepath)[1] == ".gz":
        pkl_file = gzip.open(filepath, 'rb')
    else:
        pkl_file = open(filepath, 'rb')
    obj = pickle.load(pkl_file)
    pkl_file.close()
    logger.info("Read '{}'".format(filepath))
    return obj
