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
File: python/SHE_CalibMLCore/utils.py

Created on: 09/05/17
Author: Malte Tewes & Thibault Kuntzer

General helpers
"""


import gzip
import logging
import os
import pickle
import sys

import numpy as np


logger = logging.getLogger(__name__)


def writepickle(obj, filepath, protocol=-1):
    """
    I write your python object obj into a pickle file at filepath.
    If filepath ends with .gz, I'll use gzip to compress the pickle.
    Leave protocol = -1 : I'll use the latest binary protocol of pickle.
    """
    if os.path.splitext(filepath)[1] == ".gz":
        pkl_file = gzip.open(filepath, 'wb')
    else:
        pkl_file = open(filepath, 'wb')

    pickle.dump(obj, pkl_file, protocol)
    pkl_file.close()
    logger.info("Wrote %s" % filepath)


def readpickle(filepath, encoding="latin1"):
    """
    I read a pickle file and return whatever object it contains.
    If the filepath ends with .gz, I'll unzip the pickle file.
    """

    # Inspecting sys.modules, when it complains that it does not find tenbilac objects (see __init__.py)
    # for (key, value) in sys.modules.items():
    #    if ("tenbilac" in key) or ("SHE_" in key):
    #        print("{} : {}".format(key, value))

    logger.debug("Reading '{}'...".format(filepath))
    if os.path.splitext(filepath)[1] == ".gz":
        pkl_file = gzip.open(filepath, 'rb')
    else:
        pkl_file = open(filepath, 'rb')
    obj = pickle.load(pkl_file, encoding=encoding)
    logger.warning("Encoding hack is in place!")
    pkl_file.close()
    logger.info("Read %s" % filepath)
    return obj


def sigma_clip_plus(arr, maxdev, get_indices=False):
    """
    Removes iteratively all points that are upwards of 

        mean(arr) + maxdev * std(arr)

    :param arr: a numpy array (1D) containing the data points.
    :param maxdev: maximum allowed deviation
    :param get_indices: if True returns arr, keys otherwise just keys
    """

    lenp = 1
    lena = 0

    keys = np.arange(len(arr))
    arr = np.asarray(arr)

    while lenp > lena:
        _thr = np.mean(arr) + maxdev * np.std(arr)
        lenp = np.size(arr)
        sel = arr < _thr
        keys = keys[sel]
        arr = arr[sel]
        lena = np.size(arr)

        if lena == 1:
            break

    if get_indices:
        return arr, keys
    else:
        return arr
