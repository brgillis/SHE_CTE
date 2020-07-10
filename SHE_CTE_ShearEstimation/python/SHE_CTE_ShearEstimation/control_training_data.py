""" @file control_training_data.py

    Created 23 July 2018

    Classes and functions related to loading KSB and REGAUSS training data.
"""

__updated__ = "2020-07-10"

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

from SHE_PPT import products
from SHE_PPT.logging import getLogger
from SHE_PPT.pipeline_utility import get_conditional_product
from SHE_PPT.table_formats.she_ksb_training import tf as ksbt_tf
from SHE_PPT.table_formats.she_lensmc_measurements import tf as lmcm_tf
from SHE_PPT.table_formats.she_regauss_training import tf as regt_tf
from SHE_PPT.table_utility import is_in_format
from astropy.table import Table


class ControlTraining(object):
    """An object which contains needed summary values of training data for the KSB and REGAUSS methods.
    """

    def __init__(self, training_product_filename=None, workdir="."):

        logger = getLogger(__name__)
        logger.debug("Entering ControlTraining __init__")

        p = get_conditional_product(training_product_filename, workdir)

        if p is None:
            logger.warn("No training data provided; using default shape noise of 0.")
            self.e1_var = 0.
            self.e2_var = 0.
        else:
            # Read the table stored in the data container of the product
            training_data_filename = p.get_data_filename()
            qualified_training_data_filename = os.path.join(workdir, training_data_filename)
            t = Table.read(qualified_training_data_filename)

            if is_in_format(t, ksbt_tf):  # For KSB
                self.e1_var = t[ksbt_tf.e1].var()
                self.e2_var = t[ksbt_tf.e2].var()
            elif is_in_format(t, regt_tf):  # For REGAUSS
                self.e1_var = t[regt_tf.e1].var()
                self.e2_var = t[regt_tf.e2].var()
            elif is_in_format(t, lmcm_tf):  # For LensMC
                self.e1_var = t[lmcm_tf.g1].var()
                self.e2_var = t[lmcm_tf.g2].var()
            else:
                raise ValueError("Unrecognized table format for training data in file: " + qualified_training_data_filename)

        logger.debug("Exiting ControlTraining __init__")
        return

    @property
    def e1_var(self):
        """The variance of the ellipticity measurements, component 1.
        """
        return self._e1_var

    @e1_var.setter
    def e1_var(self, value):
        """Set the variance of the ellipticity measurements, component 1.
        """
        self._e1_var = value
        return

    @property
    def e2_var(self):
        """The variance of the ellipticity measurements, component 2.
        """
        return self._e2_var

    @e2_var.setter
    def e2_var(self, value):
        """Set the variance of the ellipticity measurements, component 2.
        """
        self._e2_var = value
        return


def load_control_training_data(training_data_filename, workdir="."):
    """Loading method for getting KSB/REGAUSS/LensMC training data.
    """

    return ControlTraining(training_data_filename, workdir)
