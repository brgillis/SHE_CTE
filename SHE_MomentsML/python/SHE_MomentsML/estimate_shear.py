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
File: python/SHE_MomentsML/estimate_shear.py

Created on: 08/24/17
Author: Malte Tewes
"""

import logging
import os

from ElementsKernel.ProjectCommonRoutines import getAuxPathFile
from SHE_MomentsML import meas_run
from SHE_MomentsML import ml_run
from SHE_MomentsML import params
from SHE_MomentsML import summarize_tables

from SHE_PPT.table_formats.shear_estimates import initialise_shear_estimates_table, tf as setf

stamp_size = 256
x_buffer = -5
y_buffer = -5

# from ElementsKernel.Auxiliary import getAuxiliaryPath # This is for LODEEN 2, it seems
logger = logging.getLogger(__name__)


def estimate_shear(frame_stack, training_data=None, calibration_data=None, workdir=".", *args, **kwargs):
    """Top-level function to run the MomentsML aka MegaLUT shape measurement

    The estimation happens in 3 steps:
    - measure features by looping over sources in the detections table for each exposure
    - predict shear estimates via the previously trained machine learning
    - summarize results per galaxy

    Parameters
    ----------
    frame_stack : SHEFrameStack
        All the data (image, detections table, psfs)
    training_data : MomentsMLParams
        All the method parameters
    calibration_data : None
        Currently unused; needed for proper interface
    workdir : string
        Working directory

    Returns
    -------
    astropy.table.Table
        A table containing one row for each galaxy, in the shear_estimates table format.
    """

    if training_data is None:  # we load some default configuration

        # https://euclid.roe.ac.uk/projects/elements/wiki/Release40
        # https://euclid.roe.ac.uk/projects/elements/wiki/Release521
        # https://euclid.roe.ac.uk/projects/elements/wiki/FAQ

        default_params_path = getAuxPathFile("SHE_MomentsML/default_params")
        training_data = params.MomentsMLParams.read_dir(default_params_path)

        #training_data = params.SHE_MomentsML.params.MomentsMLParams.read_pickle("test.pkl")

    # Measuring features
    feature_tables = meas_run.meas_frame_stack(frame_stack,
                                         params=training_data)

    # Running machine-learning predictions on these feature tables
    pred_tables = [ml_run.predict(feature_table, training_data) for feature_table in feature_tables]

    # Grouping results and collecting statistics
    pred_table = summarize_tables.join(pred_tables)

    # And return a Table
    return pred_table
