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
File: python/SHE_MomentsML/meas_run.py

Created on: 09/07/17
Author: Malte Tewes
"""


import logging

from SHE_MomentsML import meas_stamp
from SHE_MomentsML import utils_table

from SHE_PPT.table_formats import detections as detf


logger = logging.getLogger(__name__)

stamp_size = 128


def meas_frame_stack(frame_stack, params):
    """Measures features on one frame stack (and related PSFs), by looping over the sources from detections_table

    Parameters
    ----------
    frame_stack : SHEframe_stack
        The object containing a science image, detections table, and PSF
    params : MomentsMLParams
        The object containing all the method parameters


    Returns
    -------
    astropy Table
        A copy of the detections_table, with measurements

    """

    #logger.info("Measuring image {} with params {}".format(frame_stack, params.meas))

    # exit()
    

    logger.info("Starting the measurement of {} sources on image with {} exposures...".format(
        len(frame_stack.detections_catalogue), len(frame_stack.exposures) ))
    logger.info("Using params: {}".format(params.meas.items("setup")))

    
    # The following loops seem very sub-optimal. This is a quick and dirty solution, it will need much more attention and better code!
    # Prepare one meas_table for each exposure
    meas_tables = []
    for exposure in frame_stack.exposures:

        # We'll work on a copy of the detections table
        meas_table = utils_table.make_masked_copy(frame_stack.detections_catalogue)

        # Adding the required columns:
        meas_stamp.prep_galsim_adamom(meas_table)
        meas_stamp.prep_skystats(meas_table)
        meas_stamp.prep_snr(meas_table)

        meas_tables.append(meas_table)

    # And looping over the catalog
    for row in frame_stack.detections_catalogue:

        #print(row.index)
        #print(row)
        stamp_stack = frame_stack.extract_galaxy_stack(row[detf.tf.ID],
                                                     width=params.meas.getint("setup", "stampsize"))

        # Skip if no data for this galaxy
        if stamp_stack.is_empty():
            continue

        # Loop over exposures
        for (stamp_exposure_index, stamp_exposure) in enumerate(stamp_stack.exposures):
            #stamp = stamp_stack.exposures[exposure_index]
            if stamp_exposure is None:
                continue
            src = meas_tables[stamp_exposure_index][row.index]
            meas_stamp.galsim_adamom(stamp_exposure, src)
            meas_stamp.skystats(stamp_exposure, src)
            meas_stamp.snr(stamp_exposure, src, gain=params.meas.getfloat("setup", "gain"))

    
    #for meas_table in meas_tables:
    #        print(meas_table["adamom_x", "adamom_y", "adamom_flux", "adamom_snr", "skystats_med"])
    #exit()
    
    return meas_tables
