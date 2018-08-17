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
File: python/SHE_MomentsML/ml_run.py

Created on: 09/07/17
Author: user
"""

import logging

import astropy.table

import SHE_CalibMLCore

from . import utils_table


logger = logging.getLogger(__name__)


def predict(cat, params):
    """Computes machine learning predictions from the input the catalog

    Parameters
    ----------
    cat : Astropy Table
        The input catalog
    params : MomentsMLParams
        The object containing all the method parameters

    """

    logger.info("Starting the ML precitions for {} sources".format(len(cat)))
    logger.info("Using params: {}".format(params.ml.items("setup")))

    # Working on a copy:
    outcat = utils_table.make_masked_copy(cat)

    # Sequentially makign one prediction after the other (they might build upon each other!)
    for (dataconfig, tenbilac) in params.tenbilac_conflist:

        # We get the required data config in form of python lists:
        inputlabels = list(eval(dataconfig.get("data", "inputlabels")))
        predlabels = list(eval(dataconfig.get("data", "predlabels")))

        logger.info("Predicting '{}' with '{}'...".format(predlabels, inputlabels))

        # Check that the predictions do not yet exist
        for predlabel in predlabels:
            if predlabel in outcat.colnames:
                raise RuntimeError("The predlabel '{}' already exists in the catalog!".format(predlabel))

        # Preparing the inputs
        inputsdata = utils_table.get_3D_data(outcat, inputlabels)

        # And call Tenbilac prediction
        preddata = tenbilac.predict(inputsdata)

        # We add this data to the outcat.
        # An explicit loop, to highlight that we care very much about the order (to get targetlabels right)
        for (i, predlabel) in enumerate(predlabels):
            logger.info("Adding predictions '{}' to catalog...".format(predlabel))
            data = preddata[:, i, :].transpose()
            assert data.ndim == 2  # Indeed this is now always 2D.
            if data.shape[1] == 1:  # If we have only one realization, just make it a 1D numpy array.
                data = data.reshape((data.size))
                assert data.ndim == 1

            newcol = astropy.table.MaskedColumn(data=data, name=predlabel)
            outcat.add_column(newcol)

    return outcat
