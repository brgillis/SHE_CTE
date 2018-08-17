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
File: tests/python/layer_test.py

Created on: 09/05/17
Author: Malte Tewes
"""

import py.test

from SHE_CalibMLCore import act
from SHE_CalibMLCore.layer import Layer
import numpy as np


class Test_layer(object):
    """

    """

    def test_run(self):

        layer = Layer(ni=2, nn=3, actfct=act.iden)

        layer.weights[0, 0] = 1     # First index is neuron, second is input (neuron from preceding layer)
        layer.weights[0, 1] = 2
        layer.weights[1, 1] = 10
        layer.weights[2, 0] = 3
        layer.weights[2, 1] = 0
        layer.biases[1] = 42

        #assert layer.report() == """Layer 'None', mode sum, ni 2, nn 3, actfct iden:
        #output 0 = iden ( input * [1. 2.] + 0.0 )
        #output 1 = iden ( input * [ 0. 10.] + 42.0 )
        #output 2 = iden ( input * [3. 0.] + 0.0 )"""

        # We create data to test the 3D case.
        # 2 features, 2 cases, 4 realizations:
        # Case 1 has features close to (1, -1), case 2 has approximatively (4, -4)
        # Order of indices: realization, feature, case
        input = np.array([
            [[1, 4], [-1, -4]],
            [[1.1, 4], [-1, -4]],
            [[1, 4], [-1, -4.1]],
            [[1, 4], [-1, -4]]
        ])

        output = layer.run(input)  # Order of indices: realization, neuron, case
        assert output.shape == (4, 3, 2)

        assert np.allclose(output[0], np.array([[-1, -4], [32, 2], [3, 12]]))
        assert np.allclose(output[1], np.array([[-0.9, -4], [32, 2], [3.3, 12]]))
        assert np.allclose(output[2], np.array([[-1, -4.2], [32, 1], [3, 12]]))
        assert np.allclose(output[3], output[0])

        assert layer.nparams() == 9

        layer.setzero()
        output = layer.run(input)
        assert np.allclose(output, 0)

        layer.addnoise()
        output = layer.run(input)

        layer.setidentity()
        output = layer.run(input)
        assert np.allclose(output[0], np.array([[1, 4], [-1, -4], [0, 0]]))
