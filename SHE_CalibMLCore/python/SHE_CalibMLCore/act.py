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
File: python/SHE_CalibMLCore/act.py

Created on: 09/05/17
Author: Malte Tewes

Activation functions
"""


import numpy as np

import logging
logger = logging.getLogger(__name__)


def sig(x):
    return 1.0 / (1.0 + np.exp(-x))    # The actual sigmoid


def sigi(x):
    return 1.0 / (1.0 + np.exp(-x * 4.0))    # A sigmoid with a central slope of 1


def sige(x):
    return 2.0 / (1.0 + np.exp(-x)) - 1.0  # Making it even


def tanh(x):
    return np.tanh(x)


def iden(x):
    return x


def relu(x):
    return np.amax([x, np.zeros_like(x)], axis=0)


def relu6(x):
    return np.amin([np.amax([x, np.zeros_like(x)], axis=0), 6 * np.ones_like(x)], axis=0)


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    x = np.linspace(-5, 5, 1000)
    acts = [sig, sigi, sige, tanh, iden, relu]

    for act in acts:
        plt.plot(x, act(x), label=act.__name__)

    plt.xlabel(r"$x$")
    plt.ylabel(r"$f(x)$")
    plt.title("Activation functions")
    plt.ylim(-1.2, 1.2)
    plt.legend()
    plt.grid()
    plt.show()
