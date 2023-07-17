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
:file: tests/python/clustering_test.py

:date: 27/06/23
:author: Gordon Gibb

"""
import pytest

import numpy as np
import functools

from SHE_CTE_PipelineUtility.clustering import spatially_batch, find_groups



def test_find_groups():
    """Tests find_groups by generating clouds of points and making sure the correct number of groups
        are identified in them"""

    # baseline number of points in each direction
    nx = ny = 128

    x = []
    y = []

    # generate regular grid of points with separation 1, this is our base set of points for the tests
    for i in range(nx):
        x += [i for j in range(ny)]
        y += [j for j in range(ny)]

    # identify groups in this with separation of 0.9 - should find 0 groups as all points have a separation of >= 1
    xx, yy, groups = find_groups(x, y, sep=0.9)

    # make sure there are 0 groups
    n_groups = groups.max() + 1
    assert (n_groups == 0)

    # now add two points in the middle of squares of points, creating two groups of 5 objects
    x2, y2 = x + [1.5, 10.5], y + [1.5, 10.5]

    xx, yy, groups = find_groups(x2, y2, sep=0.9)

    # make sure there are now two groups
    n_groups = groups.max() + 1
    assert (n_groups == 2)

    counts = np.bincount(groups + 1)
    # make sure number of un-grouped objects is correct
    assert (counts[0] == nx * ny - 2 * 4)
    # make sure the count in each group is correct
    assert (counts[1] == 5)
    assert (counts[2] == 5)

    # now add a chain of nx-1 points in the x-direction to connect a whole line of points into a single group
    x_chain = x + [i + 0.5 for i in range(nx - 1)]
    y_chain = y + [20 for i in range(nx - 1)]

    xx, yy, groups = find_groups(x_chain, y_chain, sep=0.9)

    # ensure we have one group
    n_groups = groups.max() + 1
    assert (n_groups == 1)

    counts = np.bincount(groups + 1)
    # make sure number of un-grouped objects is correct
    assert (counts[0] == nx * ny - nx)
    # make sure the count in each group is correct
    assert (counts[1] == nx * 2 - 1)

    # extreme test case: have separation large enough that all objects are collected into one group
    xx, yy, groups = find_groups(x, y, sep=1.1)

    # ensure we have one group
    n_groups = groups.max() + 1
    assert (n_groups == 1)

def test_spatially_batch():

    # create a unit square of n random points
    n = 10000
    x = np.random.random(n)
    y = np.random.random(n)

    with pytest.raises(ValueError):
        # negative batchsize makes no sense
        spatially_batch(x,y, target_batchsize = -1)

    # target_batchsize=0 - should create one batch with all objects
    batches = spatially_batch(x, y, target_batchsize=0)
    assert len(batches) == 1
    assert len(batches[0]) == n
    
    # target_batchsize of 20 - should get between n/20 and 2*(n/20) batches
    batches = spatially_batch(x, y, target_batchsize=20)
    num_batches_expected = n/20
    assert num_batches_expected < len(batches) < num_batches_expected*2
    
    # make sure we haven't lost any objects
    n_objs = functools.reduce(lambda x, y: x + len(y), batches, 0)
    assert n == n_objs


