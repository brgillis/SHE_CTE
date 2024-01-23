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
:file: tests/python/utils_test.py

:date: 2024/01/23
:author: @rrollins

"""

from hypothesis import given
from hypothesis.strategies import integers, iterables
from itertools import islice, tee
import pytest
import sys

from SHE_CTE_PipelineUtility import utils


class NotIterable:
    pass


@given(iterables(integers()), integers(min_value=1, max_value=sys.maxsize))
def test_batched_integers(iterable, batch_size):
    i1, i2 = tee(iterable)
    for batch in utils.batched(i1, batch_size):
        assert batch == tuple(islice(i2, batch_size))


@given(integers(min_value=1))
def test_batched_not_iterable(batch_size):
    with pytest.raises(TypeError, match='object is not iterable'):
        batched = utils.batched(NotIterable(), batch_size)
        next(batched)


@given(integers(max_value=0))
def test_batched_invalid_batch_size(batch_size):
    with pytest.raises(ValueError, match='batch_size must be at least one'):
        batched = utils.batched(None, batch_size)
        next(batched)
