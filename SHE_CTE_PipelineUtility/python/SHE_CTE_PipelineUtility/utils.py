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
:file: python/SHE_CTE_PipelineUtility/utils.py

:date: 2024/01/23
:author: @rrollins

"""

from itertools import islice


def batched(iterable, batch_size):
    """
    Take an input iterator and yield tuples of sequential elements of a given size

    batched('ABCDEFG', 3) --> ('A', 'B', 'C'), ('D', 'E', 'F'), ('G',)

    :param iterable: An input iterable object
    :param batch_size: The number of elements in each batch
    """

    if batch_size < 1:
        raise ValueError('batch_size must be at least one')
    it = iter(iterable)
    while batch := tuple(islice(it, batch_size)):
        yield batch
