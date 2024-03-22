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
:file: python/SHE_CTE_PipelineUtility/clustering.py

:date: 27/06/23
:author: Gordon Gibb

"""

import numpy as np

from scipy.spatial import KDTree
import networkx as nx

import logging
logger = logging.getLogger(__name__)


def find_groups(x, y, sep):
    """Identifies groups of objects (objects closer to each other than sep)

    This function assumes a Cartesian coordinate system

    Inputs:
     - x: x positions of the objects
     - y: y positions of the objects
     - sep: the maximum separation between objects to be considered for grouping

    Returns:
     - x_group: the x positions of all objects (grouped objects positions are now the group's centre of mass)
     - y_group: the y positions of all objects (grouped objects positions are now the group's centre of mass)
     - group_ids: array containing the group id for each object (-1 means object is not in a group >=0 means object belongs to a group)
    """
    
    # ensure x and y are arrays
    x = np.asarray(x)
    y = np.asarray(y)

    group_ids = np.full(len(x), -1, dtype=int)
    x_group = np.copy(x)
    y_group = np.copy(y)

    # separation is None or zero, so by definition all objects do not belong to a group
    if not sep:
        logger.info("Separation is zero, so no grouping carried out")
        return x_group, y_group, group_ids

    # construct k-d tree, and get the pairs of all points closer that sep
    xy = np.asarray([x, y]).T

    tree = KDTree(xy)
    pairs = tree.query_pairs(sep)

    # now construct a graph from these points, and determine how many discrete clusters there are
    # e.g. (a, b) and (b, c) form a cluster (via the shared point b), but (a, b) and (c, d) do not
    graph = nx.Graph()
    graph.add_edges_from(pairs)

    # NOTE: nx.connected_components returns a generator, which returns sets.
    # We need a list of lists (which can be used as array indices) so convert this generator to list of lists
    groups = [list(inds) for inds in nx.connected_components(graph)]

    logger.info("Found %d pairs of objects which form %d discrete groups", len(pairs), len(groups))

    ns = []

    # give each group an id (starting from 0) and adjust the positions of objects within each group to that group's
    # centre
    for i, inds in enumerate(groups):
        ns.append(len(inds))
        group_ids[inds] = i
        x_group[inds] = x[inds].mean()
        y_group[inds] = y[inds].mean()

    # calculate and display some statistics
    n_groups = len(groups)
    n_grouped = np.sum(group_ids > -1)
    frac = n_grouped / len(group_ids)

    logger.info("In total %d objects belong to groups: %f%% of all objects", n_grouped, frac * 100)
    if n_grouped:
        
        logger.info("Min group size: %d, max group size: %d, mean group size: %f", np.min(ns), np.max(ns), np.mean(ns))

    return x_group, y_group, group_ids


def spatially_batch(x, y, target_batchsize):
    """Spatially batches objects into batches of with a target batchsize.

    Uses a k-d tree to partition the objects in space, then traverses the tree's nodes, defining a batch as a node
    with <= target_batchsize objects in it.

    NOTE: The returned batchsize will be between target_batchsize/2 and target_batchsize.
    NOTE: If target_batchsize is zero, then all objects are placed into one batch

    Inputs:
     - x: the x positions of objects (assumes Cartesian coordinates)
     - y: the y positions of objects (assumes Cartesian coordinates)
     - target_batchsize: the batch size to aim for.

    Returns:
     - batches: a list of lists of the indices of objects belonging to each batch"""

    def _divide_to_subbatches(node, max_batchsize, batches):
        if node.children <= max_batchsize or node.split_dim == -1:
            # Node either has fewer than max_batchsize or the node is a leaf.
            # All of this node's objects are considered a batch
            if node.children > 0:
                batches.append(node.indices)
        else:
            # recursively call this on the child nodes
            _divide_to_subbatches(node.lesser, max_batchsize, batches)
            _divide_to_subbatches(node.greater, max_batchsize, batches)

    n_objs = len(x)

    logger.info("Target batchsize = %d", target_batchsize)

    if target_batchsize < 0:
        raise ValueError("Target batchsize must be non-zero and positive")

    elif target_batchsize == 0:
        logger.info("Placing all objects into one batch of size %d", n_objs)
        all_inds = np.asarray([i for i in range(n_objs)])
        return [all_inds]

    xy = np.asarray([x, y]).T
    tree = KDTree(xy)

    root_node = tree.tree._node

    # NOTE: _divide_to_subbatches populates batches
    batches = []
    _divide_to_subbatches(root_node, target_batchsize, batches)

    ns = [len(b) for b in batches]
    act_batchsize = np.mean(ns)

    logger.info("Divided %d objects into %d batches of mean size %f", len(x), len(batches), np.mean(ns))
    logger.info(
        "Minumum batchsize: %d, maximum batchsize: %d, median batchsize: %d", np.min(ns), np.max(ns), np.median(ns)
    )

    return batches
