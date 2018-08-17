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
File: python/SHE_MomentsML/utils_table.py

Created on: 09/07/17
Author: Malte Tewes
"""


import astropy.table
import copy
import numpy as np

import logging
logger = logging.getLogger(__name__)


def get_info(cat, txt=True):
    """Returns a new table "describing" the content of the table cat.
    
    Parameters
    ----------
    cat : astropy Table
        The table wich should be described
    txt : bool
        If True, some text is returned. If False, an actual Table is returned.
    
    """
    colnames = cat.colnames
    dtypes = [cat[colname].dtype for colname in colnames]
    ndims = [cat[colname].ndim for colname in colnames]
    shapes = [cat[colname].shape for colname in colnames]
    infotable = astropy.table.Table([colnames, dtypes, ndims, shapes], names=("colname", "dtype", "ndim", "shape"))
    
    infotable.sort("colname")
    infotable.meta = cat.meta
    
    if txt:
        
        lines = infotable.pformat(max_lines=-1, max_width=-1)
        lines.append("")
        lines.append("Number of rows: {}".format(len(cat)))
        lines.append("Number of columns: {}".format(len(cat.colnames)))
        lines.append("Metadata: {}".format(str(infotable.meta.items())))
        
        return "\n".join(lines)
    else:
        return infotable
    
    
    

def make_masked_copy(cat):
    """Returns a masked copy of a table"""
    # Convert the table to a masked table 
    # A bit strange: reading the doc I feel that this conversion is not needed.
    # But without it, it just doesn't result in a masked table once the masked columns are appended.
    return astropy.table.Table(copy.deepcopy(cat), masked=True)
 
    
def add_cols(cat, colnames, dtype, prefix="", masked=True):
    """Adds columns to the table, potentially masked, initializing them with all values masked
    
    This can be done to prepare a catalog to take some measurements.
    """
    
    cat_length = len(cat)
    for colname in colnames:
        if masked:
            new_col = astropy.table.MaskedColumn(name=prefix+colname,
                                                 dtype=dtype,
                                                 length=cat_length
                                                 )
        else:
            new_col = astropy.table.Column(name=prefix+colname,
                                           data=np.zeros(cat_length, dtype=dtype),
                                           dtype=dtype)

        cat.add_column(new_col)
            
    # We want to mask all entries. They will get unmasked when values will be attributed.
    if masked:
        for colname in colnames:
            cat[prefix+colname].mask = [True] * cat_length # "True" means masked
        
    
    
def get_3D_data(catalog, colnames):
    """Builds a 3D numpy array (typically for Tenbilac input) from columns of an astropy catalog.
    
    Ensure that all columns get the same shape.
    The 3D output array has shape (realization, feature, case), even if there is only one realization.
    """
    
    logger.debug("Getting 3D data from columns {}".format(colnames))
                 
    if len(colnames) == 0:
        raise RuntimeError("No colnames to get data from!")
    
    # Check for exotic catalogs (do they even exist ?)
    for colname in colnames:
        if not catalog[colname].ndim in [1, 2]:
            raise RuntimeError("Can only work with 1D or 2D columns")
    
    # Let's check the depths of the 2D colums to see what size we need.
    nreas = list(set([catalog[colname].shape[1] for colname in colnames if catalog[colname].ndim == 2]))
    #logger.info("We have the following nreas different from one in there: {}".format(nreas))
    if len(nreas) > 1:
        raise RuntimeError("The columns have incompatible depths!")

    if len(nreas) == 0:
        nrea = 1
        logger.info("For each column, only one realization is available.")
        
    else:
        nrea = nreas[0]
        logger.info("Extracting data from {0} realizations...".format(nrea))
        nrea = nreas[0]
    
    if "ngroup" in catalog.meta:
        if nrea != catalog.meta["ngroup"]:
            raise RuntimeError("Something very fishy: depth is not ngroup!")

    # And now we get the data:
    
    readycols = []
    for colname in colnames:
                
        col = np.ma.array(copy.deepcopy(catalog[colname])) # The deepcopy seems required, at least with Astropy 1.1.1
                
        if col.ndim == 2:
            pass
            
        elif col.ndim == 1:
            # This column has only one realization, and we have to "duplicate" it nrea times...
            col = np.tile(col, (nrea, 1)).transpose()
                    
        else:
            raise RuntimeError("Weird column dimension")
                                
        readycols.append(col)
        
    outarray = np.rollaxis(np.ma.array(readycols), 2)
    
    assert outarray.ndim == 3
    assert outarray.shape[1] == len(colnames)

    return outarray
    
    
    