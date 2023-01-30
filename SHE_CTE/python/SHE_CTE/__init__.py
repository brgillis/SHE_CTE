""" @file __init__.py

    Created 22 April 2019

    SHE_CTE package, for modules general to SHE_CTE
"""

__updated__ = "2021-08-19"

# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA


# Get the version from the compiled file created by Elements
from SHE_CTE_VERSION import SHE_CTE_VERSION_STRING

__version__ = SHE_CTE_VERSION_STRING

import re
SHE_CTE_RELEASE_STRING = re.match(r"[0-9]{1,2}\.[0-9]{1,2}", SHE_CTE_VERSION_STRING).group()

from pkgutil import extend_path

__path__ = extend_path(__path__, __name__)
