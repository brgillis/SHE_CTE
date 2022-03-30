""" @file executor.py

    Created 22 Feb 2022

    Class to handle primary execution of SHE_CTE executables.
"""

__updated__ = "2022-02-22"

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

from dataclasses import dataclass

from SHE_PPT.executor import LogOptions, SheExecutor
from . import __version__


@dataclass
class CteLogOptions(LogOptions):
    """ Subclass of LogOptions which overrides defaults for project_name and project_version.
    """
    project_name: str = "SHE_CTE"
    project_version: str = __version__


class SheCteExecutor(SheExecutor):
    """ Subclass of SheExecutor which overrides attribute types.
    """

    # Attributes with different types from base class
    log_options: CteLogOptions
