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
:file: tests/python/generate_detector_listfile.py

:date: 2023/03/18
:author: Richard Rollins

"""

import json
from pathlib import Path
from itertools import product
import pytest

from SHE_CTE_PipelineUtility.GenerateDetectorListfile import defineSpecificProgramOptions, mainMethod

NISP_CCD_IDS = [f"{y+1}-{x+1}" for y, x in product(range(4), range(4))]

VIS_CCD_IDS = [f"{y+1}-{x+1}" for y, x in product(range(6), range(6))]

VIS_QUADRANT_IDS = [f"{y+1}-{x+1}.{q}" for y, x, q in product(range(6), range(6), ("E", "F", "G", "H"))]


@pytest.mark.parametrize(
    "instrument,detectors",
    [('NISP', NISP_CCD_IDS), ('VIS', VIS_CCD_IDS), ("VIS", VIS_QUADRANT_IDS)]
)
def test_generate_detector_listfile(tmp_path_factory, instrument, detectors):
    parser = defineSpecificProgramOptions()
    arg_string = [
        '--workdir', str(tmp_path_factory.mktemp('workdir')),
        '--logdir', str(tmp_path_factory.mktemp('logdir')),
        '--instrument', instrument,
        '--detectors', 'listfile.json',
        ]
    if len(detectors) == 144:
        arg_string.append("--quadrants")
    args = parser.parse_args(arg_string)
    mainMethod(args)
    detector_listfile = Path(args.workdir, args.detectors)
    assert detector_listfile.is_file()
    with open(detector_listfile) as filestream:
        she_cte_ccd_ids = json.load(filestream)
    assert sorted(she_cte_ccd_ids) == sorted(detectors)
