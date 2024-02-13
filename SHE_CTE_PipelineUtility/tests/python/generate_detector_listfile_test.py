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
import pytest

from SHE_CTE_PipelineUtility.GenerateDetectorListfile import defineSpecificProgramOptions, mainMethod

NISP_CCD_IDS = [
    '1-1', '1-2', '1-3', '1-4',
    '2-1', '2-2', '2-3', '2-4',
    '3-1', '3-2', '3-3', '3-4',
    '4-1', '4-2', '4-3', '4-4',
    ]

VIS_CCD_IDS = [
    '1-1', '1-2', '1-3', '1-4', '1-5', '1-6',
    '2-1', '2-2', '2-3', '2-4', '2-5', '2-6',
    '3-1', '3-2', '3-3', '3-4', '3-5', '3-6',
    '4-1', '4-2', '4-3', '4-4', '4-5', '4-6',
    '5-1', '5-2', '5-3', '5-4', '5-5', '5-6',
    '6-1', '6-2', '6-3', '6-4', '6-5', '6-6',
    ]

VIS_QUADRANT_IDS = [
    '1-1.E', '1-1.F', '1-1.G', '1-1.H', '1-2.E', '1-2.F', '1-2.G', '1-2.H',
    '1-3.E', '1-3.F', '1-3.G', '1-3.H', '1-4.E', '1-4.F', '1-4.G', '1-4.H',
    '1-5.E', '1-5.F', '1-5.G', '1-5.H', '1-6.E', '1-6.F', '1-6.G', '1-6.H',
    '2-1.E', '2-1.F', '2-1.G', '2-1.H', '2-2.E', '2-2.F', '2-2.G', '2-2.H',
    '2-3.E', '2-3.F', '2-3.G', '2-3.H', '2-4.E', '2-4.F', '2-4.G', '2-4.H',
    '2-5.E', '2-5.F', '2-5.G', '2-5.H', '2-6.E', '2-6.F', '2-6.G', '2-6.H',
    '3-1.E', '3-1.F', '3-1.G', '3-1.H', '3-2.E', '3-2.F', '3-2.G', '3-2.H',
    '3-3.E', '3-3.F', '3-3.G', '3-3.H', '3-4.E', '3-4.F', '3-4.G', '3-4.H',
    '3-5.E', '3-5.F', '3-5.G', '3-5.H', '3-6.E', '3-6.F', '3-6.G', '3-6.H',
    '4-1.E', '4-1.F', '4-1.G', '4-1.H', '4-2.E', '4-2.F', '4-2.G', '4-2.H',
    '4-3.E', '4-3.F', '4-3.G', '4-3.H', '4-4.E', '4-4.F', '4-4.G', '4-4.H',
    '4-5.E', '4-5.F', '4-5.G', '4-5.H', '4-6.E', '4-6.F', '4-6.G', '4-6.H', 
    '5-1.E', '5-1.F', '5-1.G', '5-1.H', '5-2.E', '5-2.F', '5-2.G', '5-2.H',
    '5-3.E', '5-3.F', '5-3.G', '5-3.H', '5-4.E', '5-4.F', '5-4.G', '5-4.H',
    '5-5.E', '5-5.F', '5-5.G', '5-5.H', '5-6.E', '5-6.F', '5-6.G', '5-6.H',
    '6-1.E', '6-1.F', '6-1.G', '6-1.H', '6-2.E', '6-2.F', '6-2.G', '6-2.H',
    '6-3.E', '6-3.F', '6-3.G', '6-3.H', '6-4.E', '6-4.F', '6-4.G', '6-4.H',
    '6-5.E', '6-5.F', '6-5.G', '6-5.H', '6-6.E', '6-6.F', '6-6.G', '6-6.H'
]


@pytest.mark.parametrize("instrument,detectors", [('NISP', NISP_CCD_IDS), ('VIS', VIS_CCD_IDS), ("VIS", VIS_QUADRANT_IDS)])
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
