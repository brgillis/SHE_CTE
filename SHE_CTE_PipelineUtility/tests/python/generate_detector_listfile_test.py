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


@pytest.mark.parametrize("instrument,detectors", [('NISP', NISP_CCD_IDS), ('VIS', VIS_CCD_IDS)])
def test_generate_detector_listfile(tmp_path_factory, instrument, detectors):
    parser = defineSpecificProgramOptions()
    arg_string = [
        '--workdir', str(tmp_path_factory.mktemp('workdir')),
        '--logdir', str(tmp_path_factory.mktemp('logdir')),
        '--instrument', instrument,
        '--detectors', 'listfile.json',
        ]
    args = parser.parse_args(arg_string)
    mainMethod(args)
    detector_listfile = Path(args.workdir, args.detectors)
    assert detector_listfile.is_file()
    with open(detector_listfile) as filestream:
        she_cte_ccd_ids = json.load(filestream)
    assert sorted(she_cte_ccd_ids) == sorted(detectors)
