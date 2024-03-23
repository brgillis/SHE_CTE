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
:file: python/SHE_CTE_PipelineUtility/SHE_CTE_GenerateDetectorListfile.py

:date: 2023/03/18
:author: Richard Rollins

Executable that writes a listfile of all VIS or NISP CCD IDs. Used as a
pipeline split task for processing individual detectors in parallel.

"""

from argparse import ArgumentParser, ArgumentTypeError
import itertools
import json
from pathlib import Path

from ElementsKernel.Logging import getLogger
from EL_CoordsUtils.telescope_coords import nisp_det_specs, vis_det_specs

from SHE_PPT.argument_parser import dir_path

instrument_specs = {
    'NISP': nisp_det_specs,
    'VIS': vis_det_specs,
    }

QUADRANT_NAMES = ("E", "F", "G", "H", )

def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    parser = ArgumentParser()
    parser.add_argument('--workdir', type=dir_path, help='Work directory')
    parser.add_argument('--logdir', type=dir_path, help='Log directory')
    parser.add_argument('--instrument', type=str, choices=instrument_specs.keys(), help='Name of instrument')
    parser.add_argument('--detectors', type=str, help='Output JSON listfile of CCDIDs')
    parser.add_argument('--quadrants', action="store_true", help="Create VIS quadrants rather than CCDs")

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger('SHE_CTE_GenerateDetectorListfile')

    logger.info('#')
    logger.info('# Entering SHE_CTE_GenerateDetectorListfile mainMethod()')
    logger.info('#')

    ndet_x = instrument_specs[args.instrument].ndet_x
    ndet_y = instrument_specs[args.instrument].ndet_y
    if args.quadrants and args.instrument == "VIS":
        ccd_ids = [f'{1+y}-{1+x}.{q}' for x, y, q in itertools.product(range(ndet_x), range(ndet_y), QUADRANT_NAMES)]
    else:
        ccd_ids = [f'{1+y}-{1+x}' for x, y in itertools.product(range(ndet_x), range(ndet_y))]
    filename = Path(args.workdir, args.detectors)
    logger.info('# Writing %s %s CCDIDs to %s', len(ccd_ids), args.instrument, str(filename))
    with open(filename, 'w') as filestream:
        json.dump(ccd_ids, filestream)

    logger.info('#')
    logger.info('# Exiting SHE_CTE_GenerateDetectorListfile mainMethod()')
    logger.info('#')
