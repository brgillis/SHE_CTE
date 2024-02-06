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
:file: python/SHE_CTE_PipelineUtility/FilterListfile.py

:date: 2024/02/05
:author: @rrollins

"""

from argparse import ArgumentParser, ArgumentTypeError
from pathlib import Path

from ElementsKernel.Logging import getLogger
from SHE_PPT.file_io import read_listfile, write_listfile


def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program

    @details See the Elements documentation for more details.
    @return An  ArgumentParser.
    """

    def dir_path(path):
        if not Path(path).is_dir():
            raise ArgumentTypeError(f"{path} is not a valid path")
        return path

    parser = ArgumentParser()
    parser.add_argument('--workdir', type=dir_path, help='Name of the working directory')
    parser.add_argument('--logdir', type=str, help='Log directory')
    parser.add_argument('--input_listfile', type=str, help='Listfile SHE Object ID list data products')
    parser.add_argument('--output_listfile', type=str, help='Listfile SHE Object ID list data products')
    parser.add_argument('--flags', type=str, help='Listfile of flags for empty batches')

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger('SHE_CTE_FilterListfile')

    logger.info('#')
    logger.info('# Entering SHE_CTE_FilterListfile mainMethod()')
    logger.info('#')

    inptut = read_listfile(args.input_listfile, workdir=args.workdir)
    flags = read_listfile(args.flags, workdir=args.workdir)
    output = [i for i, f in zip(inptut, flags) if f]
    write_listfile(args.output_listfile, output, workdir=args.workdir)

    logger.info('#')
    logger.info('# Exiting SHE_CTE_FilterListfile mainMethod()')
    logger.info('#')

    return


if __name__ == '__main__':
    parser = defineSpecificProgramOptions()
    args = parser.parse_args()
    mainMethod(args)
