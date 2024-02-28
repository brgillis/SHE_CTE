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
:file: python/SHE_CTE_PipelineUtility/FlagEmptyBatches.py

:date: 2024/02/05
:author: @rrollins

"""

from argparse import ArgumentParser, ArgumentTypeError
from pathlib import Path

from ElementsKernel.Logging import getLogger
from SHE_PPT.file_io import read_listfile, read_xml_product, write_listfile


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
    parser.add_argument('--mer_catalog_listfile', type=str, help='Listfile MER catalog data products')
    parser.add_argument('--flags', type=str, help='Listfile of flags for empty batches')

    return parser


def mainMethod(args):
    """
    @brief The "main" method.

    @details This method is the entry point to the program. In this sense, it is
    similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger('SHE_CTE_FlagEmptyBatches')

    logger.info('#')
    logger.info('# Entering SHE_CTE_FlagEmptyBatches mainMethod()')
    logger.info('#')

    mer_files = read_listfile(args.mer_catalog_listfile, workdir=args.workdir)
    batch_flags = [bool(read_xml_product(f, args.workdir).Data.QualityParams.ObjectCount.value()) for f in mer_files]
    write_listfile(args.flags, batch_flags, workdir=args.workdir)

    logger.info('#')
    logger.info('# Exiting SHE_CTE_FlagEmptyBatches mainMethod()')
    logger.info('#')

    return


if __name__ == '__main__':
    parser = defineSpecificProgramOptions()
    args = parser.parse_args()
    mainMethod(args)
