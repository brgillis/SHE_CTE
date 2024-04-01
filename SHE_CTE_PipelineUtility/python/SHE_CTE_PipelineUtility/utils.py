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
:file: python/SHE_CTE_PipelineUtility/utils.py

:date: 2024/01/23
:author: @rrollins

"""

from itertools import islice

from astropy.io.fits import HDUList, PrimaryHDU
from SHE_PPT.file_io import save_product_metadata
from ST_DM_FilenameProvider.FilenameProvider import FileNameProvider

from SHE_CTE import SHE_CTE_RELEASE_STRING


def batched(iterable, batch_size):
    """
    Take an input iterator and yield tuples of sequential elements of a given size

    batched('ABCDEFG', 3) --> ('A', 'B', 'C'), ('D', 'E', 'F'), ('G',)

    :param iterable: An input iterable object
    :param batch_size: The number of elements in each batch
    """

    if batch_size < 1:
        raise ValueError('batch_size must be at least one')
    it = iter(iterable)
    while batch := tuple(islice(it, batch_size)):
        yield batch


def ceiling_division(numerator, denominator):
    """
    Returns the smallest integer greater-than or equal-to the ratio of two arguments.

    :param numerator
    :param denominator

    :return value: Integer equivalent to ceil(numerator/denominator) without floating point errors.
    """
    return -((-numerator) // denominator)


def write_data_product(type_name, instance_id, data_product, directory):
    """
    Write an XML data product to an xml file with an automatically generated valid filename.

    :param type_name: Argument passed to get_allowed_filename.
    :param instance_id: Argument passed to get_allowed_filename.
    :param data_product: XML data product object.
    :param directory: Path object for the directory to which the data product will be written.

    :return filename: Valid filename (relative to directory) to which the data product was written.
    """

    filename = FileNameProvider().get_allowed_filename(
        type_name=type_name,
        instance_id=instance_id,
        extension='.xml',
        release=SHE_CTE_RELEASE_STRING,
        processing_function='SHE',
    )
    save_product_metadata(data_product, directory / filename)
    return filename


def write_fits_hdus(type_name, instance_id, hdus, directory):
    """
    Write a tuple of FITS HDUs to a FITS file with an automatically generated valid filename.

    :param type_name: Argument passed to get_allowed_filename.
    :param instance_id: Argument passed to get_allowed_filename.
    :param hdus: Tuple of (non-primary) FITS HDUs.
    :param directory: Path object for the directory to which the FITS file will be written.

    :return filename: Valid filename (relative to directory) to which the FITS file was written.
    """

    filename = FileNameProvider().get_allowed_filename(
        type_name=type_name,
        instance_id=instance_id,
        extension='.fits',
        release=SHE_CTE_RELEASE_STRING,
        processing_function='SHE',
        )
    HDUList([PrimaryHDU(), *hdus]).writeto(directory / filename)
    return filename


def write_fits_table(type_name, instance_id, table, directory):
    """
    Write an astropy Table to a FITS file with an automatically generated valid filename.

    :param type_name: Argument passed to get_allowed_filename.
    :param instance_id: Argument passed to get_allowed_filename.
    :param table: Astropy Table.
    :param directory: Path object for the directory to which the FITS file will be written.

    :return filename: Valid filename (relative to directory) to which the FITS file was written.
    """

    filename = FileNameProvider().get_allowed_filename(
        type_name=type_name,
        instance_id=instance_id,
        extension='.fits',
        release=SHE_CTE_RELEASE_STRING,
        processing_function='SHE',
        )
    table.write(directory / filename)
    return filename
