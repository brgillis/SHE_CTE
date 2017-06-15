""" @file estimate_shears.py

    Created 27 Mar 2017

    Primary execution loop for measuring galaxy shapes from an image file.

    ---------------------------------------------------------------------

    Copyright (C) 2017 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from astropy.io import fits
from astropy.table import Table
import galsim

from icebrgpy.logging import getLogger

from SHE_CTE_ShearEstimation import magic_values as mv
from SHE_SIM_galaxy_image_generation import magic_values as sim_mv

from SHE_CTE_ShearEstimation.estimate_shear import estimate_shear
from SHE_CTE_ShearEstimation.extract_stamps import extract_stamps
from SHE_CTE_ShearEstimation.output_shear_estimates import output_shear_estimates


def find_value(args_value, name, label, detections_table, galaxies_hdulist):
    if args_value is not None:
        value = args_value
    else:
        try:
            value = galaxies_hdulist[0].header[label]
        except KeyError as _e1:
            try:
                value = detections_table.meta[label]
            except KeyError as _e2:
                raise KeyError("No " + name + " value available.")
    return value

def estimate_shears_from_args(kwargs):
    """
    @brief
        Perform shear estimation, given arguments from the command-line.
    
    @param kwargs <dict>
    
    @return None
    """

    logger = getLogger(mv.logger_name)
    
    logger.debug("Entering estimate_shears_from_args")
    
    # Load the detections table
    detections_table = Table.read(kwargs["detections_table_file_name"])
    
    # Load the galaxies and psf images
    galaxies_hdulist = fits.open(kwargs["galaxies_image_file_name"])
    psf_hdulist = fits.open(kwargs["psf_image_file_name"])
    
    # Load the P(e) table if available
    p_of_e_table_file_name = kwargs["p_of_e_table_file_name"]
    if p_of_e_table_file_name is not None:
        p_of_e_table = Table.read(p_of_e_table_file_name)
        e_half_step = (p_of_e_table["E_LOW"][1] - p_of_e_table["E_LOW"][0])/2.
        shape_noise_var = (((p_of_e_table["E_LOW"]+e_half_step)**2 * p_of_e_table["E_COUNT"]).sum() /
                           p_of_e_table["E_COUNT"].sum())
    else:
        shape_noise_var = 0.06
    
    # Get the gain, subtracted sky level, and read noise from the galaxies image
    # if they weren't passed at the command-line
    gain = find_value(kwargs["gain"], "gain", sim_mv.fits_header_gain_label, 
                      detections_table, galaxies_hdulist)
    subtracted_sky_level = find_value(kwargs["subtracted_sky_level"], "subtracted_sky_level",
                                      sim_mv.fits_header_subtracted_sky_level_label, 
                                      detections_table, galaxies_hdulist)
    read_noise = find_value(kwargs["read_noise"], "read_noise", sim_mv.fits_header_read_noise_label, 
                            detections_table, galaxies_hdulist)
    
    # Get a list of galaxy postage stamps
    gal_stamps = extract_stamps(detections_table,
                                galaxies_hdulist,
                                sim_mv.detections_table_gal_xc_label,
                                sim_mv.detections_table_gal_yc_label,)
    
    # Get a list of PSF postage stamps
    psf_stamps = extract_stamps(detections_table,
                                psf_hdulist,
                                sim_mv.detections_table_psf_xc_label,
                                sim_mv.detections_table_psf_yc_label,)
    
    # Estimate the shear for each stamp
    for i, gal_stamp, psf_stamp in zip(range(len(gal_stamps)), gal_stamps, psf_stamps):
        if i % 10 == 0:
            logger.info("Measuring shear for galaxy " + str(i) + "/" + str(len(gal_stamps)) + ".")
        shear_estimate = estimate_shear(galaxy_image=gal_stamp.image,
                                        psf_image=psf_stamp.image,
                                        method=kwargs["method"],
                                        gain=gain,
                                        subtracted_sky_level=subtracted_sky_level,
                                        read_noise=read_noise,
                                        shape_noise_var=shape_noise_var)
        gal_stamp.shear_estimate = shear_estimate
        
    logger.info("Finished estimating shear. Outputting results to " + kwargs["output_file_name"] + ".")
        
    # Output the measurements
    output_shear_estimates(stamps=gal_stamps,
                           output_file_name=kwargs["output_file_name"], 
                           galaxies_image_file_name=kwargs["galaxies_image_file_name"])
    
    logger.debug("Exiting estimate_shears_from_args")
    
    return