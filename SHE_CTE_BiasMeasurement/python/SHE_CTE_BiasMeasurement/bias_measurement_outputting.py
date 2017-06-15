""" @file bias_measurement_outputting.py

    Created 10 Apr 2017

    Output bias measurements in an appropriate table

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

from astropy.table import Table

def output_bias_measurement(bias_measurement, output_file_name, output_format):
    
    # Set up a Table which we'll write out
    output_table = Table(names=["dimension","m","m_err","c","c_err","mc_covar"],
                         dtype=[int,float,float,float,float,float])
    
    output_table.add_row({"dimension":0,
                          "m":bias_measurement.get_m(),
                          "m_err":bias_measurement.get_m_err(),
                          "c":bias_measurement.get_c(),
                          "c_err":bias_measurement.get_c_err(),
                          "mc_covar":bias_measurement.get_mc_covar()})
    output_table.add_row({"dimension":1,
                          "m":bias_measurement.m1,
                          "m_err":bias_measurement.m1_err,
                          "c":bias_measurement.c1,
                          "c_err":bias_measurement.c1_err,
                          "mc_covar":bias_measurement.m1c1_covar})
    output_table.add_row({"dimension":2,
                          "m":bias_measurement.m2,
                          "m_err":bias_measurement.m2_err,
                          "c":bias_measurement.c2,
                          "c_err":bias_measurement.c2_err,
                          "mc_covar":bias_measurement.m2c2_covar})
    
    output_table.write(output_file_name, format=output_format, overwrite=True)