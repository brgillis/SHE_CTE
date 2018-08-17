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
File: python/SHE_MomentsML/play.py

Created on: 08/25/17
"""

import argparse
import ElementsKernel.Logging as log

from ElementsKernel.ProjectCommonRoutines import getAuxPathFile

import numpy as np
import SHE_PPT.she_image
#import SHE_PPT.vis_helper
import SHE_PPT.she_image_checkplot
import SHE_MomentsML.params
import SHE_MomentsML.estimate_shear
import SHE_MomentsML.utils_table
import SHE_MomentsML.utils_io
import SHE_CalibMLCore
import astropy
import galsim
import matplotlib
import os

from SHE_PPT.table_formats import detections as detf

logger = log.getLogger('play')
log.setLevel("DEBUG")

def defineSpecificProgramOptions():
    """
    """
    parser = argparse.ArgumentParser()
    return parser


def version():
    logger.info("Numpy {}".format(np.__version__))
    logger.info("Astropy {}".format(astropy.__version__))
    logger.info("Galsim {}".format(galsim.__version__))
    logger.info("Matplotlib {}".format(matplotlib.__version__))
 
 
def play_ml():
    l = SHE_CalibMLCore.layer.Layer(2, 2)
    print(l)
    
    
def play_img():
    img = SHE_PPT.she_image.SHEImage(np.random.randn(100).reshape((10,10)))
    print(img)
    img.offset = (3, 4)
    print(img)
    img.write_to_fits("testoffset.fits", clobber=True)
    im2 = SHE_PPT.she_image.SHEImage.read_from_fits("testoffset.fits")
    print(img.offset)
    print(im2.offset)
    

def play_run():
    
    stack = SHE_PPT.vis_helper.read_shestack_from_gst_json("/home/user/Data/SHE_SIM/quick_run/output_files.json", ccdid="0")
    
    img = stack.exposures[0].science_image
    print(img)
    img = img.extract_stamp(128, 128, 256)
    checkplot = SHE_PPT.she_image_checkplot.Checkplot(img, scale=1)
    checkplot.show()


def play_params():
    
    #SHE_MomentsML.params.MomentsMLParams.make_pickle_from_dir("/home/user/Work/Projects/SHE_MomentsML/demoparams", "test.pkl")
    
    #params = SHE_MomentsML.params.MomentsMLParams.read_from_dir("/home/user/Work/Projects/SHE_MomentsML/demoparams")

    #params = SHE_MomentsML.params.MomentsMLParams.read_pickle("test.pkl")
    paramsdir = getAuxPathFile("SHE_MomentsML/default_params")
    print(paramsdir)
    
    params = SHE_MomentsML.params.MomentsMLParams.read_dir(paramsdir)
    
    print(params)



def play_test_frame_stack():
   
    frame_stack_path = getAuxPathFile('SHE_PPT/test_she_frame_stack_simple.bin')
    frame_stack = SHE_MomentsML.utils_io.read_pickle(frame_stack_path)
    print(frame_stack)
    
    """
    for row in frame_stack.detections_catalogue:
        stamp_stack = frame_stack.extract_galaxy_stack(row[detf.tf.ID], width=100)

        for exposure in stamp_stack.exposures:
            print(exposure) # This is a SHEFrame
            checkplot = SHE_PPT.she_image_checkplot.Checkplot(exposure, scale=1, z1="auto", z2="auto")
            checkplot.show()
    """
   
    output_cat = SHE_MomentsML.estimate_shear.estimate_shear(frame_stack, training_data=None, calibration_data=None, workdir=".")
    print(SHE_MomentsML.utils_table.get_info(output_cat))
    print(output_cat)
    


    
def play_estimate_shear():
    """
    """
    #stack = SHE_PPT.vis_helper.read_shestack_from_gst_json("/home/user/Data/SHE_SIM/quick_run/output_files.json", ccdid="0")
    #method_data = SHE_MomentsML.params.MomentsMLParams.read_pickle("test.pkl")
    #output_cat = SHE_MomentsML.estimate_shear.estimate_shear(stack, method_data)
    
    #output_cat = SHE_MomentsML.estimate_shear.estimate_shear(stack)

    #print(SHE_MomentsML.utils_table.get_info(output_cat))
    #print(output_cat)
    
    
    """
    img = stack.exposures[0].science_image
    checkplot = SHE_PPT.she_image_checkplot.Checkplot(img, scale=1, z1="auto", z2="auto")
    checkplot.draw_g_ellipses(output_cat, x="adamom_x", y="adamom_y", g1="adamom_g1", g2="adamom_g2", sigma="adamom_sigma")
    output_cat["s_size"] = 2.0 * output_cat["adamom_sigma"]
    checkplot.draw_g_ellipses(output_cat, x="adamom_x", y="adamom_y", g1="pre_s1", g2="pre_s2", sigma="s_size", color="green")
    checkplot.show()
    """
    
    
    

def mainMethod(args):
    """
    """

    logger = log.getLogger('play')

    version()
    
    play_test_frame_stack()
    #play_img()
    
    #play_params()
    #play_run()
    #play_estimate_shear()
    
   
