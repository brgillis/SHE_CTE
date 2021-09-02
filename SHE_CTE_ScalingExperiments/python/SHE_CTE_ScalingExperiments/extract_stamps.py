""" @file extract_stamps.py

    Created 26 Aug 2021
    

    Extracts postage stamps for a batch.
"""

__updated__ = "2021-08-26"

# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA

import time
import datetime
import os
import json
import sys
from copy import copy

import numpy as np

import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

import h5py


from .utils import read_config, get_qualified_filename, check_input, create_dummy_output
from SHE_PPT.logging import getLogger

logger = getLogger(__name__)





def extract_stamp_hdf5(f, xp, yp):
    """returns 400x400 postage stamps of the various images centred on (xp,yp)"""
    
    ny, nx = f["SCI"].shape

    if xp >= nx or xp < 0 or yp >= ny or yp < 0:
        logger.warn("Object lies outwith this CCD's boundary: (%d, %d)",xp,yp)
    
    #determine coordinates for the maximum and minimum extents of the postage stamp
    i_min = max(xp-200,0)
    i_max = min(xp+200,nx-1)

    j_min = max(yp-200,0)
    j_max = min(yp+200,ny-1)

    #extract postage stamps
    sci = f["SCI"][j_min:j_max, i_min:i_max]
    flg = f["FLG"][j_min:j_max, i_min:i_max]
    rms = f["RMS"][j_min:j_max, i_min:i_max]
    wgt = f["WGT"][j_min:j_max, i_min:i_max]
    bkg = f["BKG"][j_min:j_max, i_min:i_max]
    seg = f["SEG"][j_min:j_max, i_min:i_max]
    
    #if the postage stamps extracted are smaller than 400x400, zero pad so that they are
    if sci.shape[0] != 400 or sci.shape[1] != 400:
        #coordinates of lower left corner of the provided stamp in a 400x400 stamp, S.T. (xc, yc) is in the centre
        xc = 200 - (xp - i_min)
        yc = 200 - (yp - j_min)

        sci, flg, rms, wgt, bkg, seg = _zero_pad(xc, yc, sci, flg, rms, wgt, bkg, seg)
    

    #multiply return values by 2 in order to match extract_stamp_fits
    return sci*2, flg*2, rms*2, wgt*2, bkg*2, seg*2


def extract_stamp_fits(f, xp, yp):
    """returns 400x400 postage stamps of the various images centred on (xp,yp)"""
    
    ny, nx = f["SCI"].data.shape

    if xp >= nx or xp < 0 or yp >= ny or yp < 0:
        logger.warn("Object lies outwith this CCD's boundary (%d, %d)",xp,yp)

    #determine coordinates for the maximum and minimum extents of the postage stamp
    i_min = max(xp-200,0)
    i_max = min(xp+200,nx-1)

    j_min = max(yp-200,0)
    j_max = min(yp+200,ny-1)
   
    #extract postage stamps
    sci = f["SCI"].data[j_min:j_max, i_min:i_max]
    flg = f["FLG"].data[j_min:j_max, i_min:i_max]
    rms = f["RMS"].data[j_min:j_max, i_min:i_max]
    wgt = f["WGT"].data[j_min:j_max, i_min:i_max]
    bkg = f["BKG"].data[j_min:j_max, i_min:i_max]
    seg = f["SEG"].data[j_min:j_max, i_min:i_max]

    #if the postage stamps extracted are smalelr than 400x400, zero pad so that they are
    if sci.shape[0] != 400 or sci.shape[1] != 400:
        #coordinates of lower left corner of the provided stamp in a 400x400 stamp, S.T. (xc, yc) is in the centre
        xc = 200 - (xp - i_min)
        yc = 200 - (yp - j_min)

        sci, flg, rms, wgt, bkg, seg = _zero_pad(xc, yc, sci, flg, rms, wgt, bkg, seg)
    
    #multiply return values by 2 to ensure the arrays are actually fetched from the memmapped array
    return sci*2, flg*2, rms*2, wgt*2, bkg*2, seg*2


def _zero_pad(xc, yc, sci, flg, rms, wgt, bkg, seg):
    """If postage stamp extracted from file is smaller than 400x400 (e.g. near the edge of the CCD),
       then zero pad this stamp such that it is 400x400 """

    ny, nx = sci.shape

    sci_padded = np.zeros((400,400),dtype=np.float32)
    sci_padded[yc:yc+ny,xc:xc+nx] = sci[:,:]

    flg_padded = np.zeros((400,400),dtype=np.int32)
    flg_padded[yc:yc+ny,xc:xc+nx] = flg[:,:]
    
    rms_padded = np.zeros((400,400),dtype=np.float32)
    rms_padded[yc:yc+ny,xc:xc+nx] = rms[:,:]

    wgt_padded = np.zeros((400,400),dtype=np.float32)
    wgt_padded[yc:yc+ny,xc:xc+nx] = wgt[:,:]

    bkg_padded = np.zeros((400,400),dtype=np.float32)
    bkg_padded[yc:yc+ny,xc:xc+nx] = bkg[:,:]

    seg_padded = np.zeros((400,400),dtype=np.int64)
    seg_padded[yc:yc+ny,xc:xc+nx] = seg[:,:]
    

    return sci_padded, flg_padded, rms_padded, wgt_padded, bkg_padded, seg_padded


class CCD():
    """Object for a CCD which contains its wcs and file object"""
    def __init__(self, exposure_number, exposures, ccd_id, hdf5, chunked,workdir):
        """Sets up everything we need, but does not open the file yet, to save having open file handles when not needed"""
        self.exposure_number = int(exposure_number)
        self.ccd_id = ccd_id

        filename = exposures[exposure_number][ccd_id]
        
        qualified_filename = get_qualified_filename(workdir,filename)

        self.qualified_filename = qualified_filename

        self.chunked = chunked
        self.workdir = workdir
        self.hdf5 = hdf5
        
        #Don't open yet..., only open when needed
        self.file = None
        self.wcs = None
        self.is_open = False

    def open_file(self):
        """Opens the file"""
        if self.is_open:
            logger.warn("CCD file %f is already open",self.qualified_filename)
            return

        if  not os.path.exists(self.qualified_filename):
            raise FileNotFoundError('CCD file "%s" not found')

        if self.hdf5:
            if self.chunked:
                self.file = h5py.File(self.qualified_filename,"r",rdcc_nbytes = 200*200*4*60, rdcc_nslots=5987)
            else:
                self.file = h5py.File(self.qualified_filename,"r")

            self.wcs = WCS(self.file["SCI"].attrs["header"])
        else:
            self.file = fits.open(self.qualified_filename)

            self.wcs = WCS(self.file["SCI"].header)

        logger.info("Opened file %s for CCD %s for exposure %s",self.qualified_filename,self.ccd_id,self.exposure_number)

        self.is_open = True

    def __del__(self):
        logger.info("Closing %s",self.qualified_filename)
        self.file.close()

    def __repr__(self):
        return "CCD(%s)"%self.qualified_filename
                


def get_files_for_exposures(batch,exposures, hdf5, chunked,workdir):
    """Returns a dict of a dict of CCDs needed for each exposure
       e.g Dict[exposure] = (Dict[ccd_id] = CCD)
    """
    
    exposures_d={}
    nfiles=0

    for exp in ["1","2","3","4"]:
        files_d={}

        filelist = []
        for obj in batch:
            #get the list of all the CCDs for this exposure (may contain duplicates)
            try:
                filelist.append(obj[exp])
            #The object may not be present in this exposure, id so the key will not be present
            except KeyError:
                pass
        
        #remove duplicates
        filelist = list(set(filelist))
        
        #now create the CCD objects for each file
        for ccd_id in filelist:
            ccd = CCD(exp, exposures, ccd_id, hdf5, chunked,workdir)

            files_d[ccd_id] = ccd
        

        exposures_d[exp] = files_d
        nfiles += len(filelist)

    return exposures_d, nfiles


def extract_stamps_for_batch(batch, files,hdf5, tpause):
    """Extacts stamps for a batch of objects, with a pause between objects to simulate estimating the object's shear"""

    for i in range(len(batch)):
        obj = batch[i]
        ra = obj["ra"]
        dec = obj["dec"]

        logger.info("Extracting stamps for object with (ra, dec) = (%f,%f)",ra,dec)
        t0 = time.time()
        for exp in ["1","2","3","4"]:
            try:
                ccd_id = obj[exp]
            except KeyError:
                #The object is not present in this exposure, move onto the next exposure
                continue
            
            ccd = files[exp][ccd_id]
            
            #open file if not already open
            if not ccd.is_open:
                ccd.open_file()

            x, y = ccd.wcs.world_to_pixel(SkyCoord(ra, dec, frame="icrs", unit="deg"))

            if hdf5:
                extract_stamp_hdf5(ccd.file,int(x),int(y))
            else:
                extract_stamp_fits(ccd.file,int(x),int(y))
        t1 = time.time()
        logger.info("Stamp extraction took %fs",t1-t0)
        
        #sleep for some time to simulate doing some compute
        logger.info("Sleeping for %2.1fs to simulate shear estimation compute time",tpause[i])
        time.sleep(tpause[i])



def extract_stamps_from_args(args):

    t0 = time.time()
    tstart = datetime.datetime.now()
    logger.info("Started at %s",tstart)

    #get pipeline_config
    config = read_config(get_qualified_filename(args.workdir,args.pipeline_config))

    hdf5 = bool(config["HDF5"])
    chunked = bool(config["chunked"])
    dry_run = bool(config["dry_run"])
    mean_compute_time = config["mean_compute_time"]

    if dry_run:
        #check all the input files exist, create dummy output files, then exit
        logger.info("Dry run:")

        check_input(args.workdir,args.batch_file,"batch_file")
        check_input(args.workdir,args.exposures,"exposures")
        create_dummy_output(args.workdir,args.timing_info,"timing_info")

        return

    #load the exposures
    exposures_json = get_qualified_filename(args.workdir,args.exposures)
    with open(exposures_json,"r") as f:
        exposures = json.load(f)

    #load the batch
    batch_json = get_qualified_filename(args.workdir,args.batch_file)
    with open(batch_json,"r") as f:
        batch = json.load(f)

    #determine batch number, and seed the RNG with this
    fname = os.path.split(args.batch_file)[-1]
    num = int(fname.split("_")[1].split(".")[0])
    rng = np.random.RandomState(num)

    logger.info("Seeding RNG with seed %d",num)
    
    #determine the time each object should pause for to simulate compute time
    #pauses = rng.noncentral_chisquare(2,mean_compute_time,len(batch))
    pauses = rng.rayleigh(mean_compute_time,len(batch))
    tpause = np.sum(pauses)

    #get file descriptors for each file we need to read from
    files, num_files = get_files_for_exposures(batch,exposures,hdf5,chunked,args.workdir)
    
    #extract the stamps
    extract_stamps_for_batch(batch,files,hdf5,pauses)

    num_objects = len(batch)

    walltime = time.time()-t0

    timing_info = {
        "num_objects": num_objects,
        "num_files": num_files,
        "walltime": walltime,
        "tstart": str(tstart),
        "compute_time": tpause
    }

    outfile = get_qualified_filename(args.workdir,args.timing_info)
    logger.info("Writing timing info for this batch to %s",outfile)
    with open(outfile,"w") as f:
        json.dump(timing_info,f,indent=4)


   

    t1=time.time()
    logger.info("Completion time = %s",datetime.datetime.now())
    logger.info("Total walltime = %fs",t1-t0)