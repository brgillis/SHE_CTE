""" @file split_fits.py

    Created 25 Aug 2021
    

    Takes Euclid VIS FITS files and splits them into smaller (FITS or HDF5 files), one per CCD.
"""

__updated__ = "2021-08-25"

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


import os
import sys
import json
import time
import datetime

import numpy as np
import matplotlib.pyplot as plt

import h5py

import astropy.io.fits as fits
from astropy.wcs import WCS

from SHE_PPT.logging import getLogger

from .utils import read_config, get_qualified_filename, check_input, create_dummy_output

logger = getLogger(__name__)

#size of the HDF5 chunk. This should probably be 
chunks = (200,200)

#dummy segmentation mask
mask = np.zeros((4136, 4096),np.int64)


def extract_ccd_to_fits(exposure_no, ccd_no, hdul_det, hdul_wgt, hdul_bkg, workdir):
    """Takes various arrays belonging to the specified CCD and writes them to a FITS file"""
    ccd_id = "CCDID "+ccd_no
    
    #Determine the name of the file to write to
    filename = os.path.join(workdir,"CCDs/FITS","%03d_"%(exposure_no)+ccd_no+".fits")

    if os.path.exists(filename):
        logger.warning("FITS file %s exists. Skipping.",filename)
        return filename
    
    #get the HDUs we need from the fits files
    sci = hdul_det[ccd_id+".SCI"]
    sci.header["EXTNAME"] = "SCI"
    flg = hdul_det[ccd_id+".FLG"]
    flg.header["EXTNAME"] = "FLG"
    rms = hdul_det[ccd_id+".RMS"]
    rms.header["EXTNAME"] = "RMS"
    wgt = hdul_wgt[ccd_id]
    wgt.header["EXTNAME"] = "WGT"
    bkg = hdul_bkg[ccd_id]
    bkg.header["EXTNAME"] = "BKG"

    primaryHDU = fits.PrimaryHDU()
    hdul = fits.HDUList([primaryHDU,sci,flg,rms,wgt,bkg])
    
    #create HDU for dummy segmentation mask
    msk = fits.ImageHDU(name="SEG",data=mask)
    hdul.append(msk)

    
    logger.info("Beginning writing CCD to %s",filename)

    hdul.writeto(filename,overwrite=True)

    logger.info("Completed writing CCD to %s",filename)
    
    #make sure to remove the open file handles to the memmapped HDUs
    del(sci.data)
    del(flg.data)
    del(rms.data)
    del(wgt.data)
    del(bkg.data)

    return filename


def _make_dataset(hdf_file, name, hdu, chunks=True, compression=False):
    """Makes a dataset in a HDF5 file"""

    if compression is True:
        compression = "gzip"
    else:
        compression = None
    
    #create the dataset
    logger.info("Begin writing dataset %s", name)
    dataset=hdf_file.create_dataset(name,data=hdu.data,chunks=chunks,compression=compression)
    logger.info("End writing dataset %s", name)

    logger.info("Begin writing dataset %s metadata", name)

    #put the fits header into a "header" attribute
    dataset.attrs["header"] = str(hdu.header)

    #loop over all the items in the header, and set these as atributes of the dataset
    for key, value in hdu.header.items():
        #sometimes we have blank header items... if this is the case, skip them
        if key == "":
            continue
        dataset.attrs[key] = value

    logger.info("End writing dataset %s metadata", name)
    
    #make sure we delete the reference to the data from the HDU
    del(hdu.data)




def extract_ccd_to_hdf5(exposure_no, ccd_no, hdul_det, hdul_wgt, hdul_bkg ,workdir, chunks=True, compression = False):
    """Takes various arrays belonging to the CCD ccd_no and writes them to a HDF5 file"""
    ccd_id = "CCDID "+ccd_no
    
    #determine filename to write this CCD to
    if chunks is not None:
        chunkname="Chunks"
    else:
        chunkname = "NoChunks"
    
    filename = os.path.join(workdir,"CCDs/HDF5","%03d_"%(exposure_no)+ccd_no+"_"+chunkname+".hdf5")

    if os.path.exists(filename):
        logger.warning("HDF5 file %s exists. Skipping.",filename)
        return filename

    #get the HDUs we need from the fits files
    sci = hdul_det[ccd_id+".SCI"]
    flg = hdul_det[ccd_id+".FLG"]
    rms = hdul_det[ccd_id+".RMS"]
    wgt = hdul_wgt[ccd_id]
    bkg = hdul_bkg[ccd_id]


    with h5py.File(filename,"w") as f:

        _make_dataset(f,"SCI",sci,chunks=chunks,compression=compression)
        _make_dataset(f,"FLG",flg,chunks=chunks,compression=compression)
        _make_dataset(f,"RMS",rms,chunks=chunks,compression=compression)
        _make_dataset(f,"WGT",wgt,chunks=chunks,compression=compression)
        _make_dataset(f,"BKG",bkg,chunks=chunks,compression=compression)
        
        #add dummy segmentation mask
        logger.info("Begin writing dummy segmentation map dataset")
        msk = f.create_dataset("SEG", data = mask, chunks=chunks)
        logger.info("End writing dummy segmentation map dataset")
        header = str(fits.ImageHDU(name="SEG",data = mask).header)
        msk.attrs["header"]=header

    logger.info("Written CCD to %s",filename)

    return filename

        

def split_fits_from_args(args):
    """Main entry point. For an exposure, takes its fits files and converts them to smaller fits or HDF5 files, one per CCD"""

    t0 = time.time()
    tstart = datetime.datetime.now()
    logger.info("Start time = %s",tstart)

    workdir = args.workdir
    
    #get relevant config options from the config file
    config = read_config(get_qualified_filename(workdir,args.pipeline_config))

    hdf5 = bool(config["HDF5"])
    chunked = bool(config["chunked"])
    memmap = bool(config["memmap"])
    dry_run = bool(config["dry_run"])
    compression = bool(config["compression"])

    logger.info("Use HDF5 outputs? %s",hdf5)
    if hdf5:
        logger.info("Chunked HDF5? %s",chunked)
        logger.info("HDF5 compression? %s",compression)

    logger.info("Memmapped input FITS? %s",memmap)

    if dry_run:
        #check all the input files exist, create dummy output files, then exit
        logger.info("Dry run:")

        check_input(args.workdir,args.input_fits_json,"input_fits_json")
        create_dummy_output(args.workdir,args.output_json,"output_json")

        return


    if chunked:
        chunks = (200,200)
    else:
        chunks = None

    #get the input file
    input_fits_json = args.input_fits_json

    with open(get_qualified_filename(workdir,input_fits_json), "r") as f:
        input_fits = json.load(f)
    
    #exposure number
    exp = int(input_fits["exp"])
    
    #files for the detector images, weights and backgrounds
    det_file = get_qualified_filename(workdir,input_fits["det"])
    wgt_file = get_qualified_filename(workdir,input_fits["wgt"])
    bkg_file = get_qualified_filename(workdir,input_fits["bkg"])
    
    
    #open these files
    hdul_det = fits.open(det_file,memmap=memmap)
    hdul_wgt = fits.open(wgt_file,memmap=memmap)
    hdul_bkg = fits.open(bkg_file,memmap=memmap)

    #make directories for output files
    try:
        os.mkdir(os.path.join(workdir,"CCDs"))
    except FileExistsError:
        logger.warn("Directory %s already exists",os.path.join(workdir,"CCDs"))
    if hdf5:
        try:
            os.mkdir(os.path.join(workdir,"CCDs","HDF5"))
        except FileExistsError:
            logger.warn("Directory %s already exists",os.path.join(workdir,"CCDs","HDF5"))
    else:
        try:
            os.mkdir(os.path.join(workdir,"CCDs","FITS"))
        except FileExistsError:
            logger.warn("Directory %s already exists",os.path.join(workdir,"CCDs","FITS"))

    
    files_dict = {}
    #loop over all the detectors and create new files for them
    for i in range(1,7):
        for j in range(1,7):

            ccd_no = "%1d-%1d"%(i,j)
            logger.info("Exposure %d, CCD %s:",exp, ccd_no)
            
            #create HDF5 files
            if hdf5:
                filename=extract_ccd_to_hdf5(exp, ccd_no, hdul_det, hdul_wgt, hdul_bkg,workdir,chunks,compression)
            #create FITS files
            else:
                filename=extract_ccd_to_fits(exp, ccd_no, hdul_det, hdul_wgt, hdul_bkg, workdir)

            files_dict[ccd_no] = filename

    hdul_det.close()
    hdul_wgt.close()
    hdul_bkg.close()

    final_dict = {}
    final_dict[exp] = files_dict

    with open(get_qualified_filename(workdir,args.output_json),"w") as f:
        json.dump(final_dict,f,indent=4)

    t1 = time.time()

    d = {
        "tstart": str(tstart),
        "walltime": t1-t0,
        "exposure": exp
    }

    timing_info = get_qualified_filename(args.workdir,args.timing_info)
    logger.info("Writing timing information to %s",timing_info)
    with open(timing_info,"w") as f:
        json.dump(d,f, indent=4)


            
    logger.info("Converted exposure %s in %f seconds",exp,t1-t0)
    logger.info("Completion time = %s",datetime.datetime.now())



    
