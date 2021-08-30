""" @file extract_objects.py

    Created 25 Aug 2021
    

    Extracts objects from the MER table(s) for our observation.
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
from copy import deepcopy
import json
import time
import datetime

from multiprocessing import Pool
from functools import partial


import numpy as np

import astropy.io.fits as fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import h5py


from SHE_PPT.logging import getLogger

from .utils import read_config, get_qualified_filename, check_input, create_dummy_output


logger=getLogger(__name__)

NUMPROCS = 8



def get_stacked_pixel_coords(localWCS, stackedWCS, local_x, local_y):
    """Given Pixel coordinates in CCD pixels, returns the coordinates on the stacked frame"""
    #convert to world coordinates
    world = localWCS.pixel_to_world(local_x,local_y)
    #convert to stacked frame coordinates
    xp, yp = stackedWCS.world_to_pixel(world)
    return float(xp), float(yp)


def get_catalogue(catalogue_list):
    """concatenates all the MER catalogues into one"""

    tlist=[]
    
    #read in each catalogue, appending them to the list tlist
    for f in catalogue_list:
        logger.info("Reading %s",f)
        t = deepcopy(Table.read(f))
        tlist.append(t)
    
    #merge all the catalogues into one
    catalogue = vstack(tlist)

    #make sure remove close all the table file objects
    del(tlist)

    return catalogue


def get_ccd_wcs(filename, workdir, hdf5):
    """Gets the WCS for a ccd from a given exposure"""
    
    qualified_filename = get_qualified_filename(workdir,filename)

    logger.info("Getting WCS for %s",qualified_filename)

    if hdf5:
        with h5py.File(qualified_filename) as f:
            header = f["SCI"].attrs["header"]
            ny, nx = f["SCI"].shape
    else:
        with fits.open(qualified_filename) as f:
            header = f["SCI"].header
            ny, nx = f["SCI"].data.shape
            #make sure we close the file pointer to this
            del(f["SCI"].data)

    return WCS(header), nx, ny



def get_stacked_wcs(stacked_file):
    """Gets the WCS for the stacked frame"""
    logger.info("Getting WCS for stacked file %s",stacked_file)
    with fits.open(stacked_file) as f:
        header = f["SCI"].header
    
    return WCS(header)


def get_ccd_boundaries(CCDs,workdir,stackedWCS,hdf5):
    """Returns a dictonary that for each exposure and CCD contains the coordinates of the lower left and upper right pixels of the CCD in the stacked image pixel coodinates"""
    
    #return dictionary
    d = {}

    for exp in ["1","2","3","4"]:
        #dictionary for this exposure
        dexp={}

        #for each ccd, extract its coordinates in the stacked imge
        for ccd_name, filename in CCDs[exp].items():
            wcs, nx, ny =get_ccd_wcs(filename,workdir,hdf5=hdf5)

            lowerx, lowery = get_stacked_pixel_coords(wcs,stackedWCS,0,0)
            upperx, uppery = get_stacked_pixel_coords(wcs,stackedWCS,nx-1,ny-1)

            dexp[ccd_name]={}
            dexp[ccd_name]["lower"] = (lowerx, lowery)
            dexp[ccd_name]["upper"] = (upperx, uppery)
                

        d[exp] = dexp

    return d

def get_object_ccds(x,y,ccdBoundaries):
    """For an object's stacked image x and y pixel coordinates, determine which CCD it is in for each exposure"""
    
    #dictionary that will contain {exposure_number: ccd_id}
    d={}
    
    #loop over all 4 exposures
    for exp in ["1","2","3","4"]:
        #for each CCD, see if the object is in it
        for ccd_name, ccd_info in ccdBoundaries[exp].items():
            xmin, ymin = ccd_info["lower"]
            xmax, ymax = ccd_info["upper"]

            if x > xmin and x < xmax and y > ymin and y < ymax:
                #this is the CCD that the object is in
                d[exp] = ccd_name
                break

    return d



def get_objects_in_fov_from_catalogue(catalogue, stackedWCS, ccdBoundaries):
    """Goes through all the objects in the catalogue, and returns a list of those in the FOV, where the list contains information on which CCD and exposure each can be found on """

    print("Getting list of objects in the FOV and their exposures/CCDs")
    
    # Offsets of each exposure are 
    #
    #    Exp1
    #      Exp2
    #        Exp3
    #          Exp4
    #
    
    # work out the maximum pixel range on the stacked image where there is an image.
    # We don't care for objects outwith this as they are not imaged in any CCD
    xmin, _ = ccdBoundaries["1"]["1-1"]["lower"]
    xmax, _ = ccdBoundaries["4"]["6-6"]["upper"]
    _, ymin = ccdBoundaries["4"]["1-1"]["lower"]
    _, ymax = ccdBoundaries["1"]["6-6"]["upper"]

    objects = []

    for i in range(len(catalogue)):
        index = i
        row = catalogue[i]
        ra = row["RIGHT_ASCENSION"]
        dec = row["DECLINATION"]

        sky = SkyCoord(ra, dec, frame="icrs", unit="deg")
        x, y = stackedWCS.world_to_pixel(sky)

        if x > xmin and x < xmax:
            if y > ymin and y < ymax:
                #object is in an exposure

                obj = get_object_ccds(x, y, ccdBoundaries)
                if len(obj) == 0:
                    #object is not in any CCDs
                    continue
                obj["ra"] = ra
                obj["dec"] = dec
                obj["index"] = index
                obj["x"] = int(x)
                obj["y"] = int(y)
                #print(obj)
                objects.append(obj)
    
    logger.info("Exctacted %d objects",len(objects))
    
    return objects



def get_objects_in_fov(all_objects, stackedWCS, ccdBoundaries):
    """Goes through all the objects in the catalogue, and returns a list of those in the FOV, where the list contains information on which CCD and exposure each can be found on """

    # Offsets of each exposure are 
    #
    #    Exp1
    #      Exp2
    #        Exp3
    #          Exp4
    #
    
    # work out the maximum pixel range on the stacked image where there is an image.
    # We don't care for objects outwith this as they are not imaged in any CCD
    xmin, _ = ccdBoundaries["1"]["1-1"]["lower"]
    xmax, _ = ccdBoundaries["4"]["6-6"]["upper"]
    _, ymin = ccdBoundaries["4"]["1-1"]["lower"]
    _, ymax = ccdBoundaries["1"]["6-6"]["upper"]

    objects = []

    for item in all_objects:
        index = item["index"]
        ra = item["ra"]
        dec = item["dec"]

        sky = SkyCoord(ra, dec, frame="icrs", unit="deg")
        x, y = stackedWCS.world_to_pixel(sky)

        if x > xmin and x < xmax:
            if y > ymin and y < ymax:
                #object is in an exposure

                obj = get_object_ccds(x, y, ccdBoundaries)
                if len(obj) == 0:
                    #object is not in any CCDs
                    continue
                obj["ra"] = ra
                obj["dec"] = dec
                obj["index"] = index
                obj["x"] = int(x)
                obj["y"] = int(y)
                #print(obj)
                objects.append(obj)

    pid = os.getpid()
    
    logger.info("(PID: %d) Exctacted %d objects for batch",pid,len(objects))
    
    return objects


def extract_objects_from_args(args):

    t0 = time.time()
    logger.info("Started at %s",datetime.datetime.now())

    #extract the options we need from the pipeline_config
    config = read_config(get_qualified_filename(args.workdir,args.pipeline_config))
    hdf5 = bool(config["HDF5"])
    dry_run = bool(config["dry_run"])

    if dry_run:
        #check all the input files exist, create dummy output files, then exit
        logger.info("Dry run:")

        check_input(args.workdir,args.catalogue_listfile,"catalogue_listfile")
        check_input(args.workdir,args.stacked_image,"stacked_image")
        check_input(args.workdir,args.exposures,"exposures")

        create_dummy_output(args.workdir,args.output_objects_list,"output_objects_list")
        create_dummy_output(args.workdir,args.combined_catalogue,"combined_catalogue")

        return
    
    #get the input files and extract their contents where necessary
    catalogue_listfile = get_qualified_filename(args.workdir,args.catalogue_listfile)
    logger.info("Catalogue listfile = %s",catalogue_listfile)
    with open(catalogue_listfile,"r") as f:
        catalogue_list = json.load(f)
    
    exposures = get_qualified_filename(args.workdir,args.exposures)
    logger.info("Exposures = %s",exposures)
    with open(exposures,"r") as f:
        CCDs= json.load(f)

    stacked_image = get_qualified_filename(args.workdir,args.stacked_image)
    logger.info("Stacked_image = %s",stacked_image)


    catalogue = get_catalogue(catalogue_list)

    stackedWCS = get_stacked_wcs(stacked_image)

    ccdBoundaries = get_ccd_boundaries(CCDs,args.workdir,stackedWCS,hdf5)
    
    #create list of object ra, dec and their row number within the catalogue
    all_objects = []
    for index, row in enumerate(catalogue):
        d={}
        d["index"] = index
        d["ra"] = row["RIGHT_ASCENSION"]
        d["dec"] = row["DECLINATION"]
        all_objects.append(d)
    
    #now split these into batches to be processed in parallel
    #  (5 batches per process to compromise load balance and too many batches)
    numbatches = 5*NUMPROCS
    batchsize = len(all_objects)//numbatches
    
    batches=[]
    for i in range(numbatches):
        batch = all_objects[i*batchsize:(i+1)*batchsize]
        batches.append(batch)

    lastbatch = all_objects[numbatches*batchsize:]
    if len(lastbatch) > 0:
        batches.append(lastbatch)

    logger.info("Extracting objects in FOV in parallel using %d processes:",NUMPROCS)
    
    #create multiprocessing pool
    with Pool(NUMPROCS) as p:
        #do this in parallel, returns a list of objects per batch
        results=p.map(partial(get_objects_in_fov,stackedWCS=stackedWCS,ccdBoundaries=ccdBoundaries),batches)
    
    #flatten batched results into one list
    objects = []
    for result in results:
        objects += result

    logger.info("Extracted %d objects in total",len(objects))


    #objects=get_objects_in_fov_from_catalogue(catalogue, stackedWCS, ccdBoundaries)


    output_objects_list = get_qualified_filename(args.workdir,args.output_objects_list)
    logger.info("Writing object list to %s",output_objects_list)
    with open(output_objects_list,"w") as f:
        json.dump(objects,f,indent=4)
    
    combined_catalogue = get_qualified_filename(args.workdir,args.combined_catalogue)
    logger.info("Writing combined catalogue to %s",combined_catalogue)
    catalogue.write(combined_catalogue,overwrite=True)

    t1 = time.time()
    logger.info("Stopped at %s",datetime.datetime.now())
    logger.info("Total runtime = %fs",t1-t0)





