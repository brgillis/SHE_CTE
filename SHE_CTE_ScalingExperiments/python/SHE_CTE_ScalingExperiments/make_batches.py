""" @file make_batches.py

    Created 26 Aug 2021
    

    Batches the objects for extracting postage stamps.
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

from .utils import read_config, get_qualified_filename, check_input, create_dummy_output_listfile
from SHE_PPT.logging import getLogger

logger = getLogger(__name__)




def batch_spatially(objects, batchsize):
    """Crude spatial batching. Split the image into many rectangular regions, and sort objects into these"""
    logger.info("Batching spatially")

    t = len(objects)

    #number of batches needed for batchsize objects per batch
    n_2dbatches = t/batchsize
    #number of batches in each direction
    n_1dbatches = int(np.sqrt(n_2dbatches))

    logger.info("Total number of batches = %d", n_2dbatches)
    logger.info("Number of batches in each direction = %d", n_1dbatches)
    
    #determine max/min extent of the objects in x and y
    xmin = np.inf
    xmax = 0
    ymin = np.inf
    ymax = 0
    for obj in objects:
        if obj["x"] < xmin:
            xmin = np.floor(obj["x"])
        if obj["x"] > xmax:
            xmax = np.ceil(obj["x"])
        if obj["y"] < ymin:
            ymin = np.floor(obj["y"])
        if obj["y"] > ymax:
            ymax = np.ceil(obj["y"])

    #determine the spatial extent of each batch in x and y
    dx = (xmax-xmin)/(n_1dbatches-1)
    dy = (ymax-ymin)/(n_1dbatches-1)

    logger.info("Spatial extent of each batch in pixels: (x, y) = (%d,%d)",dx,dy)
    
    #array for counting the number of objects in each batch
    count = np.zeros((n_1dbatches, n_1dbatches),dtype=np.int32)
    
    #make a 2d array of lists for storing the objects in each batch
    batches2d = [[ [] for _ in range(n_1dbatches)] for _ in range(n_1dbatches)]

    #sort the objects into 2d batches
    for obj in objects:
        #get the coordinates of the object
        x = obj["x"]
        y = obj["y"]
        
        #get the 2-d batch number
        i = int((x-xmin)//dx)
        j = int((y-ymin)//dy)

        count[j,i] +=1

        batches2d[j][i].append(obj)

    mean_num = np.mean(count)
    median_num = np.median(count)

    logger.info("Mean/Median number of objects per batch = %f/%f",mean_num, median_num)
    logger.info("Max/Min number of objects per batch = %f,%f",count.max(),count.min())
    logger.info("Standard deviation of batch count = %f",np.std(count))

    #flatten bins into a 1d list
    batches=[]

    for row in batches2d:
        batches += row

    return batches


def batch_by_order_in_list(objects,batchsize):
    """Batch objects by the order they appear in the objects list"""

    logger.info("Batching objects by the order they appear in the list")
    numbatches = len(objects)//batchsize
    
    batches=[]
    for i in range(numbatches):
        batch = objects[i*batchsize:(i+1)*batchsize]
        batches.append(batch)

    lastbatch = objects[numbatches*batchsize:]
    if len(lastbatch) > 0:
        batches.append(lastbatch)

    logger.info("Split objects into %d batches",len(batches))

    return batches

    
    

       


def make_batches_from_args(args):
    t0 = time.time()
    logger.info("Started at %s",datetime.datetime.now())

    #get required config options
    config = read_config(get_qualified_filename(args.workdir,args.pipeline_config))
    maxbatches = config["maxbatches"]
    batchsize = config["batchsize"]
    spatial_batching = bool(config["spatial_batching"])
    dry_run = bool(config["dry_run"])

    if dry_run:
        #check all the input files exist, create dummy output files, then exit
        logger.info("Dry run:")

        check_input(args.workdir,args.objects_list,"objects_list")
        if maxbatches <=0:
            maxbatches=1

        create_dummy_output_listfile(args.workdir,args.batch_listfile,"batch_listfile",maxbatches)

        return
    
    #read in the objects
    objects_file = get_qualified_filename(args.workdir,args.objects_list)
    with open(objects_file,"r") as f:
        objects = json.load(f)

    #first batch all the objects
    if spatial_batching:
        batches = batch_spatially(objects,batchsize)
    else:
        batches = batch_by_order_in_list(objects,batchsize)
    
    #set the random seed for reproducable results
    random_state = np.random.RandomState(seed=1)

    # shuffle the batches array so that if we want a small number of batches we aren't 
    # going to be selecting the first few batches which may be atypical as they are on
    # the edge of the image
    random_state.shuffle(batches)
    
    # if maxbatches < 0 we assume we want all the batches
    if maxbatches < 0:
        maxbatches = len(batches)
    
    #write all the batches to separate files
    try:
        os.mkdir(os.path.join(args.workdir,"batches"))
    except FileExistsError:
        logger.warn("Batch subdirectory already exists")
    
    batchedfiles = []

    logger.info("Writing %d batch files...",maxbatches)
    for i in range(maxbatches):
        filename = "batch_%06d.json"%i
        relative_filename = os.path.join("batches",filename)
        batchedfiles.append(relative_filename)
        #write batch
        with open(get_qualified_filename(args.workdir,relative_filename),"w") as f:
            json.dump(batches[i],f,indent=4)

    batch_listfile = get_qualified_filename(args.workdir,args.batch_listfile)
    logger.info("Writing listfile of the batches to %s",batch_listfile)
    with open(batch_listfile,"w") as f:
        json.dump(batchedfiles,f,indent=4)

    t1 = time.time()
    logger.info("Stopped at %s",datetime.datetime.now())
    logger.info("Total time taken = %fs",t1-t0)
    

    
        
