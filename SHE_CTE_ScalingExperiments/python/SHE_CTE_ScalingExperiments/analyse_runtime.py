""" @file analyse_runtime.py

    Created 27 Aug 2021
    

    Generates runtime plots and info for the stamp extraction.
"""

__updated__ = "2021-08-27"

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
import json
import datetime

import numpy as np
import matplotlib.pyplot as plt

from SHE_PPT.logging import getLogger

from .utils import read_config, get_qualified_filename, check_input, create_dummy_output

logger = getLogger(__name__)

def str_to_datetime(time_str):
    """Converts a datetime string (from datetime) back into a datetime object"""
    return datetime.datetime.strptime(time_str,'%Y-%m-%d %H:%M:%S.%f')

def load_timing_json(workdir,timing_file):
    """Loads a timing file and converts its start time string into a datetime"""
    fname = get_qualified_filename(workdir,timing_file)
    with open(fname,"r") as f:
        timing_info = json.load(f)
    timing_info["tstart"] = str_to_datetime(timing_info["tstart"])
    return timing_info




def analyse_runtime_from_args(args):

    workdir = args.workdir
    pipeline_config = args.pipeline_config

    config = read_config(get_qualified_filename(workdir,pipeline_config))

    dry_run = bool(config["dry_run"])

    if dry_run:
        #check all the input files exist, create dummy output files, then exit
        logger.info("Dry run:")

        check_input(args.workdir,args.stamp_timing_listfile,"stamp_timing_listfile")
        check_input(args.workdir,args.split_timing_listfile,"split_timing_listfile")
        
        create_dummy_output(args.workdir,args.results,"results")

        return

    #dictionary to hold the results
    results={}

    #get statistics on SplitFits

    #read listfile
    split_listfile = get_qualified_filename(args.workdir, args.split_timing_listfile)
    with open(split_listfile,"r") as f:
        split_filelist = json.load(f)
    
    #get the individual files to read
    files=[]
    for item in split_filelist:
        if type(item) is list:
            files.append(get_qualified_filename(args.workdir,item[0]))
        else:
            files.append(get_qualified_filename(args.workdir,item))
    
    #read them in
    split_timings = []
    for f in files:
        with open(f,"r") as fd:
            split_timings.append(json.load(fd))
    
    split_stamp_walltimes = []
    exposures=[]
    split_tstarts = []
    
    #extract lists of the stamp_walltimess, exposure numbers and start times
    for item in split_timings:
        split_stamp_walltimes.append(item["walltime"])
        exposures.append(int(item["exposure"]))
        split_tstarts.append(str_to_datetime(item["tstart"]))

    split_stamp_walltimes = np.asarray(split_stamp_walltimes)
    split_tstarts = np.asarray(split_tstarts)
    exposures = np.asarray(exposures)
    
    #calculate the stop time (start time + stamp_walltimes)
    split_tstops=[]
    for i in range(len(split_stamp_walltimes)):
        split_tstops.append(split_tstarts[i]+datetime.timedelta(seconds=split_stamp_walltimes[i]))
    split_tstops = np.asarray(split_tstops)

    split_meantime = split_stamp_walltimes.mean()

    total_splitfits_walltime = (np.max(split_tstops)-np.min(split_tstarts)).total_seconds()



    
    #get statistics for ExtractStamps
    
    #open listfile from the SHE_ExtractStamps runs
    listfile = get_qualified_filename(workdir,args.stamp_timing_listfile)

    with open(listfile,"r") as f:
        filelist = json.load(f)

    files = []
    for f in filelist:
        files.append(f[0])

    stamp_timings = []
    
    #load in the results from the listfiles
    for f in files:
        stamp_timings.append(load_timing_json(workdir,f))
    
    #convert to arrays
    stamp_walltimes = []
    num_objects = []
    num_files = []
    stamp_tstarts = []
    compute_time = []
    for timing in stamp_timings:
        stamp_walltimes.append(timing["walltime"])
        num_objects.append(timing["num_objects"])
        num_files.append(timing["num_files"])
        stamp_tstarts.append(timing["tstart"])
        compute_time.append(timing["compute_time"])

    stamp_walltimes = np.asarray(stamp_walltimes)
    compute_time = np.asarray(compute_time)
    num_objects = np.asarray(num_objects)
    num_files = np.asarray(num_files)
    stamp_tstarts = np.asarray(stamp_tstarts)
    
    #calculate the finish times of each task (start time + stamp_walltimes)
    stamp_tstops=[]
    for i in range(len(stamp_walltimes)):
        stamp_tstops.append(stamp_tstarts[i] + datetime.timedelta(seconds=stamp_walltimes[i]))
    stamp_tstops = np.asarray(stamp_tstops)

    meanstamptime = stamp_walltimes.mean()
    
    iotime = stamp_walltimes-compute_time

    iocpuh = np.sum(iotime)/3600
    stampcpuh = np.sum(stamp_walltimes)/3600.
    
    meaniotime = iotime.mean()
    iotimemin = iotime.min()
    iotimemax = iotime.max()
    iotimestd = iotime.std()


    #calculate number of concurrent jobs as a function of time
    dt = datetime.timedelta(seconds=1)

    njobs = []
    times = []

    tstart = np.min(stamp_tstarts)
    tstop = np.max(stamp_tstops)

    t = tstart + dt
    while (t < tstop):
        n=0
        for i in range(len(stamp_tstarts)):
            if stamp_tstarts[i] < t and stamp_tstops[i] > t:
                n+=1
        njobs.append(n)
        times.append(t)

        t+= dt

    njobs = np.asarray(njobs)



    total_pipeline_walltime = (np.max(stamp_tstops) - np.min(split_tstarts)).total_seconds()
    
    #print/record total pipeline runtime
    logger.info("Total Pipeline Runtime = %fs",total_pipeline_walltime)
    results["Total runtime for whole pipeline"] = total_pipeline_walltime

    #print/record SplitFits info
    logger.info("Total SplitFits Walltime = %fs",total_splitfits_walltime)
    results["Total SplitFits Walltime"] = total_splitfits_walltime
    
    logger.info("Mean walltime for SplitFits tasks = %fs",split_meantime)
    results["Mean SplitFits task walltime"] = split_meantime

    splitcpuh = np.sum(split_stamp_walltimes)/3600.
    logger.info("Total CPUh for SplitFits tasks = %f",splitcpuh)
    results["Total SplitFits CPUh"] = splitcpuh


    
    #print/record ExtractStamps info

    total_stamp_walltimes = (np.max(stamp_tstops)-np.min(stamp_tstarts)).total_seconds()

    results["Total ExtractStamps walltime"]=total_stamp_walltimes
    logger.info("Total ExtractStamps walltime = %fs",total_stamp_walltimes)

    logger.info("Mean walltime for ExtractStamps tasks = %fs",meanstamptime)
    logger.info("Mean I/O time per ExtractStamps task = %fs",meaniotime)
    logger.info("Min/Max I/O time per ExtractStamps task = %f/%f s",iotimemin,iotimemax)
    logger.info("Standard Deviation in ExtractStamps I/O time = %fs",iotimestd)
    logger.info("CPUh for ExtractStamps I/O =  %f",iocpuh)
    logger.info("Total CPUh for ExtractStamps = %f",stampcpuh)
    
    results["Mean walltime for ExtractStamps tasks"]=meanstamptime
    results["Mean ExtractStamps task I/O time"]=meaniotime
    results["Min ExtractStamps task I/O time"]=iotimemin
    results["Max ExtractStamps task I/O time"]=iotimemax
    results["Standard Deviation of ExtractStamps task I/O time"]=iotimestd
    
    results["ExtractStamps I/O CPUh"] = iocpuh
    results["Total ExtractStamps CPUh"] = stampcpuh

    logger.info("Mean number of running ExtractStamps tasks at any given time = %f",njobs.mean())
    logger.info("Max number of running ExtractStamps tasks at any given time = %d",njobs.max())

    results["Mean number of running ExtractStamps tasks at any given time"] = float(njobs.mean())
    results["Max number of running ExtractStamps tasks at any given time"] = int(njobs.max())
    
    plotsdir=os.path.join(workdir,"plots")
    try:
        os.mkdir(plotsdir)
    except(FileExistsError):
        logger.info("Plotting directory %s already exists",plotsdir)
        pass

    #make various plots
    
    plt.plot(num_objects, iotime,".")
    plt.xlabel("Number of Objects per batch")
    plt.ylabel("I/O time taken per task (s)")
    plt.title("ExtractStamps I/O time")
    filename = os.path.join(plotsdir,"iotime_objects.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["I/O time Vs Num Objects"] = filename
    plt.clf()

    plt.plot(num_files, iotime,".")
    plt.xlabel("Number of files accessed per batch")
    plt.ylabel("I/O time taken per task (s)")
    plt.title("ExtractStamps I/O time")
    filename = os.path.join(plotsdir,"iotime_files.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["I/O time Vs Num Files"] = filename
    plt.clf()

    plt.plot(stamp_tstarts,iotime,".")
    plt.xlabel("Start time")
    plt.ylabel("I/O time taken per task (s)")
    plt.title("ExtractStamps I/O time")
    filename = os.path.join(plotsdir,"start_time_stamp.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["Start Time Vs I/O time"] = filename
    plt.clf()

    for i in range(len(stamp_walltimes)):
        tstart = stamp_tstarts[i]
        tstop = stamp_tstops[i]
        plt.plot([tstart,tstop],[i,i],color="blue")
    plt.xlabel("Time")
    plt.ylabel("Task number")
    plt.title("Timeline of ExtractStamps tasks")
    filename = os.path.join(plotsdir,"duration_batchno.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["Total duration Vs Batch Number"] = filename
    plt.clf()


    plt.plot(times,njobs)
    plt.xlabel("Time")
    plt.ylabel("Number of running tasks")
    plt.title("Number of Running ExtractStamps tasks with time")
    filename = os.path.join(plotsdir,"num_running_tasks.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["Number of running tasks with time"] = filename
    plt.clf()


    for i in range(len(split_stamp_walltimes)):
        tstart = split_tstarts[i]
        tstop = split_tstops[i]
        plt.plot([tstart,tstop],[exposures[i],exposures[i]],color="blue")
    plt.xlabel("Time")
    plt.ylabel("Exposure")
    plt.title("Timeline of SplitFits tasks")
    filename = os.path.join(plotsdir,"duration_split.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["Total duration Vs Exposure"] = filename
    plt.clf()


    plt.bar(exposures,split_stamp_walltimes)
    plt.xlabel("Exposure")
    plt.ylabel("stamp_walltimes (s)")
    plt.title("Runtime of SplitFits per exposure")
    filename=os.path.join(plotsdir,"split_stamp_exposure_walltimes.png")
    plt.savefig(filename)
    results["split_stamp_walltimess for Exposures"] = filename

    #add pipeline_config to the results dict
    results["pipeline_config"] = config




    #write results to file
    output_json = get_qualified_filename(workdir,args.results)
    logger.info("Writing results to %s",output_json)
    with open(output_json,"w") as f:
        json.dump(results,f,indent=4)
    
    #also write it to the plots directory
    with open(os.path.join(plotsdir,"results.json"),"w") as f:
        json.dump(results,f,indent=4)
    












