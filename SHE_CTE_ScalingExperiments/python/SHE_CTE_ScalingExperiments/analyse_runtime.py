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
    walltime = []
    num_objects = []
    num_files = []
    start_time = []
    compute_time = []
    for timing in stamp_timings:
        walltime.append(timing["walltime"])
        num_objects.append(timing["num_objects"])
        num_files.append(timing["num_files"])
        start_time.append(timing["tstart"])
        compute_time.append(timing["compute_time"])

    walltime = np.asarray(walltime)
    compute_time = np.asarray(compute_time)
    num_objects = np.asarray(num_objects)
    num_files = np.asarray(num_files)
    start_time = np.asarray(start_time)
    
    #calculate the finish times of each task (start time + walltime)
    finish_times=[]
    for i in range(len(walltime)):
        finish_times.append(start_time[i] + datetime.timedelta(seconds=walltime[i]))
    finish_times = np.asarray(finish_times)
    
    iotime = walltime-compute_time

    iocpuh = np.sum(iotime)/3600
    stampcpuh = np.sum(walltime)/3600.
    
    meaniotime = iotime.mean()
    iotimemin = iotime.min()
    iotimemax = iotime.max()
    iotimestd = iotime.std()

    results={}

    logger.info("Mean I/O time per task = %fs",meaniotime)
    logger.info("Min/Max time = %f/%f s",iotimemin,iotimemax)
    logger.info("Standard Deviation in time = %fs",iotimestd)
    logger.info("CPUh for I/O =  %f",iocpuh)
    logger.info("Total CPUh = %f",stampcpuh)
    
    

    total_walltime = (np.max(finish_times)-np.min(start_time)).total_seconds()

    logger.info("Walltime for all the tasks = %fs",total_walltime)


    results["Mean stamp IO time"]=meaniotime
    results["Min stamp IO time"]=iotimemin
    results["Max stamp IO time"]=iotimemax
    results["Standard Deviation of stamp time"]=iotimestd
    results["Total stamp Walltime"]=total_walltime
    results["Stamp IO CPUh"] = iocpuh
    results["Total stamp CPUh"] = stampcpuh


    #calculate number of concurrent jobs as a function of time
    dt = datetime.timedelta(seconds=1)

    njobs = []
    times = []

    tstart = np.min(start_time)
    tstop = np.max(finish_times)

    t = tstart + dt
    while (t < tstop):
        n=0
        for i in range(len(start_time)):
            if start_time[i] < t and finish_times[i] > t:
                n+=1
        njobs.append(n)
        times.append(t)

        t+= dt

    njobs = np.asarray(njobs)
    
    logger.info("Mean number of running jobs at any given time = %f",njobs.mean())
    logger.info("Max number of running jobs at any given time = %d",njobs.max())

    results["Mean number of running jobs at any given time"] = float(njobs.mean())
    results["Max number of running jobs at any given time"] = int(njobs.max())
    
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

    plt.plot(start_time,iotime,".")
    plt.xlabel("Start time")
    plt.ylabel("I/O time taken per task (s)")
    plt.title("ExtractStamps I/O time")
    filename = os.path.join(plotsdir,"start_time_stamp.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["Start Time Vs I/O time"] = filename
    plt.clf()

    for i in range(len(walltime)):
        tstart = start_time[i]
        tstop = finish_times[i]
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

    results["pipeline_config"] = config

    
    #get statistics on SplitFits
    split_listfile = get_qualified_filename(args.workdir, args.split_timing_listfile)
    with open(split_listfile,"r") as f:
        split_filelist = json.load(f)
    
    files=[]
    for item in split_filelist:
        if type(item) is list:
            files.append(get_qualified_filename(args.workdir,item[0]))
        else:
            files.append(get_qualified_filename(args.workdir,item))

    split_timings = []
    for f in files:
        with open(f,"r") as fd:
            split_timings.append(json.load(fd))
    
    walltimes = []
    exposures=[]
    tstarts = []

    for item in split_timings:
        walltimes.append(item["walltime"])
        exposures.append(int(item["exposure"]))
        tstarts.append(str_to_datetime(item["tstart"]))

    walltimes = np.asarray(walltimes)
    tstarts = np.asarray(tstarts)
    exposures = np.asarray(exposures)

    meantime = walltimes.mean()

    logger.info("Mean time for SplitFits = %fs",meantime)
    results["Mean SplitFits time"] = meantime

    tstops=[]
    for i in range(len(walltimes)):
        tstops.append(tstarts[i]+datetime.timedelta(seconds=walltimes[i]))
    
    tstops = np.asarray(tstops)

    splitcpuh = np.sum(walltimes)/3600.
    logger.info("CPUh for SplitFits = %f",splitcpuh)
    results["Total split CPUh"] = splitcpuh
    results["Total CPUh"] = splitcpuh + stampcpuh






    for i in range(len(walltimes)):
        tstart = tstarts[i]
        tstop = tstops[i]
        plt.plot([tstart,tstop],[exposures[i],exposures[i]],color="blue")
    plt.xlabel("Time")
    plt.ylabel("Exposure")
    plt.title("Timeline of SplitFits tasks")
    filename = os.path.join(plotsdir,"duration_split.png")
    logger.info("Writing plot to %s",filename)
    plt.savefig(filename)
    results["Total duration Vs Exposure"] = filename
    plt.clf()


    plt.bar(exposures,walltimes)
    plt.xlabel("Exposure")
    plt.ylabel("Walltime (s)")
    plt.title("Runtime of SplitFits per exposure")
    filename=os.path.join(plotsdir,"exposure_walltimes.png")
    plt.savefig(filename)
    results["Walltimes for Exposures"] = filename




    




    #write results to file
    output_json = get_qualified_filename(workdir,args.results)
    logger.info("Writing results to %s",output_json)
    with open(output_json,"w") as f:
        json.dump(results,f,indent=4)












