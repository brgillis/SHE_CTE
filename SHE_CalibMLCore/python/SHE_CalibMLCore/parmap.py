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
File: python/SHE_CalibMLCore/parmap.py

Created on: 09/05/17

http://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class

"""


from multiprocessing import Process, Pipe, Queue
import multiprocessing


def spawn(f):
    """
    Helper function for `parmap` which prepares the worker `f` to be spawned.
    """
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def parmap(f, X, ncpu):
    """
    This is an alternative to multiprocessing.Pool to avoid the limitations of the package (pickling stuff...)

    .. note:: see http://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class

    .. note:: It is very possible that multiprocessing.Pool is fixed in python3 
    """
    pipe = [Pipe() for x in X]
    processes = [Process(target=spawn(f), args=(c, x)) for x, (p, c) in zip(X, pipe)]
    numProcesses = len(processes)
    processNum = 0
    outputList = []
    while processNum < numProcesses:
        endProcessNum = min(processNum + ncpu, numProcesses)
        for proc in processes[processNum:endProcessNum]:
            proc.start()
        # It is crucial to call recv() before join() to avoid deadlocks !
        for proc, c in pipe[processNum:endProcessNum]:
            outputList.append(proc.recv())
        for proc in processes[processNum:endProcessNum]:
            proc.join()

        processNum = endProcessNum
    return outputList


if __name__ == '__main__':
    print(parmap(lambda i: i * 2, [1, 2, 3, 4, 6, 7, 8], 3))
