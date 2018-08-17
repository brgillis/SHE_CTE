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
File: python/SHE_CalibMLCore/train.py

Created on: 09/05/17
Author: Malte Tewes

Training a network happens here
"""


import numpy as np
from datetime import datetime
import os
import copy

import logging
logger = logging.getLogger(__name__)

from . import layer
from . import utils
from . import err
from . import act
#from . import plot
from . import regul
from . import opt


class Training:
    """
    Holds together everthing related to the process of training a Net (or a WNet): the training data and the network.
    """

    def __init__(self, net, dat, errfctname="msb", regulweight=None, regulfctname=None, itersavepath=None, saveeachit=False, autoplotdirpath=".", autoplot=False, logpath=None, trackbiases=False, verbose=False, name=None):
        """

        Sets up
        - housekeeping lists
        - error function

        :param dat: can be None, if you want a dummy container for your Net.

        :param trackbiases: If True, will track and plot the evolution of biases as function of input feature values.
            Could get massive and slow down the training. Also the plot is rather slow.
        :type trackbiases: bool

        :param logpath: If not None, I will log to this file, and not propagate my logger at all.
            Unfortunately this does not work well when using multiple threads/processes,
            as logs from every thread get written into the logfiles.

        :param saveeachit: If True, writes the full status and history to disk at each iteration.
            If False, only some snapshots are written, but they still contain the full history!
        :type saveeachtit: bool


        """

        # We start by dealing with the logger.
        if logpath is not None:
            fh = logging.FileHandler(logpath, delay=True)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter("PID %(process)d: %(levelname)s: %(name)s(%(funcName)s): %(message)s"))
            logger.addHandler(fh)
            logger.propagate = False
            # Even if the file is only openend at the first message (delay=True in FileHandler above),
            # we do have to make sure the directories exist, so that logging can be done at any time.
            logpathdir = os.path.split(logpath)[0]
            if not os.path.exists(logpathdir):
                os.makedirs(logpathdir)

        self.name = name

        self.set_net(net)
        self.set_dat(dat)

        # Let's check compatibility between those two!
        if dat != None:
            assert net.ni == self.dat.getni()
            if errfctname in ["mse", "msb", "msrb", "msre", "msbw", "msbwsignorm"]:
                assert net.no == self.dat.getno()
            elif errfctname in ["msbwnet"]:
                assert net.no == 2 * self.dat.getno()
            else:
                logger.warning("Unknown error function, will blindly go ahead...")

        # Setting up the cost function
        self.errfctname = errfctname
        self.errfct = eval("err.{0}".format(self.errfctname))
        if regulfctname is not None:
            if regulweight is None:
                raise ValueError("Regularisation is set, but no weight")
            else:
                if regulweight == 0:
                    logger.warning("Regularisation is set, but the weight is zero! Let's boldy go ahead anyhow.")
                elif regulweight < 0:
                    logger.warning("Regularisation is set, but the weight is negative! Let's boldy go ahead anyhow.")
                self.regulfctname = regulfctname
                self.regulfct = eval("regul.{0}".format(self.regulfctname))
                self.regullam = regulweight
        else:
            self.regullam = None

        # We initialize some counters for the optimization:
        self.optit = 0  # The iteration counter
        self.optcall = 0  # The cost function call counter
        self.optitcall = 0  # Idem, but gets reset at each new iteration
        self.opterr = np.inf  # The current cost function value on the training set

        # And some lists describing the optimization:
        self.opterrs = []  # The cost function value on the training set at each (!) call

        self.optitparams = []  # A copy of the network parameters at each iteration
        self.optiterrs_train = []  # The cost function value on the training set at each iteration
        self.optiterrs_val = []  # The cost function value on the validation set at each iteration

        self.optitcalls = []  # The cost function call counter at each iteration
        self.optittimes = []  # Time taken for iteration, in seconds
        self.optbatchchangeits = []  # Iteration counter (in fact indices) when the batche gets changed

        # Should we save all the biases at each snapshot (e.g. minibatch-change) (and make the related plots) ?
        self.trackbiases = trackbiases
        self.biassnaps_train = []
        self.biassnaps_val = []
        self.biassnaps_it = []  # Iteration counter of snapshot

        self.verbose = verbose
        self.itersavepath = itersavepath
        self.saveeachit = saveeachit

        self.autoplotdirpath = autoplotdirpath
        self.autoplot = autoplot

        if dat is not None:  # Otherwise it's just to create a dummy container, no need to report
            logger.info("Done with setup of {self}".format(self=self))

        # And let's test this out before we start, so that it fails fast in case of a problem.
        # But we only do so if the file does not yet exist, to avoid overwriting an existing training.
        # While this idea was OK for debugging, it's just weird to save anything to disk at this stage, so let's remove it.
        # This should be checked when starting the training.
        # if self.itersavepath is not None:
        #    if not os.path.exists(self.itersavepath):
        #        self.save(self.itersavepath)

    def set_dat(self, dat):
        """
        Allows to add or replace training data (e.g. when reading a self.save()...)
        """
        self.dat = dat
        self.set_datstr()

    def set_datstr(self):
        """
        We do this so that we can still display this string on plots if dat has been removed.
        """
        self.datstr = str(self.dat)

    def set_net(self, net):
        """
        Replaces the network object
        """
        self.net = net
        self.params = self.net.get_params_ref()  # Fast connection to the network parameters
        self.paramslice = slice(None)  # By default, all params are free to be optimized

    def set_paramslice(self, **kwargs):
        """
        The paramslice allows to specify which params you want to be optimized.
        This is relevant for instance when training a MultNet or a WNet.
        We use a slice of this. Indexing with a boolean array ("mask") would seem nicer, but fancy indexing does not preserve
        the references. Hence using a slice is a good compromise for speed.

        The Net object knows how to get such a slice, given the kwargs and the pecularities of the Net.
        """

        self.paramslice = self.net.get_paramslice(**kwargs)
        logger.info("Paramslice is set with kwargs '{}' : {}/{} params are free to be optimized.".format(
            kwargs, len(self.params[self.paramslice]), self.net.nparams())
        )

    def __str__(self):
        """
        A short spaceless automatic description
        """
        if self.dat is not None:  # If available, we use the live dat.
            datstr = str(self.dat)
        else:  # Otherwise, we use datstr as backup solution.
            # To ensure that it works also if datstr is missing (backwards compatibility).
            datstr = getattr(self, "datstr", "Ouch!")

        autotxt = "{ef}({self.net}, {datstr})".format(ef=self.get_costfctname(), self=self, datstr=datstr)
        return autotxt

    def get_costfctname(self):
        """
        Returns a string with a description of the cost function
        """
        txterrfct = "{self.errfctname}".format(self=self)
        if hasattr(self, 'regullam') and self.regullam is not None:
            txterrfct += "+{self.regullam:1.0e}x{self.regulfctname}".format(self=self)
        return txterrfct

    def title(self):
        """
        Returns the name and string, typically nicer for plots.
        """

        if self.name is not None:
            return "Training '{name}': {auto} ({it} it, {tmin:.1f} min)".format(name=self.name, auto=str(self), it=self.optit, tmin=np.sum(self.optittimes) / 60.0)
        else:
            return "{auto} ({it} it, {tmin:.1f} min)".format(auto=str(self), it=self.optit, tmin=np.sum(self.optittimes) / 60.0)

    def takeover(self, othertrain):
        """
        This copies the net and all the progress counters and logs from another train object into the current one.
        Useful if you want to carry on a training with a new train object, typically with different settings and
        maybe on different data.
        """

        logger.info("Setting up training '{}' to take over the work from '{}'...".format(self.name, othertrain.name))

        if not self.net.nparams() == othertrain.net.nparams():
            raise RuntimeError("Other network is not compatible, this is fishy!")

        self.set_net(othertrain.net)

        # Copying the counter positions:
        self.optit = othertrain.optit
        self.optcall = othertrain.optcall
        self.optitcall = 0  # Gets reset at each new iteration
        self.opterr = othertrain.opterr

        # And some lists describing the optimization:
        self.opterrs = othertrain.opterrs[:]

        self.optitparams = othertrain.optitparams[:]
        self.optiterrs_train = othertrain.optiterrs_train[:]
        self.optiterrs_val = othertrain.optiterrs_val[:]

        self.optitcalls = othertrain.optitcalls[:]
        self.optittimes = othertrain.optittimes[:]
        try:
            self.optbatchchangeits = othertrain.optbatchchangeits[:]
        except AttributeError:
            self.optbatchchangeits = []

        self.biassnaps_train = othertrain.biassnaps_train
        self.biassnaps_val = othertrain.biassnaps_val
        self.biassnaps_it = othertrain.biassnaps_it

        logger.info("Done with the takeover")

    def save(self, filepath, keepdata=False, resetcache=True):
        """
        Saves the training progress into a pkl file
        As the training data is so massive, by default we do not save it!
        Note that this might be done at each iteration!
        """
        if resetcache is True:
            self.net.resetcache()
        if keepdata is True:
            logger.info("Writing training to disk and keeping the data...")
            utils.writepickle(self, filepath)
        else:
            tmptraindata = self.dat
            self.set_datstr()
            self.dat = None
            utils.writepickle(self, filepath)
            self.dat = tmptraindata

    def plotpath(self, plotname, suffix="_optitXXXXX", dirpath=None):
        """
        Prepares a path to be used by checkplots.
        """

        if dirpath is None:
            dirpath = self.autoplotdirpath

        if suffix == "_optitXXXXX":
            suffix = "_optit{0:05d}".format(self.optit)

        return os.path.join(dirpath, plotname + suffix + ".png")

    def makeplots(self, **kwargs):
        """
        Saves a bunch of default checkplots into the specified directory.
        Can typically be called at the end of training, or after iterations.

        kwargs are passed to plotpath.
        """

        logger.info("Making and writing plots...")
        plot.sumevo(self, filepath=self.plotpath("sumevo", **kwargs))
        plot.outdistribs(self, filepath=self.plotpath("outdistribs", **kwargs))
        plot.errorinputs(self, filepath=self.plotpath("errorinputs", **kwargs))
        self.plot_netviz(**kwargs)
        if self.trackbiases:
            plot.biasevo(self, filepath=self.plotpath("biasevo", **kwargs))

        logger.info("Done with plots")

    def plot_netviz(self, filepath='None', **kwargs):
        """
        Slightly different from the other plots as it takes a Net, not a Training.
        """

        #write_mode = False
        # for l in self.net.layers:
        #    if not l.mode == "sum":
        #        write_mode = True
        #
        # if write_mode:
        #    title = r"$\to$".join([r"${m}^{{{n}}}\mathrm{{{actfct}}}$".format(n=l.nn, m=plot.get_symbol(l.mode, False), actfct=l.actfct.__name__) for l in self.net.layers])
        # else:
        #    title = r"$\to$".join([r"$\mathrm{{{n}/{actfct}}}$".format(n=l.nn, actfct=l.actfct.__name__) for l in self.net.layers])
        #title += r"$;\ \mathrm{{{}}}$".format(self.get_costfctname())
        #title += r"$;\ \mathrm{{it:{:05d}}}$".format(self.optit)
        #title += r"$;\ \mathrm{{{tmin:.1f}min}}$".format(tmin=np.sum(self.optittimes)/60.0)

        title = str(self) + ", {} it, {:.1f} min".format(self.optit, np.sum(self.optittimes) / 60.0)

        if filepath is 'None':
            filepath = self.plotpath("netviz")
        plot.netviz(self.net, filepath=filepath, title=title, **kwargs)

    def start(self):
        """
        Called a the beginning of a training 
        """
        logger.info("Starting an optimization cycle...")
        self.testcost()
        self.iterationstarttime = datetime.now()
        self.optitcall = 0

        # If this is the very first start, we prepare plots to viz the initial contidions:
        if self.optit == 0 and self.autoplot:
            self.plot_netviz()

    def end(self):
        """
        Called at the end of a training (each minibatch) depening on the algo.
        This is also a moment to save what we have to disk.
        """
        self.optitcall = 0
        logger.info("Cumulated training time: {0:.2f} s".format(np.sum(self.optittimes)))
        logger.info("Optimization cycle finished.")
        if self.trackbiases:
            self.savebiasdetails()
        if self.itersavepath != None:
            if not self.saveeachit:  # Indeed, otherwise no need to save it at this stage!
                self.save(self.itersavepath)
        if self.autoplot:
            self.makeplots()

    def callback(self, *args):
        """
        Function called by the optimizer after each "iteration".
        Print out some info about the training progress,
        saves status of the counters,
        and optionally writes the network itself to disk.
        """
        # print args
        # exit()

        self.optit += 1
        now = datetime.now()
        secondstaken = (now - self.iterationstarttime).total_seconds()
        callstaken = self.optitcall

        # Not sure if it is needed to update the params (of if the optimizer
        # already did it), but it cannot harm and is fast:
        self.params[self.paramslice] = args[0]  # Updates the network parameters

        self.optittimes.append(secondstaken)
        self.optiterrs_train.append(self.opterr)
        self.optitcalls.append(self.optcall)
        self.optitparams.append(copy.deepcopy(self.params))  # We add a copy of the current params

        # Now we evaluate the cost on the validation set:
        valerr = self.valcost()
        self.optiterrs_val.append(valerr)

        valerrratio = valerr / self.opterr

        mscallcase = 1000.0 * float(secondstaken) / (float(callstaken) *
                                                     self.dat.getntrain())  # Time per call and training case

        logger.info("Iter. {self.optit:4d}, {ef} train = {self.opterr:.6e}, val = {valerr:.6e} ({valerrratio:4.1f}), {time:.4f} s for {calls} calls ({mscallcase:.4f} ms/cc)".format(
            self=self, ef=self.get_costfctname(), time=secondstaken, valerr=valerr, valerrratio=valerrratio, calls=callstaken, mscallcase=mscallcase))

        if self.itersavepath != None:
            if self.saveeachit:
                self.save(self.itersavepath)

        # We reset the iteration counters:
        self.iterationstarttime = datetime.now()
        self.optitcall = 0

        # And now we take care of getting a new batch
        # self.randombatch()

    def cost(self, p):
        """
        The "as-fast-as-possible" function to compute the training error based on parameters p.
        This gets called repeatedly by the optimizers.
        """

        self.params[self.paramslice] = p  # Updates the network parameters

        # Compute the outputs
        outputs = self.net.run(self.dat.traininputs)  # This is not a masked array!

        # And now evaluate the error (cost) function.
        if self.dat.trainoutputsmask is not None:
            outputs = np.ma.array(outputs, mask=self.dat.trainoutputsmask)
        err = self.errfct(outputs, self.dat.traintargets, auxinputs=self.dat.trainauxinputs)

        if self.regullam is not None:
            # TODO: Is there a faster way to do this?
            weights = np.concatenate(([(l.weights).flatten() for l in self.net.layers]))
            err += self.regullam * self.regulfct(weights)

        self.opterr = err
        self.optcall += 1
        self.optitcall += 1
        self.opterrs.append(err)

        if self.verbose:
            logger.debug(
                "Iteration {self.optit:4d}, call number {self.optcall:8d}: cost = {self.opterr:.8e}".format(self=self))
            logger.debug("\n" + self.net.report())

        return err

    def currentcost(self):
        return self.cost(p=self.params[self.paramslice])

    def testcost(self):
        """
        Calls the cost function and logs some info.
        """

        logger.info("Testing cost function calls...")
        starttime = datetime.now()
        err = self.currentcost()
        endtime = datetime.now()
        took = (endtime - starttime).total_seconds()
        logger.info("On the training set:   {took:.4f} seconds, {self.errfctname} = {self.opterr:.8e}".format(
            self=self, took=took))
        starttime = datetime.now()
        err = self.valcost()
        endtime = datetime.now()
        took = (endtime - starttime).total_seconds()
        logger.info("On the validation set: {took:.4f} seconds, {self.errfctname} = {err:.8e}".format(
            self=self, took=took, err=err))

    def valcost(self):
        """
        Evaluates the cost function on the validation set.
        """
        outputs = self.net.run(self.dat.valinputs)  # This is not a masked array!
        if self.dat.valoutputsmask is not None:
            outputs = np.ma.array(outputs, mask=self.dat.valoutputsmask)

        err = self.errfct(outputs, self.dat.valtargets, auxinputs=self.dat.valauxinputs)

        return err

    def savebiasdetails(self):
        """
        Appends bias details to some internal lists.
        For plotting purposes only.
        A bit ugly, all this duplication :-/
        Note that we use the FULL training data, not only the current batch.
        """
        logger.info("Saving bias details...")
        fulltrainoutputs = self.net.run(self.dat.fulltraininputs)  # This is not a masked array!
        if self.dat.fulltrainoutputsmask is not None:
            fulltrainoutputs = np.ma.array(fulltrainoutputs, mask=self.dat.fulltrainoutputsmask)
        fulltrainbiases = np.mean(fulltrainoutputs, axis=0) - self.dat.fulltraintargets

        valoutputs = self.net.run(self.dat.valinputs)  # This is not a masked array!
        if self.dat.valoutputsmask is not None:
            valoutputs = np.ma.array(valoutputs, mask=self.dat.valoutputsmask)
        valbiases = np.mean(valoutputs, axis=0) - self.dat.valtargets

        # Those 2 arrays have shape (neuron, cases)
        self.biassnaps_train.append(fulltrainbiases)
        self.biassnaps_val.append(valbiases)
        self.biassnaps_it.append(self.optit)  # We record the iteration counter

    def opt(self, algo="brute", mbsize=None, mbfrac=0.1, mbloops=10, **kwargs):
        """
        General optimization of the network parameters by an algorithm, with support for minibatches.
        Just set mbfrac and mbloops to 1 if you don't want minibatches.
        Using "algo" and not "method" so that "method" can be a kwarg passed to scipy.optimize.fmin...
        """

        optfct = eval("opt.{0}".format(algo))  # Turn the string into an acutal function

        for loopi in range(mbloops):
            if mbloops > 1:
                logger.info("Starting minibatch loop {loopi} of {mbloops}...".format(loopi=loopi + 1, mbloops=mbloops))
            self.dat.random_minibatch(mbsize=mbsize, mbfrac=mbfrac)
            self.optbatchchangeits.append(self.optit)  # We record this minibatch change

            self.start()
            optfct(self, **kwargs)  # Call the optfct, with the training object as first argument

            logger.info("Done with optimization")
            self.end()
