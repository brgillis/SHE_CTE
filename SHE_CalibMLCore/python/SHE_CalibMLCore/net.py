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
File: python/SHE_CalibMLCore/net.py

Created on: 09/05/17
Author: Malte Tewes


Net represent a "simple" network.

"""


import numpy as np

import logging
logger = logging.getLogger(__name__)

from . import layer
from . import utils
from . import act
from . import data


class Net():
    """
    Object representing a network made out of one or several hidden layers.
    """

    def __init__(self, ni, nhs, no=1, onlyid=False, actfctname="tanh", oactfctname="iden", multactfctname="iden", name=None, inames=None, onames=None):
        """
        :param ni: Number of input features
        :param nhs: Numbers of neurons in hidden layers. Times -1 if you want product units in that layer.
        :type nhs: tuple
        :param no: Number of ouput neurons
        :param onlyid: Set this to true if you want identity activation functions on all layers
            (useful for debugging). Note that this overrides both actfctname and oactfctname!

        :param actfctname: name of activation function for hidden layers
        :param oactfctname: idem for output layer

        :param name: if None, will be set automatically
        :type name: string

        :param inames: a list of names (strings) for the input nodes, to be used, e.g., in checkplots.
            These names have a purely decorative function, and are optional.
        :param onames: idem, for the ouptut nodes.


        """

        self.ni = ni
        self.nhs = nhs
        self.no = no
        self.name = name

        # We take care of the inames and onames:
        if inames is None:
            self.inames = ["i_" + str(i) for i in range(self.ni)]
        else:
            self.inames = inames
            if len(self.inames) != self.ni:
                raise RuntimeError("Your number of inames is wrong")
        if onames is None:
            self.onames = ["o_" + str(i) for i in range(self.no)]
        else:
            self.onames = onames
            if len(self.onames) != self.no:
                raise RuntimeError("Your number of onames is wrong")

        # Note that we do not save this. Layers might evolve dynamically in future!
        iniarch = np.array([self.ni] + self.nhs + [self.no])

        actfct = eval("act.{0}".format(actfctname))  # We turn the string actfct option into an actual function
        oactfct = eval("act.{0}".format(oactfctname))  # idem
        multactfct = eval("act.{0}".format(multactfctname))  # idem

        self.layers = []  # We build a list containing only the hidden layers and the output layer
        for (i, nh) in enumerate(self.nhs):
            if nh > 0:
                self.layers.append(layer.Layer(ni=iniarch[i], nn=nh, actfct=actfct, name="h" + str(i), mode="sum"))
            elif nh < 0:  # a "mult"-layer
                self.layers.append(layer.Layer(ni=iniarch[i], nn=-nh,
                                               actfct=multactfct, name="h" + str(i), mode="mult"))
            else:
                raise ValueError("Cannot have 0 hidden nodes")
        # Adding the output layer:
        self.layers.append(layer.Layer(ni=self.nhs[-1], nn=no, actfct=oactfct, name="o"))

        if onlyid:  # Then all layers get the Id activation function:
            for l in self.layers:
                l.actfct = act.iden

        logger.info("Built " + str(self))

    def __str__(self):
        """
        A short string describing the network
        """
        archtxt = str(self.ni) + "|" + "|".join(["{n}/{modecode}{actfct}".format(n=l.nn,
                                                                                 modecode=l.modecode(), actfct=l.actfct.__name__) for l in self.layers])
        autotxt = "[{archtxt}={nparams}]".format(archtxt=archtxt, nparams=self.nparams())

        if self.name is None:
            return autotxt
        else:
            return "'{name}' {autotxt}".format(name=self.name, autotxt=autotxt)

    def report(self):
        """
        Returns a text about the network parameters, useful for debugging.
        """
        txt = ["=" * 120, str(self)]
        for l in self.layers:
            txt.append(l.report())
        txt.append("=" * 120)
        return "\n".join(txt)

    def resetcache(self):
        """
        Purges the cache of every layer (e.g., to save memory before writing to file)
        """
        for l in self.layers:
            l.resetcache()

    def save(self, filepath, resetcache=True):
        """
        Saves self into a pkl file
        """
        if resetcache:
            self.resetcache()
        utils.writepickle(self, filepath)

    def nparams(self):
        """
        Returns the number of parameters of the network
        """
        return sum([l.nparams() for l in self.layers])

    def get_params_ref(self):
        """
        Get a single 1D numpy array containing references to all network weights and biases.
        Note that each time you call this, you loose the "connection" to the ref from any previous calls.

        Warning: if you modify this code, pay attention to update the paramslice-related stuff as well!
        """

        ref = np.empty(self.nparams())
        ind = 0

#        #First layer first, weights and biases.
#            for l in self.layers:
#                ref[ind:ind+(l.nn*l.ni)] = l.weights.flatten() # makes a copy
#                ref[ind+(l.nn*l.ni):ind+l.nparams()] = l.biases.flatten() # makes a copy
#                l.weights = ref[ind:ind+(l.nn*l.ni)].reshape(l.nn, l.ni) # a view
#                l.biases = ref[ind+(l.nn*l.ni):ind+l.nparams()] # a view
#                ind += l.nparams()

        # Starting at the end, biases before weights

        for l in self.layers[::-1]:

            ref[ind:ind + l.nn] = l.biases.flatten()  # makes a copy
            ref[ind + l.nn:ind + l.nparams()] = l.weights.flatten()  # makes a copy
            l.biases = ref[ind:ind + l.nn]  # a view
            l.weights = ref[ind + l.nn:ind + l.nparams()].reshape(l.nn, l.ni)  # a view
            ind += l.nparams()

        # Note that such tricks do not work, as indexing by indices creates copies:
        #indices = np.arange(self.nparams())
        # np.random.shuffle(indices)
        # return ref[indices]

        assert ind == self.nparams()
        return ref

    def get_paramlabels(self):
        """
        Returns a list with labels describing the params. This is for humans and plots.
        Note that plots might expect these labels to have particular formats.
        """

        paramlabels = []
        ind = 0

        # if schema == 2:
        for l in self.layers[::-1]:

            paramlabels.extend(l.nn * ["layer-{l.mode}-{l.name}_bias".format(l=l)])
            paramlabels.extend(l.nn * l.ni * ["layer-{l.mode}-{l.name}_weight".format(l=l)])

        assert len(paramlabels) == self.nparams()

        return paramlabels

    def addnoise(self, **kwargs):
        """
        Adds random noise to all parameters.
        """

        logger.info("Adding noise to network parameters ({})...".format(str(kwargs)))

        for l in self.layers:
            l.addnoise(**kwargs)

    def setidentity(self, **kwargs):
        """
        Adjusts the network parameters so to approximatively get an identity relation
        between the ith output and the ith input (for each i in the outputs).

        This should be a good starting position for "calibration" tasks. Example: first
        input feature is observed galaxy ellipticity e1, and first output is shear g1.
        """

        for l in self.layers:
            if l.nn < self.no or self.ni < self.no:
                raise RuntimeError("Network is too small for setting identity!")

        for l in self.layers:
            l.setidentity(**kwargs)

        logger.info("Done with setting weights to identity.")

    def run(self, inputs):
        """
        Propagates input through the network "as fast as possible".
        This works for 1D, 2D, and 3D inputs, see layer.run().
        Note that this forward-running does not care about the fact that some of the inputs might be masked!
        In fact it **ignores** the mask and will simply compute unmasked outputs.
        Use predict() if you have masked inputs and want to "propagate" the mask appropriately.
        """

        # Normally we should go ahead and see if it fails, but in this particular case it's more helpful to test ahead:

        if inputs.ndim == 3:
            if inputs.shape[1] != self.ni:
                raise ValueError("Inputs with {ni} features (shape = {shape}) are not compatible with {me}".format(
                    ni=inputs.shape[1], shape=inputs.shape, me=str(self)))
        elif inputs.ndim == 2:
            if inputs.shape[0] != self.ni:
                raise ValueError("Inputs with {ni} features (shape = {shape}) are not compatible with {me}".format(
                    ni=inputs.shape[0], shape=inputs.shape, me=str(self)))

        output = inputs
        for l in self.layers:
            output = l.run(output)
        return output

    def predict(self, inputs):
        """
        We compute the outputs from the inputs using self.run, but here we do take care of the potential mask.

        This is never used during the training phase.

        :param inputs: a (potentially masked) 3D array

        :returns: a 3D array, appropriately masked

        """

        logger.info("Predicting with input = {intype} of shape {inshape}".format(
            intype=str(type(inputs)), inshape=str(inputs.shape)))

        if inputs.ndim != 3:
            raise ValueError("Sorry, I only accept 3D input")

        (inputs, outputsmask) = data.demask(inputs, no=self.no)

        # We can simply run the network with the unmasked inputs:

        logger.info("Running the actual predictions...")
        outputs = self.run(inputs)

        # And now mask these outputs, if required:

        if outputsmask is not None:
            outputs = np.ma.array(outputs, mask=outputsmask)

        return outputs
