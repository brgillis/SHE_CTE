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
File: python/SHE_CalibMLCore/plot.py

Created on: 09/05/17
Author: Malte Tewes & Thibault Kuntzer

Plots directly related to tenbilac objects

"""


import numpy as np
import itertools
import re
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.lines
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.cm
#from mpl_toolkits.axes_grid1 import make_axes_locatable

from . import err
from . import net


import logging
logger = logging.getLogger(__name__)


def errevo(ax, train, showtimes=True):
    """
    Plots the evolution of the error curve during the training, iteration per iteration.

    :param ax: matplotlib axes
    :param showtimes: If True, some training times are written on the curve, in minutes.
    """

    logger.info("Preparing error evolution plot for {train}".format(train=str(train)))
    # Preparint the data:
    optiterrs_train = np.array(train.optiterrs_train)
    optiterrs_val = np.array(train.optiterrs_val)
    optits = np.arange(len(train.optitparams))
    optbatchchangeits = getattr(train, "optbatchchangeits", [])

    # The cpu durations
    optittimes = np.array(train.optittimes)
    cumoptittimes = np.cumsum(optittimes) / 60.0  # minutes
    assert cumoptittimes.size == optittimes.size

    if optittimes.size > 10:
        labelstep = int(float(optittimes.size) / 10.0)
        timeindices = range(labelstep, optittimes.size, labelstep)
    else:
        timeindices = []

    for optbatchchangeit in optbatchchangeits:
        ax.axvline(optbatchchangeit, color="gray")

    ax.plot(optits, optiterrs_train, ls="-", color="black", label="Training batch")
    ax.plot(optits, optiterrs_val, ls="--", color="red", label="Validation set")

    if showtimes:
        for i in timeindices:
            ax.annotate("{0:.1f}".format(cumoptittimes[i]), xy=(
                optits[i], optiterrs_val[i]), xytext=(0, 10), textcoords='offset points')

    ax.set_yscale('log')
    ax.set_xlabel("Iteration")
    ax.set_xlim((optits[0], optits[-1]))
    ax.set_ylabel("Cost function value ({0})".format(train.errfctname))
    ax.legend()
    ax.set_title(train.title())


def paramsevo(ax, train):
    """
    Plots the evolution of the actual network parameters.

    """

    optits = np.arange(len(train.optitparams))

    # We have a normal Net:
    mynet = train.net
    optitparams = np.array(train.optitparams)

    paramlabels = mynet.get_paramlabels()
    multmode = False

    assert optitparams.shape[1] == mynet.nparams()
    for paramindex in range(mynet.nparams()):

        label = paramlabels[paramindex]
        if label.endswith("_bias"):
            ls = "--"
        elif label.endswith("_weight"):
            ls = "-"

        layermode = re.match("layer-(.*)-(.*)_(.*)", label).group(1)
        layername = re.match("layer-(.*)-(.*)_(.*)", label).group(2)

        if layername == "o":
            color = "black"
        elif layermode == "sum":
            color = "blue"
        elif layermode == "mult":
            multmode = True
            color = "green"
        else:
            raise ValueError("Layer mode {} unknown".format(layermode))

        pla = ax.plot(optits, optitparams[:, paramindex], ls=ls, color=color)

    ax.set_xlabel("Iteration")
    ax.set_xlim((optits[0], optits[-1]))
    ax.set_ylabel("Network parameter value")

    # Now creating the legend

    black_patch = matplotlib.patches.Patch(color='black', label='Output layer')

    line = matplotlib.lines.Line2D([], [], color='black', marker='', ls="-", label='Weight')
    dashed = matplotlib.lines.Line2D([], [], color='black', marker='', ls="--", label='Bias')

    handles = [line, dashed, black_patch]

    if multmode:
        red_patch = matplotlib.patches.Patch(color='blue', label='$\Sigma$ Hidden layers')
        green_patch = matplotlib.patches.Patch(color='green', label='$\Pi$ Hidden layers')
        handles += [red_patch]
        handles += [green_patch]
    else:
        red_patch = matplotlib.patches.Patch(color='blue', label='Hidden layers')
        handles += [red_patch]

    ax.legend(handles=handles)


def sumevo(train, filepath=None, showtimes=True):
    """
    Visualization of the evolution of the network parameters and error during the training,
    iteration per iteration

    :param showtimes: If True, some training times are written on the curve, in minutes.
    """

    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(2, 1, 1)
    errevo(ax, train, showtimes=showtimes)

    ax = plt.subplot(2, 1, 2)
    paramsevo(ax, train)

    plt.tight_layout()
    if filepath is None:
        plt.show()
    else:
        logger.info("Writing paramscurve to {0}".format(filepath))
        plt.savefig(filepath)
    plt.close()  # Important, otherwise it's still around for the next plt.show()


def biasevo(train, filepath=None):
    """
    Viz of the evolution of the individual bias terms at each iteration (massive!)
    """

    nsnaps = len(train.biassnaps_it)
    if nsnaps < 1:
        logger.warning("No snapshots to biasevo plot!")
        return

    # Members of this list have shape (outputs, cases). Now its (iterations, outputs, cases)
    trainbiases = np.array(train.biassnaps_train)
    # Members of this list have shape (outputs, cases). Now its (iterations, outputs, cases)
    valbiases = np.array(train.biassnaps_val)

    biases = np.dstack((trainbiases, valbiases))  # shape is (iteration, output-neuron, case)

    assert nsnaps == biases.shape[0]
    ncas = biases.shape[2]

    # For colors, we get the input data
    if train.dat is None:
        logger.warning("Need dat for biasevo plot!")
        return
    meantraininputs = np.mean(train.dat.fulltraininputs, axis=0)  # shape is (feature, case)
    meanvalinputs = np.mean(train.dat.valinputs, axis=0)  # shape is (feature, case)
    meaninputs = np.hstack((meantraininputs, meanvalinputs))  # shape is still (feature, case)

    assert meaninputs.shape[1] == biases.shape[2]

    logger.info("Preparing biasevo plot with {} cases and {} snapshots...".format(ncas, nsnaps))

    fig = plt.figure(figsize=(8 * train.net.no, 5 * train.net.ni))

    for io in range(train.net.no):
        for ii in range(train.net.ni):
            ax = plt.subplot(train.net.ni, train.net.no, io * train.net.ni + 1 + ii)

            # for optbatchchangeit in optbatchchangeits:
            #    ax.axvline(optbatchchangeit, color="gray", zorder=-20)
            ax.axhline(0.0, color="gray", zorder=-20)

            cmap = matplotlib.cm.get_cmap("jet")
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(
                vmin=np.min(meaninputs[ii, :]), vmax=np.max(meaninputs[ii, :])))

            for icas in range(ncas):
                plt.plot(train.biassnaps_it, biases[:, io, icas], color=sm.to_rgba(meaninputs[ii, icas]), alpha=0.1)

            ax.set_xlabel("Iteration (only snapshots at minibatch-changes are shown)")
            ax.set_ylabel("Bias on target '{}'".format(train.net.onames[io]))

            ax.set_yscale("symlog", linthreshy=1.e-3, linscaley=1)

            # fake up the array of the scalar mappable. Urgh...
            #sm._A = []
            sm.set_array([])
            cax = plt.colorbar(sm)
            cax.set_label("Feature '{}' (mean across reas)".format(train.net.inames[ii]))

    plt.tight_layout()
    if filepath is None:
        plt.show()
    else:
        logger.info("Writing biasevo to {0}".format(filepath))
        plt.savefig(filepath)
    plt.close()  # Important, otherwise it's still around for the next plt.show()


def outdistribs(train, filepath=None):
    """
    Viz of the output 
    """
    ncol = 7
    fig = plt.figure(figsize=(4 * ncol, 3.5 * train.net.no))

    dat = train.dat
    net = train.net
    assert dat is not None

    logger.info("Computing predictions...")
    trainoutputs = np.ma.array(net.run(dat.traininputs), mask=dat.trainoutputsmask)  # masked, 3D
    valoutputs = np.ma.array(net.run(dat.valinputs), mask=dat.valoutputsmask)
    logger.info("Done")

    trainerrors = trainoutputs - dat.traintargets  # 3D - 2D = 3D
    valerrors = valoutputs - dat.valtargets

    trainbiases = np.mean(trainerrors, axis=0)  # 2D (node, case)
    valbiases = np.mean(valerrors, axis=0)

    trainstds = np.std(trainoutputs, axis=0)  # 2D
    valstds = np.std(valoutputs, axis=0)

    valmsrbterms = err.msrb(valoutputs, dat.valtargets, rawterms=True)  # 2D
    trainmsrbterms = err.msrb(trainoutputs, dat.traintargets, rawterms=True)

    for io in range(train.net.no):

        # Subplots: (lines, columns, number)
        # We collect the stuff to be plotted as 1D arrays:

        # Warning, ravel does ignore the mask, so we use flatten.

        # The simple outputs
        thiso_valoutputs = valoutputs[:, io, :].flatten().compressed()
        thiso_trainoutputs = trainoutputs[:, io, :].flatten().compressed()
        ax = plt.subplot(train.net.no, ncol, (io * ncol) + 1)
        ax.hist(thiso_valoutputs, bins=50, histtype="step", color="red", label="Validation set")
        ax.hist(thiso_trainoutputs, bins=50, histtype="step", color="black", label="Training batch")
        ax.set_yscale('log')
        ax.set_ylabel("Counts (reas and cases)")
        ax.set_xlabel("Pred. output for '{0}'".format(net.onames[io]))

        # The targets
        thiso_valtargets = dat.valtargets[io, :].flatten()
        thiso_traintargets = dat.traintargets[io, :].flatten()
        ax = plt.subplot(train.net.no, ncol, (io * ncol) + 2)
        ax.hist(thiso_valtargets, bins=20, histtype="step", color="red", label="Validation set")
        ax.hist(thiso_traintargets, bins=20, histtype="step", color="black", label="Training batch")
        ax.set_ylabel("Counts (cases)")
        ax.set_xlabel("Targets '{0}'".format(net.onames[io]))

        # The prediction errors
        thiso_valerrors = valerrors[:, io, :].flatten().compressed()
        thiso_trainerrors = trainerrors[:, io, :].flatten().compressed()
        histrange = (-1.5, 1.5)
        ax = plt.subplot(train.net.no, ncol, (io * ncol) + 3)
        ax.hist(thiso_valerrors, bins=50, histtype="step", color="red", label="Validation set", range=histrange)
        ax.hist(thiso_trainerrors, bins=50, histtype="step", color="black", label="Training batch", range=histrange)
        ax.set_ylabel("Counts (reas and cases)")
        ax.set_xlabel("Errors of pred. '{0}'".format(net.onames[io]))
        ax.set_yscale('log')
        ax.set_xlim(histrange)

        # The biases
        thiso_valbiases = valbiases[io, :].flatten().compressed()
        thiso_trainbiases = trainbiases[io, :].flatten().compressed()
        ax = plt.subplot(train.net.no, ncol, (io * ncol) + 4)
        ax.hist(thiso_valbiases, bins=50, histtype="step", color="red", label="Validation set")
        ax.hist(thiso_trainbiases, bins=50, histtype="step", color="black", label="Training batch")
        ax.set_yscale('log')
        ax.set_ylabel("Counts (cases)")
        ax.set_xlabel("Bias of pred. for '{0}'".format(net.onames[io]))

        # The stds
        thiso_valstds = valstds[io, :].flatten().compressed()
        thiso_trainstds = trainstds[io, :].flatten().compressed()
        ax = plt.subplot(train.net.no, ncol, (io * ncol) + 5)
        ax.hist(thiso_valstds, bins=50, histtype="step", color="red", label="Validation set")
        ax.hist(thiso_trainstds, bins=50, histtype="step", color="black", label="Training batch")
        ax.set_yscale('log')
        ax.set_ylabel("Counts (cases)")
        ax.set_xlabel("STD of pred. for '{0}'".format(net.onames[io]))

        # Against each other
        ax = plt.subplot(train.net.no, ncol, (io * ncol) + 6)
        assert thiso_trainstds.size == thiso_trainbiases.size  # Indeed even if they were masked, their masks should be identical
        assert thiso_valstds.size == thiso_valbiases.size
        ax.plot(thiso_valbiases, thiso_valstds, marker=".", ms=2, ls="None",
                color="red", label="Validation set", rasterized=True)
        ax.plot(thiso_trainbiases, thiso_trainstds, marker=".", ms=2, ls="None",
                color="black", label="Training batch", rasterized=True)
        ax.set_ylabel("STD of pred. for '{0}'".format(net.onames[io]))
        ax.set_xlabel("Bias of pred. for '{0}'".format(net.onames[io]))

        # The MSRB terms:
        thiso_valmsrbterms = valmsrbterms[io, :].flatten().compressed()
        thiso_trainmsrbterms = trainmsrbterms[io, :].flatten().compressed()
        ax = plt.subplot(train.net.no, ncol, (io * ncol) + 7)
        ax.hist(thiso_valmsrbterms, bins=50, histtype="step", color="red", label="Validation set")
        ax.hist(thiso_trainmsrbterms, bins=50, histtype="step", color="black", label="Training batch")
        ax.set_ylabel("Counts (cases)")
        ax.set_xlabel("Relative biases of '{0}'".format(net.onames[io]))
        try:
            ax.set_yscale('log')
        except:
            logger.warning("Probably you have only one rea, this plot failed.")
            plt.close()
            return()

    plt.tight_layout()
    if filepath is None:
        plt.show()
    else:
        logger.info("Writing outdistribs to {0}".format(filepath))
        plt.savefig(filepath)
    plt.close()  # Important, otherwise it's still around for the next plt.show()


def errorinputs(train, filepath=None, io=0):
    """
    Viz of the prediction errors as function of inputs

    :param io: the index of the ouput I should use. If you have only one neuron, this is 0.

    """
    nlines = 3

    fig = plt.figure(figsize=(3.3 * train.net.ni, 3.5 * nlines))
    plt.figtext(0.5, 1.0, train.title(), ha="center", va="top")

    dat = train.dat
    net = train.net
    assert dat is not None

    logger.info("Computing predictions...")
    trainoutputs = np.ma.array(net.run(dat.traininputs), mask=dat.trainoutputsmask)
    valoutputs = np.ma.array(net.run(dat.valinputs), mask=dat.valoutputsmask)
    logger.info("Done")

    trainerrors = trainoutputs - dat.traintargets  # 3D - 2D = 3D: (rea, label, case)
    valerrors = valoutputs - dat.valtargets  # idem

    #valmsrbterms = err.msrb(valoutputs, dat.valtargets, rawterms=True)
    #trainmsrbterms = err.msrb(trainoutputs, dat.traintargets, rawterms=True)

    for ii in range(train.net.ni):

        ax = plt.subplot(nlines, train.net.ni, ii + 1)

        ax.hist(np.ravel(dat.valinputs[:, ii, :]), bins=50, histtype="step", color="red", label="Validation set")
        ax.hist(np.ravel(dat.traininputs[:, ii, :]), bins=50, histtype="step", color="black", label="Training batch")
        ax.set_yscale('log')
        if ii == 0:
            ax.set_ylabel("Counts (reas + cases, masked = 0)")
        # if ii == int(train.net.ni/2):
        #    ax.set_title(train.title())
        #ax.set_xlabel("Input '{0}'".format(net.inames[ii]))
        ax.set_xticklabels([])  # Hide x tick labels
        ax.set_xlim(-1.0, 1.0)
        # in this plot the masked inputs are shown with a value of 0, this is desired and expected.

        # Now we viz the biases and stds of the predictions (along the reas)

        # The biases and stds on the validation set:
        valbiases = np.mean(valerrors, axis=0)[io, :]  # 1D (case)
        valstds = np.std(valoutputs, axis=0)[io, :]
        # Inflated to 2D (rea, case), with all reas having the same value.
        valbiases = np.tile(valbiases, (dat.valinputs.shape[0], 1))
        valstds = np.tile(valstds, (dat.valinputs.shape[0], 1))

        # We want to get a version of these particular inputs with a mask, selecting only one node
        valinputs = dat.valinputs[:, ii:ii + 1, :]  # 3D : (rea, 1, case)
        # Here an input is masked if any other input of that rea was masked.
        valinputs = np.ma.array(valinputs, mask=dat.valoutputsmask)
        assert valinputs.shape[1] == 1
        valinputs = valinputs[:, 0, :]  # 2D : (rea, case)

        assert valbiases.shape == valinputs.shape
        assert valstds.shape == valinputs.shape

        # For the plot, we only use points for which the inputs are unmasked
        valbiases = valbiases.flatten()
        valstds = valstds.flatten()
        valinputs = valinputs.flatten()
        valbiases = np.ma.array(valbiases, mask=valinputs.mask)  # copying the mask
        valstds = np.ma.array(valstds, mask=valinputs.mask)
        valbiases = valbiases.compressed()  # Getting rid of the masked elements
        valstds = valstds.compressed()
        valinputs = valinputs.compressed()

        assert valinputs.size == valbiases.size
        assert valinputs.size == valstds.size

        ax = plt.subplot(nlines, train.net.ni, train.net.ni + ii + 1)
        ax.plot(valinputs, valbiases, marker=".", color="red", ls="None", ms=1)
        ax.axhline(0.0, color="black", lw=1, ls="--")
        if ii == 0:
            ax.set_ylabel("Bias of pred. for '{0}'".format(net.onames[io]))
        else:
            ax.set_yticklabels([])  # Hide y tick labels
        #ax.set_xlabel("Input '{0}'".format(net.inames[ii]))
        ax.set_xlim(-1.0, 1.0)
        ax.set_xticklabels([])  # Hide x tick labels

        ax = plt.subplot(nlines, train.net.ni, 2 * train.net.ni + ii + 1)
        ax.plot(valinputs, valstds, marker=".", color="red", ls="None", ms=1)
        if ii == 0:
            ax.set_ylabel("STD of pred. for '{0}'".format(net.onames[io]))
        else:
            ax.set_yticklabels([])  # Hide y tick labels
        ax.set_xlabel("Input '{0}'".format(net.inames[ii]))
        ax.set_xlim(-1.0, 1.0)

    plt.tight_layout()
    if filepath is None:
        plt.show()
    else:
        logger.info("Writing errorinputs to {0}".format(filepath))
        plt.savefig(filepath)
    plt.close()  # Important, otherwise it's still around for the next plt.show()


def draw_link(ax, start, end, **kwargs):
    """
    Computes a Bezier curve between two points (`start` and `end`) on an axis `ax`

    :param ax: the axis to draw on
    :param start: The starting point, must be an array [Sx, Sy]
    :param end: The end point of the curve, must be an array [Ex, Ey]

    Any additional kawrgs are directly passed to `patches.PathPatch` to control the line style.
    """

    verts = [
        (start[0], start[1]),  # P0
        (start[0] + 0.5, start[1]),  # P1
        (end[0] - 0.5, end[1]),  # P2
        (end[0], end[1]),  # P3
    ]

    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4,
             ]

    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='None', **kwargs)
    ax.add_patch(patch)


def scale_bias(b, scale=75.):
    """
    Returns the scaled biased for the `plot.netviz` plot.

    :param b: the bias to scale
    :param scale: the scale parameter to be applied on the absolute value of `b`. Default=20
    """
    return scale * np.abs(b)


def scale_weight(w, scale=2.):
    """
    Returns the scaled weight for the `plot.netviz` plot.

    :param w: the weight to scale
    :param scale: the scale parameter to be applied on the absolute value of `w`. Default=2
    """
    return scale * np.abs(w)


def get_color(v, pos="orange", neg="navy"):
    """
    Returns a color code for the scalar value `v`.

    :param v: scalar value
    :param pos: name or color-code for the positive values
    :param neg: name or color-code for the negative values
    """
    if v > 0:
        c = pos
    else:
        c = neg
    return c


def get_symbol(mode, latex=True):
    """
    Returns the symbol used for displaying the layer mode
    """
    if latex:
        if mode == "sum":
            return r"$\sum$"
        elif mode == "mult":
            return r"$\prod$"
        else:
            raise NotImplemented("mode not recognised")
    else:
        if mode == "sum":
            return "\Sigma"
        elif mode == "mult":
            return "\Pi"
        else:
            raise NotImplemented("mode not recognised")


def netviz(net, title="", legend=True, filepath=None):
    """
    Draws a visualisation of the network in the style of the `Tensorflow` playground.

    :param net: the network
    :param title: The string to display as title.
    :param legend: Whether to show the legend for the tickness of the lines and points.
    :param filepath: The file path to save the data to. If `None` (default) shows the figures.
    """

    nmax = np.amax([net.ni, net.no, np.amax(net.nhs)]) * 1.

    fig = plt.figure()
    ax = fig.add_subplot(111)

    write_mode = False
    for l in net.layers:
        if not l.mode == "sum":
            write_mode = True

    biasm = "^"

    plt.title(title)

    plt_kwargs = {'marker': 's', 's': 35, 'c': 'k', 'zorder': 1}

    for ii, l in enumerate(net.layers):
        dy = nmax / 2 - (l.ni * 1.) / 2
        lnis = np.arange(l.ni) + dy

        plt.scatter(np.zeros_like(lnis) + ii, lnis, **plt_kwargs)

        if ii == 0:
            # Write the name of the features
            inames = [r"$\mathrm{{{}}}$".format(i) for i in net.inames]
            for iii, inp in enumerate(inames):
                plt.text(ii - 0.08, iii + dy, inp.replace("_", "\_"),
                         horizontalalignment='right', verticalalignment='center')

        if ii >= len(net.layers) - 1:
            flnis = np.arange(net.no)
            dyf = (nmax - (net.no * 1.)) * 0.5
        else:
            flnis = np.arange(net.layers[ii + 1].ni)
            dyf = (nmax - (net.layers[ii + 1].ni * 1.)) * 0.5

        for iw, w in enumerate(l.weights):
            # Draw the weights
            for il, link in enumerate(w):
                draw_link(ax, start=[ii, lnis[il]],
                          end=[ii + 1, flnis[iw] + dyf], lw=scale_weight(link), edgecolor=get_color(link), zorder=-1)

            # Draw the biases
            plt.scatter([ii + 1.], [flnis[iw] + dyf + 0.12], c=get_color(l.biases[iw]),
                        edgecolors="None", s=scale_bias(l.biases[iw]), marker=biasm)

        if write_mode:
            plt.text(ii + 1., flnis[0] + dyf - 0.25, get_symbol(l.mode), horizontalalignment='center', size=7)

    # Draw output
    nos = np.arange(net.no)
    plt.scatter(np.zeros_like(nos) + ii + 1, nos + dyf, **plt_kwargs)

    # Name the output
    onames = [r"$\mathrm{{{}}}$".format(i) for i in net.onames]
    for iii, inp in enumerate(onames):
        plt.annotate(inp.replace("_", "\_"), xy=(ii + 1.08, iii + dyf),
                     horizontalalignment='left', verticalalignment='center')

    # Draw legend
    if legend:
        ws = [-1., -0.5, 0.5, 1.]
        for iw, w in enumerate(ws):
            yy = nos[-1] + dyf - 1. - iw * 0.2
            plt.annotate(r"$%1.1f$" % w, xy=(ii + 0.95, yy), horizontalalignment='right', verticalalignment='center')
            plt.scatter([ii + 1], [yy], c=get_color(w), edgecolors="None", s=scale_bias(w), marker=biasm)
            plt.plot([ii + 1.07, ii + 1.47], [yy, yy], c=get_color(w), lw=scale_weight(w))

    # Taking care of a few things
    plt.xlim([-0.23 * (len(net.nhs) + 2), ii + 1.7])
    plt.tight_layout()
    plt.axis('off')
    plt.gca().invert_yaxis()

    if filepath is None:
        plt.show()
    else:
        logger.info("Writing netviz to {0}".format(filepath))
        plt.savefig(filepath)
    plt.close()  # Important, otherwise it's still around for the next plt.show()


def summaryerrevo(committee, filepath=None, ax=None):
    """
    First take at plotting the error curves of committee members, to compare their performances.

    :param committee: a list of Training objects
    """

    if ax is None:
        fig = plt.figure(figsize=(14, 10))
        ax = plt.subplot(1, 1, 1)
    else:
        fig = None

    # We sort the committee:
    committee = sorted(committee, key=lambda trainobj: trainobj.optiterrs_train[-1], reverse=True)

    logger.info("Preparing summary plot with {} members...".format(len(committee)))
    coloriter = iter(plt.cm.jet(np.linspace(0, 1, len(committee))))

    for trainobj in committee:

        # Preparint the data:
        optiterrs_train = np.array(trainobj.optiterrs_train)
        optiterrs_val = np.array(trainobj.optiterrs_val)
        optits = np.arange(len(trainobj.optitparams))
        optbatchchangeits = getattr(trainobj, "optbatchchangeits", [])

        trainerr = trainobj.optiterrs_train[-1]
        valerr = trainobj.optiterrs_val[-1]
        valerrratio = valerr / trainerr

        color = next(coloriter)
        ax.plot(optits, optiterrs_train, ls="-", color=color, label="'{}': {:.2e} ({:.1f})".format(
            trainobj.name, trainerr, valerrratio))
        ax.plot(optits, optiterrs_val, ls="--", color=color)  # No label for the validation

    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_xlabel("Iteration")
    #ax.set_xlim((optits[0], optits[-1]))
    ax.set_ylabel("Cost function value")
    ax.legend()
    # ax.set_title(train.title())

    if fig is not None:
        plt.tight_layout()
        if filepath is None:
            plt.show()
        else:
            logger.info("Writing summaryerrevo to {}".format(filepath))
            plt.savefig(filepath)
        plt.close()  # Important, otherwise it's still around for the next plt.show()
