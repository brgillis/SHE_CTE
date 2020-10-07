""" @file plot_bias_measurements.py

    Created 26 Apr 2017

    Main function to plot bias measurements.
"""

__updated__ = "2020-10-07"

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
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

from copy import deepcopy
from os.path import join

from SHE_PPT import products
from SHE_PPT.file_io import read_xml_product
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import fsolve

import matplotlib.pyplot as pyplot
import numpy as np


psf_gal_size_ratio = (0.2 / 0.3)**2


def Spline(*args, **kwargs):
    return InterpolatedUnivariateSpline(*args, k=1, **kwargs)


def get_spline_slope(x_vals, y_vals, x):
    spline = Spline(x_vals, y_vals)
    return spline.derivative()(x)


def jackknife_resampling(data):
    """jackknife_resampling, copied from astropy v4.0.1 and modified to work with 2D data"""

    n = data.shape[0]
    if n <= 0:
        raise ValueError("data must contain at least one measurement.")

    len_data = len(data[0])

    resamples = np.empty([n, n - 1, len_data])

    for i in range(n):
        indices_to_del = np.arange(len_data) + len_data * i
        resamples[i] = np.delete(data, indices_to_del).reshape(n - 1, len_data)

    return resamples


def jackknife_stats(data, statistic, confidence_level=0.95):
    """jackknife_stats, copied from astropy v4.0.1"""

    if not (0 < confidence_level < 1):
        raise ValueError("confidence level must be in (0, 1).")

    # make sure original data is proper
    n = data.shape[0]
    if n <= 0:
        raise ValueError("data must contain at least one measurement.")

    # Only import scipy if inputs are valid
    from scipy.special import erfinv

    resamples = jackknife_resampling(data)

    stat_data = statistic(data)
    jack_stat = np.apply_along_axis(statistic, 1, resamples)
    mean_jack_stat = np.mean(jack_stat, axis=0)

    # jackknife bias
    bias = (n - 1) * (mean_jack_stat - stat_data)

    # jackknife standard error
    std_err = np.sqrt((n - 1) * np.mean((jack_stat - mean_jack_stat) * (jack_stat -
                                                                        mean_jack_stat), axis=0))

    # bias-corrected "jackknifed estimate"
    estimate = stat_data - bias

    z_score = np.sqrt(2.0) * erfinv(confidence_level)
    conf_interval = estimate + z_score * np.array((-std_err, std_err))

    return estimate, bias, std_err, conf_interval


testing_data_labels = {"S": "Sky Level (ADU/pixel)",
                       "E": r"$\sigma(e)$",
                       "T": r"Truncation Radius (multiple of $r_{\rm s}$)"}

titles = {"S": "Varying Sky Background Level",
          "E": "Varying Galaxy Ellipticity Distribution Sigma",
          "T": "Varying Disk Truncation Radius", }

testing_data_labels_no_units = {"S": "Sky Level",
                                "E": r"$\sigma(e)$",
                                "T": "Truncation Radius"}

tag_template = "Ep0Pp0Sp0Tp0"

testing_variant_labels = ("m2", "m1", "p0", "p1", "p2")

measurement_key_templates = ("mDIM", "mDIM_err", "cDIM", "cDIM_err")
measurement_colors = {"m": "r", "c": "b"}

centre_index = 2

x_values = {"S": [45.52, 45.61, 45.71, 45.80,  45.90],
            "E": [0.18556590758327751, 0.21333834512347458, 0.2422044791810781,
                  0.27099491059570091, 0.30241731263684996],
            "T": [3.0, 4.0, 4.5, 6.0, 10.0],
            }

x_ranges = {"S": [45.50, 45.92],
            "E": [0.170, 0.315],
            "T": [2.5, 10.5]}

target_limit_factors = {"P": {"m": 64, "c": 32},
                        "S": {"m": 32, "c": 20},
                        "E": {"m": 24,  "c": 32},
                        "T": {"m": 24,  "c": 32}, }

y_range = (1e-6, 5e-1)

err_factor = 1

m_target = 1e-4
c_target = 5e-6

method_colors = {"KSB": "k",
                 "KSB_big": (0.5, 0.5, 0.5),
                 "REGAUSS": "r",
                 "REGAUSS_big": (1.0, 0.5, 0.5),
                 "LensMC": "b",
                 "MomentsML": "m",
                 "BFD": "g", }
method_offsets = {"KSB": 0,
                  "KSB_big": 0.005,
                  "REGAUSS": 0.01,
                  "REGAUSS_big": 0.015,
                  "LensMC": -0.01,
                  "MomentsML": -0.02,
                  "BFD": 0.02, }

target_labels = {"base": r"$0.05\times$",
                 "high": r"$1.00\times$", }
target_shapes = {"base": (3, 2, 0),
                 "high": (3, 2, 180), }

calibration_labels = ["", "_normed"]

fontsize = 12
text_size = 18


def plot_bias_measurements_from_args(args):
    """ @TODO main docstring
    """

    # Determine the qualified path to the root data folder
    if args.root_data_folder is None:
        root_data_folder = args.workdir
    elif args.root_data_folder[0] == "/":
        root_data_folder = args.root_data_folder
    else:
        # Relative to workdir in this case
        root_data_folder = join(args.workdir, args.root_data_folder)

    # Open and and keep in memory all bias measurements
    all_bias_measurements = {}

    def read_bias_measurements(tag):
        if not tag in all_bias_measurements:
            all_bias_measurements[tag] = read_xml_product(
                join(root_data_folder, args.bias_measurements_head + tag + ".xml"), workdir=root_data_folder)

    # Do a loop of reading for each property
    for testing_variant in testing_variant_labels:

        e_tag = tag_template.replace("Ep0", "E" + testing_variant)
        read_bias_measurements(e_tag)

        s_tag = tag_template.replace("Sp0", "S" + testing_variant)
        read_bias_measurements(s_tag)

        t_tag = tag_template.replace("Tp0", "T" + testing_variant)
        read_bias_measurements(t_tag)

    # Plot the biases and errors for each measurement
    for testing_data_key in testing_data_labels:

        fractional_limits = {}

        for measurement_key_template in measurement_key_templates:

            if not args.plot_error and "_err" in measurement_key_template:
                continue

            all_methods_data = {}

            for calibration_label in calibration_labels:

                # Don't do a normed plot for errors or if disabled
                if ("_err" in measurement_key_template or args.unnormed_only) and calibration_label == "_normed":
                    continue

                # We can't fully skip unnormed plots, but we still won't plot them
                do_fig = (not args.normed_only) or calibration_label == "_normed"

                # Plot regularly for dim = 0
                measurement_key = measurement_key_template.replace("DIM", "")
                measurement_key_1 = measurement_key_template.replace("DIM", "1")
                measurement_key_2 = measurement_key_template.replace("DIM", "2")

                # Set up the figure
                if do_fig:
                    fig = pyplot.figure()
                    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

                    ax = fig.add_subplot(1, 1, 1)
                    ax.set_xlabel(testing_data_labels[testing_data_key], fontsize=fontsize)
                    if calibration_label == "_normed":
                        ax.set_ylabel(r"$\Delta " + measurement_key + "$ (relative to " +
                                      testing_data_labels_no_units[testing_data_key] +
                                      "=%2.2f)" % x_values[testing_data_key][2], fontsize=fontsize)
                    else:
                        ax.set_ylabel("$" + measurement_key.replace("_err", r"_{\rm err}") + "$", fontsize=fontsize)

                xlim = deepcopy(x_ranges[testing_data_key])

                # Plot points for each method
                for method in args.methods:

                    lx = []
                    ly1 = []
                    ly2 = []
                    ly1_o = []
                    ly2_o = []

                    # Loop over each testing variant (m2, m1, etc.)
                    for i in range(len(testing_variant_labels)):

                        # Get the tag for this variant
                        testing_variant = testing_variant_labels[i]
                        tag = tag_template.replace(testing_data_key + "p0", testing_data_key + testing_variant)

                        # Set the x value
                        lx.append(x_values[testing_data_key][i] + method_offsets[method] *
                                  np.abs(x_values[testing_data_key][3] - x_values[testing_data_key][1]))

                        # Get the bias measurements for this method and testing variant
                        g1_bias_measurements, g2_bias_measurements = all_bias_measurements[tag].get_method_bias_measurements(
                            method, workdir=root_data_folder)
                        g1_bias_measurement = getattr(g1_bias_measurements, measurement_key)
                        g2_bias_measurement = getattr(g2_bias_measurements, measurement_key)

                        if g1_bias_measurement == '' or g1_bias_measurement == 'NaN':
                            g1_bias_measurement = np.nan
                        if g2_bias_measurement == '' or g2_bias_measurement == 'NaN':
                            g2_bias_measurement = np.nan

                        # If we're norming, correct bias measurements by the central value
                        if calibration_label == "_normed":
                            g1_central_bias_measurements, g2_central_bias_measurements = (
                                all_bias_measurements[tag_template].get_method_bias_measurements(method, workdir=root_data_folder))
                            g1_central_bias_measurement = getattr(g1_central_bias_measurements, measurement_key)
                            g2_central_bias_measurement = getattr(g2_central_bias_measurements, measurement_key)

                            # Correct differently for m and c
                            if "m" in measurement_key:
                                g1_bias_measurement = (1 + g1_bias_measurement) / (1 + g1_central_bias_measurement) - 1
                                g2_bias_measurement = (1 + g2_bias_measurement) / (1 + g2_central_bias_measurement) - 1
                            elif "c" in measurement_key:
                                g1_bias_measurement -= g1_central_bias_measurement
                                g2_bias_measurement -= g2_central_bias_measurement

                        ly1.append(g1_bias_measurement)
                        ly2.append(g2_bias_measurement)

                        # To calculate combined error, we also need non-error
                        if "_err" in measurement_key:
                            g1_o = getattr(g1_bias_measurements, measurement_key.replace("_err", ""))
                            g2_o = getattr(g2_bias_measurements, measurement_key.replace("_err", ""))
                        else:
                            # And for non-error values, we'll want to plot error bars too, so we need that
                            g1_o = getattr(g1_bias_measurements, measurement_key + "_err")
                            g2_o = getattr(g2_bias_measurements, measurement_key + "_err")

                        if g1_o != '' or g1_o == 'NaN':
                            ly1_o.append(g1_o)
                        else:
                            ly1_o.append(np.nan)
                        if g2_o != '' or g1_o == 'NaN':
                            ly2_o.append(g2_o)
                        else:
                            ly2_o.append(np.nan)

                    x_vals = np.array(lx, dtype=float)
                    y1_vals = np.array(ly1, dtype=float)
                    y2_vals = np.array(ly2, dtype=float)
                    y1_o_vals = np.array(ly1_o, dtype=float)
                    y2_o_vals = np.array(ly2_o, dtype=float)

                    # Determine combined y differently for error and non-error values
                    if "_err" in measurement_key:
                        y_vals = (np.abs(y1_o_vals) * y1_vals + np.abs(y2_o_vals)
                                  * y2_vals) / np.sqrt(y1_o_vals**2 + y2_o_vals**2)
                        y1_errs = None
                        y2_errs = None
                        y_errs = None
                    elif calibration_label == "_normed":
                        y_vals = np.sqrt(y1_vals**2 + y2_vals**2)
                        # Carry over errors from previous run, on unnormed data
                    else:
                        y_vals = np.sqrt(y1_vals**2 + y2_vals**2)
                        y1_errs = y1_o_vals
                        y2_errs = y2_o_vals
                        y_errs = (np.abs(y1_vals) * y1_errs + np.abs(y2_vals)
                                  * y2_errs) / np.sqrt(y1_vals**2 + y2_vals**2)

                    # Plot the values (and optionally error bars)
                    if not "_err" in measurement_key:
                        if do_fig:
                            ax.errorbar(x_vals, y_vals, y_errs, color=method_colors[method], linestyle='None')
                    else:
                        y1_vals *= err_factor
                        y2_vals *= err_factor
                        y_vals *= err_factor
                    if do_fig:
                        ax.plot(x_vals, y_vals, color=method_colors[method], marker='o')

                    if "_big" not in method:
                        label = method
                    else:
                        label = "Big " + method.replace("_big", "")

                    # Save this data for the next plot if not showing errors
                    if not "_err" in measurement_key:
                        method_data = {"x": x_vals,
                                       "y": y_vals,
                                       "y_err": y_errs,
                                       "y1": y1_vals,
                                       "y1_err": y1_errs,
                                       "y2": y2_vals,
                                       "y2_err": y2_errs,
                                       "y1_o": y1_o_vals,
                                       "y2_o": y2_o_vals,
                                       }

                        # Calculate y/y1/y2 slope at the central value and jackknife error of it
                        central_x = x_vals[centre_index]

                        def get_spline_slope_at_central_x(xy_vals):
                            x_vals, y_vals = zip(*xy_vals)
                            return get_spline_slope(x_vals, y_vals, central_x)

                        for chosen_y_vals, y_label in ((y_vals, "y_slope"),
                                                       (y1_vals, "y1_slope"),
                                                       (y2_vals, "y2_slope")):
                            xy_vals = np.array(list(zip(x_vals, chosen_y_vals)))
                            estimate, _, stderr, conf_interval = jackknife_stats(
                                xy_vals, get_spline_slope_at_central_x, 0.95)
                            method_data[y_label] = estimate
                            method_data[y_label + "_err"] = stderr
                            method_data[y_label + "_conf_interval"] = conf_interval

                        all_methods_data[(method, calibration_label)] = method_data

                # Plot the target line
                if "m" in measurement_key:
                    target = m_target
                else:
                    target = c_target

                if do_fig:

                    pyplot.title(titles[testing_data_key])

                    ax.plot(xlim, [target, target], label=None, color="k", linestyle="dashed")
                    ax.plot(xlim, [20 * target, 20 * target], label=None, color="k", linestyle="dotted")
                    ax.plot(xlim, [0, 0], label=None, color="k", linestyle="solid")

                    # Set the limits and scale
                    ax.set_xlim(xlim)
                    # ax.set_ylim(y_range) # FIXME - uncomment once we know about what the range will be
                    if not calibration_label == "_normed":
                        ax.set_yscale("log", nonposy="clip")

                    # Show the legend
                    ax.legend(loc="lower right", numpoints=1)

                    # Save and show it
                    output_filename = join(args.workdir, args.output_file_name_head + "_" +
                                           testing_data_key + "_" + measurement_key + calibration_label + "." +
                                           args.output_format)
                    pyplot.savefig(output_filename, format=args.output_format, bbox_inches="tight", pad_inches=0.05)
                    if not args.hide:
                        fig.show()
                    else:
                        pyplot.close()

                    # Plot individual indices as well
                    for index in (1, 2):

                        # Set up the figure
                        fig = pyplot.figure()

                        pyplot.title(titles[testing_data_key])

                        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

                        ax = fig.add_subplot(1, 1, 1)
                        ax.set_xlabel(testing_data_labels[testing_data_key], fontsize=fontsize)
                        if calibration_label == "_normed":
                            ax.set_ylabel(r"$\Delta " + measurement_key + r"_{\rm " + str(index) +
                                          "}$ (relative to " + testing_data_labels_no_units[testing_data_key] +
                                          "=%2.2f)" % x_values[testing_data_key][2], fontsize=fontsize)
                        else:
                            ax.set_ylabel("$" + measurement_key + "_" + str(index) + "$", fontsize=fontsize)

                        # Plot the values and error bars

                        x_spline_vals = np.linspace(x_vals[0], x_vals[-1], 100)

                        for method in args.methods:
                            ax.plot(x_vals,
                                    all_methods_data[(method, calibration_label)]["y" + str(index)],
                                    color=method_colors[method], marker='o', label=method)

                            ax.errorbar(x_vals,
                                        all_methods_data[(method, calibration_label)]["y" + str(index)],
                                        all_methods_data[(method, calibration_label)]["y" + str(index) + "_err"],
                                        color=method_colors[method], linestyle='None')

                        # Plot zero and limits
                        xlim = deepcopy(ax.get_xlim())
                        ax.plot(xlim, [target, target], label=None, color="k", linestyle="dashed")
                        ax.plot(xlim, [20 * target, 20 * target], label=None, color="k", linestyle="dotted")
                        ax.plot(xlim, [-target, -target], label=None, color="k", linestyle="dashed")
                        ax.plot(xlim, [-20 * target, -20 * target], label=None, color="k", linestyle="dotted")
                        ax.plot(xlim, [0, 0], label=None, color="k", linestyle="solid")
                        ax.set_xlim(xlim)

                        # Save and show it

                        ax.legend(loc="lower right", numpoints=1)
                        output_filename = join(args.workdir, args.output_file_name_head + "_" +
                                               testing_data_key + "_" + measurement_key + str(index) +
                                               calibration_label + "." + args.output_format)
                        pyplot.savefig(output_filename, format=args.output_format,
                                       bbox_inches="tight", pad_inches=0.05)
                        if not args.hide:
                            fig.show()
                        else:
                            pyplot.close()

            # Now plot dim 1 v dim 2

            if not args.unnormed_only:

                # Set up the figure
                fig = pyplot.figure()

                pyplot.title(titles[testing_data_key])

                fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

                ax = fig.add_subplot(1, 1, 1)
                ax.set_xlabel(r"$\Delta " + measurement_key_1.replace("_err", r"_{\rm err}") + "$", fontsize=fontsize)
                ax.set_ylabel(r"$\Delta " + measurement_key_2.replace("_err", r"_{\rm err}") + "$", fontsize=fontsize)

                if("m" in measurement_key):
                    base_target = m_target
                else:
                    base_target = c_target

                target_factor = target_limit_factors[testing_data_key][measurement_key.replace("_err", "")]
                xlim = [-1.1 * base_target, 1.1 * base_target]
                ylim = [-1.1 * base_target, 1.1 * base_target]

                # Plot points for each method
                for method in args.methods:

                    # Get the data we saved before
                    method_data = all_methods_data[(method, "")]

                    ax.plot(method_data["y1"] - method_data["y1"][2], method_data["y2"] - method_data["y2"][2],
                            color=method_colors[method], marker='o', linestyle='None')

                    xvals = method_data["y1"] - method_data["y1"][2]
                    yvals = method_data["y2"] - method_data["y2"][2]

                    if 1.1 * np.abs(xvals).max() > xlim[1]:
                        xlim[0] = -1.1 * np.abs(xvals).max()
                        xlim[1] = 1.1 * np.abs(xvals).max()
                    if 1.1 * np.abs(yvals).max() > ylim[1]:
                        ylim[0] = -1.1 * np.abs(yvals).max()
                        ylim[1] = 1.1 * np.abs(yvals).max()

                    # Calculate and plot an interpolating spline
                    y1_spline = Spline(method_data["x"], xvals)
                    y2_spline = Spline(method_data["x"], yvals)

                    if "_err" not in measurement_key:
                        def y_spline(x): return np.sqrt(y1_spline(x)**2 + y2_spline(x)**2)
                    else:
                        y1_o_spline = Spline(method_data["x"], method_data["y1_o"] - method_data["y1_o"][2])
                        y2_o_spline = Spline(method_data["x"], method_data["y2_o"] - method_data["y2_o"][2])

                        def y_spline(x): return (np.abs(y1_spline(x)) * y1_o_spline(x) + np.abs(y2_spline(x))
                                                 * y1_o_spline(x)) / np.sqrt(y1_o_spline(x)**2 + y2_o_spline(x)**2)

                    x_spline_vals = np.linspace(x_vals[0], x_vals[-1], 100)

                    y1_spline_vals = y1_spline(x_spline_vals)
                    y2_spline_vals = y2_spline(x_spline_vals)

                    if "_big" not in method:
                        label = method
                    else:
                        label = "Big " + method.replace("_big", "")

                    ax.plot(y1_spline_vals, y2_spline_vals, color=method_colors[method], marker='None',
                            label=label)

                    if args.plot_fractional_limits and "_err" not in measurement_key:

                        # Now try to solve for where it intersects the target lines
                        limit_label_base = method + "_" + measurement_key

                        for (target, target_key) in ((base_target, "base"), (20 * base_target, "high")):

                            limit_label = limit_label_base + "_" + target_key
                            intersections = {}

                            for (i, side_label) in ((1, "low"), (3, "high")):
                                guess = method_data["x"][2] + target / (method_data["y"][i] -
                                                                        method_data["y"][2]) * (method_data["x"][i] -
                                                                                                method_data["x"][2])
                                intersections[side_label] = fsolve(lambda x: y_spline(x) - target, guess)

                            fractional_limits[limit_label] = (
                                (intersections["high"] - intersections["low"]) / (2. * method_data["x"][2]))[0]

                        print(("Fraction limit on " + testing_data_labels_no_units[testing_data_key] + " for method " +
                               method + " for " + measurement_key + ": " +
                               str(fractional_limits[limit_label_base + "_base"]) + ",\t" +
                               str(fractional_limits[limit_label_base + "_high"])))

                theta_vals = np.linspace(0, 2 * np.pi, 360)

                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

                ax.plot(base_target * np.cos(theta_vals), base_target *
                        np.sin(theta_vals), label=None, color="k", linestyle="dashed",)
                ax.plot(20 * base_target * np.cos(theta_vals), 20 * base_target *
                        np.sin(theta_vals), label=None, color="k", linestyle="dotted")
                ax.plot(xlim, [0, 0], label=None, color="k", linestyle="solid")
                ax.plot([0, 0], ylim, label=None, color="k", linestyle="solid")

                # Show the legend
                ax.legend(loc="lower right", numpoints=1)

                # Label it
                ax.text(0.05, 0.95, testing_data_labels_no_units[testing_data_key],
                        horizontalalignment='left', verticalalignment='top', transform=ax.transAxes,
                        fontsize=text_size)

                # Save and show it
                output_filename = join(args.workdir, args.output_file_name_head + "_" +
                                       testing_data_key + "_" + measurement_key + "_2D." + args.output_format)
                pyplot.savefig(output_filename, format=args.output_format, bbox_inches="tight", pad_inches=0.05)
                if not args.hide:
                    fig.show()
                else:
                    pyplot.close()

        if args.plot_fractional_limits:

            # Make plots of the fractional limits for each method

            # Set up the figure
            fig = pyplot.figure()

            pyplot.title(titles[testing_data_key])

            fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel("Method", fontsize=fontsize)
            ax.set_ylabel("Allowed $\\Delta$ on " + testing_data_labels_no_units[testing_data_key], fontsize=fontsize)
            ax.set_yscale("log", nonposy="clip")
            ax.set_ylim(1e-4, 1e0)

            xticks = list(range(len(args.methods)))
            ax.set_xticks(xticks)
            xticklabels = []
            for method in args.methods:
                if "_big" not in method:
                    xticklabels.append(method)
                else:
                    xticklabels.append("Big " + method.replace("_big", ""))
            ax.set_xticklabels(xticklabels)

            for measurement_key in ("m", "c"):
                for target_key in ("high", "base"):

                    limits = []
                    for method in args.methods:
                        limits.append(fractional_limits[method + "_" + measurement_key + "_" + target_key])

                    ax.scatter(xticks, limits, label=measurement_key + " " + target_labels[target_key],
                               marker=target_shapes[target_key], color=measurement_colors[measurement_key],
                               s=256)

            ax.legend(loc="upper right", scatterpoints=1)

            # Save and show it
            output_filename = join(args.workdir, args.output_file_name_head + "_" +
                                   testing_data_key + "_fractional_limits." + args.output_format)
            pyplot.savefig(output_filename, format=args.output_format, bbox_inches="tight", pad_inches=0.05)

    return
