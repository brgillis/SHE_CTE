""" @file measurement_extraction.py

    Created 26 Apr 2017

    Main function to plot bias measurements.
"""

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

import argparse
from os.path import join

from SHE_PPT import products
from SHE_PPT.file_io import read_xml_product

from astropy.table import Table
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.optimize import fsolve


products.shear_bias_measurements.init()

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

testing_data_labels = {"P": "PSF Size",
                       "S": "Sky Level",
                       "E": "P(e)"}

tag_template = "Ep0Pp0Sp0"

testing_variant_labels = ("m2", "m1", "p0", "p1", "p2")

measurement_key_templates = ("mDIM", "mDIM_err", "cDIM", "cDIM_err")
measurement_colors = {"m": "r", "c": "b"}

x_values = {"P": [0.8, 0.9, 1.0, 1.1, 1.2],
            "S": [8.0608794667689825, 9.0670343062718768, 10.073467500059127,
                  11.07982970368346,  12.086264012414167],
            "E": [0.18556590758327751, 0.21333834512347458, 0.2422044791810781,
                  0.27099491059570091, 0.30241731263684996],
            }

x_ranges = {"P": (0.75, 1.25),
            "S": (7.5, 12.5),
            "E": (0.170, 0.315)}

target_limit_factors = {"P": {"m": 64, "c": 32},
                        "S": {"m": 32, "c": 20},
                        "E": {"m": 24,  "c": 32}, }

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

fontsize = 12
text_size = 18


def plot_bias_measurements_from_args(args):
    """ @TODO main docstring
    """

    # Determine the qualified path to the root data folder
    if args.root_data_folder[0] == "/":
        root_data_folder = args.root_data_folder
    else:
        # Relative to workdir in this case
        root_data_folder = join(args.workdir, args.root_data_folder)

    # Open and and keep in memory all bias measurements
    all_bias_measurements = {}

    def read_bias_measurements(tag):
        if not tag in all_bias_measurements:
            all_bias_measurements[tag] = read_xml_product(join(root_data_folder, args.data_folder_head + tag +
                                                               "/she_measure_bias/shear_bias_measurements.xml"))

    # Do a loop of reading for each property
    for testing_variant in testing_variant_labels:

        e_tag = tag_template.replace("Ep0", "E" + testing_variant)
        read_bias_measurements(e_tag)

        p_tag = tag_template.replace("Pp0", "P" + testing_variant)
        read_bias_measurements(p_tag)

        s_tag = tag_template.replace("Sp0", "S" + testing_variant)
        read_bias_measurements(s_tag)

    # Plot the biases and errors for each measurement
    for testing_data_key in testing_data_labels:

        fractional_limits = {}

        for measurement_key_template in measurement_key_templates:

            # Plot regularly for dim = 0
            measurement_key = measurement_key_template.replace("DIM", "")
            measurement_key_1 = measurement_key_template.replace("DIM", "1")
            measurement_key_2 = measurement_key_template.replace("DIM", "2")

            # Set up the figure
            fig = pyplot.figure()
            fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel(testing_data_labels[testing_data_key], fontsize=fontsize)
            ax.set_ylabel("$" + measurement_key.replace("_err", r"_{\rm err}") + "$", fontsize=fontsize)

            all_methods_data = {}

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
                    lx.append(x_ranges[testing_data_key][i] + method_offsets[method] *
                              np.abs(x_values[testing_data_key][3] - x_values[testing_data_key][1]))

                    g1_bias_measurements, g2_bias_measurements = all_bias_measurements[tag].get_method_bias_measurements(
                        method)
                    ly1.append(getattr(g1_bias_measurements, measurement_key))
                    ly2.append(getattr(g2_bias_measurements, measurement_key))

                    # To calculate combined error, we also need non-error
                    if "_err" in measurement_key:
                        ly1_o.append(getattr(g1_bias_measurements, measurement_key.replace("_err", "")))
                        ly2_o.append(getattr(g2_bias_measurements, measurement_key.replace("_err", "")))
                    else:
                        # And for non-error values, we'll want to plot error bars too, so we need that
                        ly1_o.append(getattr(g1_bias_measurements, measurement_key + "_err"))
                        ly2_o.append(getattr(g2_bias_measurements, measurement_key + "_err"))

                x_vals = np.array(lx)
                y1_vals = np.array(ly1)
                y2_vals = np.array(ly2)
                y1_o_vals = np.array(ly1_o)
                y2_o_vals = np.array(ly2_o)

                # Determine combined y differently for error and non-error values
                if "_err" in measurement_key:
                    y_vals = (y1_o_vals * y1_vals + y2_o_vals * y2_vals) / np.sqrt(y1_o_vals**2 + y2_o_vals**2)
                    y_errs = None
                else:
                    y_vals = np.sqrt(y1_vals**2 + y2_vals**2)
                    y_errs = (y1_vals * y1_o_vals + y2_vals * y2_o_vals) / np.sqrt(y1_vals**2 + y2_vals**2)

                # Plot the values (and optionally error bars)
                if "_err" not in measurement_key:
                    ax.errorbar(x_vals, y_vals, y_errs, color=method_colors[method], linestyle='None')
                else:
                    y1_vals *= err_factor
                    y2_vals *= err_factor
                    y_vals *= err_factor
                ax.plot(x_vals, y_vals, color=method_colors[method], marker='o', linestyle='None')

                # Calculate and plot an interpolating spline
                y1_spline = Spline(x_vals, y1_vals)
                y2_spline = Spline(x_vals, y2_vals)

                def y_spline(x): return np.sqrt(y1_spline(x)**2 + y2_spline(x)**2)

                x_spline_vals = np.linspace(x_vals[0], x_vals[-1], 100)
                y_spline_vals = y_spline(x_spline_vals)

                if "_big" not in method:
                    label = method
                else:
                    label = "Big " + method.replace("_big", "")

                ax.plot(x_spline_vals, y_spline_vals, color=method_colors[method], marker='None',
                        label=label)

                # Save this data for the next plot
                method_data = {"x": x_vals,
                               "y": y_vals,
                               "y_err": y_errs,
                               "y1": y1_vals,
                               "y2": y2_vals,
                               "y1_o": y1_o_vals,
                               "y2_o": y2_o_vals,
                               }
                all_methods_data[method] = method_data

            # Plot the target line
            if "m" in measurement_key:
                target = m_target
            else:
                target = c_target

            xlim = x_ranges[testing_data_key]

            ax.plot(xlim, [target, target], label=None, color="k", linestyle="dashed")
            ax.plot(xlim, [-target, -target], label=None, color="k", linestyle="dashed")
            ax.plot(xlim, [20 * target, 20 * target], label=None, color="k", linestyle="dotted")
            ax.plot(xlim, [-20 * target, -20 * target], label=None, color="k", linestyle="dotted")
            ax.plot(xlim, [0, 0], label=None, color="k", linestyle="solid")

            # Set the limits and scale
            ax.set_xlim(xlim)
            # ax.set_ylim(y_range) # FIXME - uncomment once we know about what the range will be
            ax.set_yscale("log", nonposy="clip")

            # Show the legend
            ax.legend(loc="lower right", numpoints=1)

            # Save and show it
            output_filename = join(args.workdir, args.output_file_name_head + "_" +
                                   testing_data_key + "_" + measurement_key + "." + args.output_format)
            pyplot.savefig(output_filename, format=args.output_format, bbox_inches="tight", pad_inches=0.05)
            if not args.hide:
                fig.show()
            else:
                pyplot.close()

            # Now plot dim 1 v dim 2

            # Set up the figure
            fig = pyplot.figure()
            fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel(r"$\Delta " + measurement_key_1.replace("_err", r"_{\rm err}") + "$", fontsize=fontsize)
            ax.set_ylabel(r"$\Delta " + measurement_key_2.replace("_err", r"_{\rm err}") + "$", fontsize=fontsize)

            # Plot the target line
            if("m" in measurement_key):
                base_target = m_target
            else:
                base_target = c_target

            # Plot points for each method
            for method in args.methods:

                # Get the data we saved before
                method_data = all_methods_data[method]

                ax.plot(method_data["y1"], method_data["y2"],
                        color=method_colors[method], marker='o', linestyle='None')

                # Calculate and plot an interpolating spline
                y1_spline = Spline(method_data["x"], method_data["y1"])
                y2_spline = Spline(method_data["x"], method_data["y2"])

                def y_spline(x): return np.sqrt(y1_spline(x)**2 + y2_spline(x)**2)

                x_spline_vals = np.linspace(x_vals[0], x_vals[-1], 100)

                y1_spline_vals = y1_spline(x_spline_vals)
                y2_spline_vals = y2_spline(x_spline_vals)

                if "_big" not in method:
                    label = method
                else:
                    label = "Big " + method.replace("_big", "")

                ax.plot(y1_spline_vals, y2_spline_vals, color=method_colors[method], marker='None',
                        label=label)

                if "_err" not in measurement_key:

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

                    print(("Fraction limit on " + testing_data_labels[testing_data_key] + " for method " +
                           method + " for " + measurement_key + ": " +
                           str(fractional_limits[limit_label_base + "_base"]) + ",\t" +
                           str(fractional_limits[limit_label_base + "_high"])))

            theta_vals = np.linspace(0, 2 * np.pi, 360)

            target_factor = target_limit_factors[testing_data_key][measurement_key.replace("_err", "")]

            xlim = (-20 * target_factor * base_target, 20 * target_factor * base_target)
            ax.set_xlim(xlim)
            # ax.set_xscale('symlog',linthreshx=target_factor*base_target)
            ylim = (-20 * target_factor * base_target, 20 * target_factor * base_target)
            ax.set_ylim(ylim)
            # ax.set_yscale('symlog',linthreshy=target_factor*base_target)

            ax.plot(base_target * np.cos(theta_vals), base_target *
                    np.sin(theta_vals), label=None, color="k", linestyle="dashed",)
            ax.plot(20 * base_target * np.cos(theta_vals), 20 * base_target *
                    np.sin(theta_vals), label=None, color="k", linestyle="dotted")
            ax.plot(xlim, [0, 0], label=None, color="k", linestyle="solid")
            ax.plot([0, 0], ylim, label=None, color="k", linestyle="solid")

            # Show the legend
            ax.legend(loc="lower right", numpoints=1)

            # Label it
            ax.text(0.05, 0.95, testing_data_labels[testing_data_key],
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

        # Make plots of the fractional limits for each method

        # Set up the figure
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)

        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel("Method", fontsize=fontsize)
        ax.set_ylabel("Allowed $\\Delta$ on " + testing_data_labels[testing_data_key], fontsize=fontsize)
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
        if not args.hide:
            fig.show()
        else:
            pyplot.close()

    return
