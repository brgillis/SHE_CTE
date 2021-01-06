""" @file DrawSensTypes.py

    Created 5 January 2021

    Main function to draw a plot explaining different types of sensitivity testing.
"""

__updated__ = "2021-01-05"

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

from SHE_PPT.logging import getLogger
from SHE_PPT.utility import get_arguments_string

import SHE_CTE
import matplotlib.pyplot as pyplot


fig_title = "Types of Sensitivity Testing"
title_fontsize = 18

fig_xlabel = r"$X$"
fig_ylabel = r"$\hat{X}$"
axis_label_fontsize = 18

# Details for axis arrows
axis_arrow_headwidth = 0.035
axis_arrow_headlength = 0.035
axis_arrow_overhang = 0.
axis_linewidth = 1.

# Details for arrows in the plot
axis_arrow_headwidth = 0.025
axis_arrow_headlength = 0.025
axis_arrow_overhang = 0.3
plot_linewidth = 0.5


def defineSpecificProgramOptions():
    """
    @brief
        Defines options for this program.

    @return
        An ArgumentParser.
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_DrawSensTypes defineSpecificProgramOptions()')
    logger.debug('#')

    parser = argparse.ArgumentParser()

    # Output data
    parser.add_argument('--sens_types_plot', type=str, default="sens_types.pdf",
                        help='Desired name of the output sensitivity types plot.')

    # Arguments
    parser.add_argument('--hide', action="store_true",
                        help="Don't show any generated plots.")

    parser.add_argument('--profile', action='store_true',
                        help='Store profiling data for execution.')

    # Arguments needed by the pipeline runner
    parser.add_argument('--workdir', type=str, default=".")
    parser.add_argument('--logdir', type=str, default=".")

    logger.debug('# Exiting SHE_CTE_DrawSensTypes defineSpecificProgramOptions()')

    return parser


def mainMethod(args):
    """
    @brief
        The "main" method for this program, to measure bias.

    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = getLogger(__name__)

    logger.debug('#')
    logger.debug('# Entering SHE_CTE_DrawSensTypes mainMethod()')
    logger.debug('#')

    exec_cmd = get_arguments_string(args, cmd="E-Run SHE_CTE " + SHE_CTE.__version__ + " SHE_CTE_DrawSensTypes",
                                    store_true=["profile", "debug", "hide"])
    logger.info('Execution command for this step:')
    logger.info(exec_cmd)

    if args.profile:
        import cProfile
        cProfile.runctx("draw_sens_types_from_args(args)", {},
                        {"draw_sens_types_from_args": draw_sens_types_from_args,
                         "args": args}, filename="measure_bias.prof")
    else:
        draw_sens_types_from_args(args)

    logger.debug('# Exiting SHE_CTE_DrawSensTypes mainMethod()')

    return


def main():
    """
    @brief
        Alternate entry point for non-Elements execution.
    """

    parser = defineSpecificProgramOptions()

    args = parser.parse_args()

    mainMethod(args)

    return


if __name__ == "__main__":
    main()


def draw_sens_types_from_args(args):
    """ @TODO main docstring
    """

    # Set up the general layout of the plot
    fig = pyplot.figure()

    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.set_xticks([])
    ax.set_yticks([])

    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.075, right=0.75, top=0.9, left=0.075)

    # Remove the default axes
    for side in ["left", "right", "bottom", "top"]:
        ax.spines[side].set_visible(False)

    # Draw arrows as replacement axes
    ax.arrow(0., 0., 1.3, 0., facecolor='k', edgecolor='k', linewidth=axis_linewidth,
             head_width=axis_arrow_headwidth, head_length=axis_arrow_headlength, overhang=axis_arrow_overhang,
             length_includes_head=True, clip_on=False)

    ax.arrow(0., 0., 0., 1., facecolor='k', edgecolor='k', linewidth=axis_linewidth,
             head_width=axis_arrow_headwidth, head_length=axis_arrow_headlength, overhang=axis_arrow_overhang,
             length_includes_head=True, clip_on=False)

    # Draw plot and axis arrows
    pyplot.title(fig_title, fontsize=title_fontsize)
    ax.set_xlabel(fig_xlabel, fontsize=axis_label_fontsize)
    ax.set_ylabel(fig_ylabel, fontsize=axis_label_fontsize)

    # Draw an arrow for the ideal line
    ax.arrow(0.999, 0.999, 0.001, 0.001, facecolor="k", edgecolor="k",
             head_width=axis_arrow_headwidth, head_length=axis_arrow_headlength, overhang=axis_arrow_overhang,
             length_includes_head=True, clip_on=False)
    ax.plot([0., 1.], [0., 1.], linestyle=":", c="k")

    # Save and show the plot
    qualified_output_filename = join(args.workdir, args.sens_types_plot)
    pyplot.savefig(qualified_output_filename, format=args.sens_types_plot.split(".")[-1],
                   bbox_inches="tight", pad_inches=0.05)

    if not args.hide:
        pyplot.show()

    return
