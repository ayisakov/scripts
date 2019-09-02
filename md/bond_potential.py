#!/usr/bin/env python3

# Usage:
# data file

import sys
import argparse
import math
from enum import Enum

def splitargs(argstr, num_args, style):
    args = argstr.split(',')
    if len(args) != num_args:
        raise Exception("Incorrect number of parameters for the {0:s} bond style".format(style))
    return args

def harmonic(argstr, distance):
    args = splitargs(argstr, 2, "harmonic")
    K = float(args[0])
    r_0 = float(args[1])
    return K * (distance - r_0)**2.

def fene(argstr, distance):
    args = splitargs(argstr, 4, "FENE")
    K = float(args[0])
    r_0 = float(args[1])
    epsilon = float(args[2])
    sigma = float(args[3])
    one = -0.5 * K * r_0**2. * math.log(1 - (distance / r_0)**2.)
    sigr = sigma / distance
    two = 4. * epsilon * (sigr**12. - sigr**6.)
    return one + two + epsilon

styles_supported = {"harmonic": harmonic, "fene": fene}

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--style", help="where STYLE is a string of the form "
                    "<style>:[style_arg_1,style_arg2,...,style_arg_n]\n"
                    "which contains the bond style, followed by style-specific arguments;"
                    " supported styles:\n"
                    "harmonic:<constant,equilibrium_distance>\n"
                    "fene:<constant,equilibrium_distance,epsilon,sigma>", type=str, required=True)
parser.add_argument("-d", "--distance", help="the separation between the bond atoms' centers of mass",
                    type=float, required=True)
args = parser.parse_args()
style_parms = args.style.split(':')
if len(style_parms) != 2:
    raise Exception("Malformed bond style string '{0:s}'".format(style_parms))
style = style_parms[0]
if style not in styles_supported.keys():
    raise Exception("Unsupported bond style {0:s}".format(style))
constants = style_parms[1]
distance = args.distance
print(styles_supported[style](constants, distance))
