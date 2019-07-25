#!/usr/bin/env python
""" 
"""
__author__ = 'Adam Duster'
__copyright__ = ''
__credits__ = ['Adam Duster']
__license__ = 'CC-BY-SA'
__version__ = '0.1'
__email__ = 'adam.duster@ucdenver.edu'
__status__ = 'Development'

import argparse
import h5py as h5
import numpy as np
import os
import shutil
import subprocess

def get_args(args=None):
    """ This is written as a default funtion to put at beginning of all Python
    scripts which require command line arguments. This uses the argparse module
    which must be declared in the main program to ensure that the object is able
    to be used by the caller
    --Adam Duster 21 June 2017
    """
    parser = argparse.ArgumentParser(description='see header of python script')
    parser.add_argument(
        '-i',
        '--input',
        help='Input xyz file',
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        help='Output xyz file',
        required=True)
    parser.add_argument(
        '-s',
        '--step',
        help='The change in theta between the two files',
        required=False,
        default=np.pi/100000,
        type=float
    )
    parser.add_argument(
        '-t',
        '--theta',
        help='The initial theta for the first geometry',
        required=False,
        default=np.pi/2,
        type=float
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help='Controls the level of output, use multipe v for more output',
        required=False,
        action='count',
        default=0)
    parser.add_argument(
        '-d',
        '--debug',
        help='Enter debug mode',
        required=False,
        action='store_true',
        default=False)
    return parser.parse_args(args)

def write_coords(c, fpath):
    s = "\n{0:1s}      {1:12.10f}      {2:12.10f}       {3:12.10f}"
    with open(fpath, 'w') as ofi:
        ofi.write('3\n')
        ofi.write(s.format('O',0.,0.,0.))
        ofi.write(s.format('H',1.,0.,0.))
        ofi.write(s.format('H',c[0],c[1],0.0))
    return
## Vars
arg_vals = None
args = get_args(arg_vals)
print("This script makes a perturbed coordinate file")
print("It adds %e to the angle of betweeen the 3-1-2 angle" % args.step)
print("The chosen theta is %e and the chosen step is %e" % (args.theta, args.step))
#
# Calculate the original coordinates
#
c = np.asarray([0.,0.])
c[0] = np.cos(args.theta)
c[1] = np.sin(args.theta)
write_coords(c, args.input)
#
# We use polar coordinates here to output the new x,y,z
#
new_theta = args.theta + args.step
c[0] = np.cos(new_theta)
c[1] = np.sin(new_theta)
write_coords(c, args.output)
