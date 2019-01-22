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
        help='Input h5py file',
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        help='Output h5py file',
        required=True)
    parser.add_argument(
        '-s',
        '--step',
        help='The change in coordinates between the two files',
        required=False,
        default=1e-6,
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

## Vars
arg_vals = None
args = get_args(arg_vals)
print("This script makes a perturbed coordinate file")
print("It adds %e to each cartesian coordinate" % args.step)
print("The following output is the original and perturbed coordinate for the" \
        "first atom")
#
# Make the perturbed coordinate file
#
key = 'cartesian_coords'
shutil.copyfile(args.input, args.output)
with h5.File(args.output, 'r+') as ofi:
    data = ofi[key]
    data[0,0,0] = ofi[key][0,0,0] + args.step
    ofi.close()
with h5.File(args.output, 'r') as ofi, h5.File(args.input, 'r') as ifi:
    print(ifi[key][0,:2,:], ofi[key][0,:2,:])
