#!/usr/bin/env python
"""
Compare the arrays for the radial and angular basis functions
in two h5py files.

Example:
./read_hdf5.py -i myfile.h5 -r refernce_file.h5
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
        help='Input file name',
        required=True)
    parser.add_argument(
        '-r',
        '--reference',
        help='Reference file name',
        required=True)
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


def main():
    arg_vals = None
    args = get_args(arg_vals)
    keys2check = ['h_radial_sym_funcs',
                  'h_angular_sym_funcs',
                  'o_radial_sym_funcs',
                  'o_angular_sym_funcs']
    ifi = h5.File(args.input, 'r')
    rfi = h5.File(args.reference, 'r')
    for key in keys2check:
        i_shape = ifi[key].shape
        r_shape = rfi[key].shape
        if (i_shape != r_shape) and (key != 'o_radial_sym_funcs'):
            print("Error for key %s" % key)
            print("Input shape {0} does not equal reference shape {1}".format(i_shape, r_shape))
            continue
        big = 0
        small = 0
        equal = 0
        if key == 'o_radial_sym_funcs':
            for j in range(i_shape[0]):
                for k in range(0,24):
                    if np.isclose(ifi[key][j,k], rfi[key][j,k]):
                        equal += 1
                    elif ifi[key][j,k] > rfi[key][j,k]:
                        big += 1
                    elif ifi[key][j,k] < rfi[key][j,k]:
                        small += 1
            for j in range(i_shape[0]):
                for k in range(25,46):
                    if np.isclose(ifi[key][j,k+2], rfi[key][j,k]):
                        equal += 1
                    elif ifi[key][j,k+2] > rfi[key][j,k]:
                        big += 1
                    elif ifi[key][j,k+2] < rfi[key][j,k]:
                        small += 1
        else:
            for j in range(i_shape[0]):
                for k in range(i_shape[1]):
                    if np.isclose(ifi[key][j,k], rfi[key][j,k]):
                        equal += 1
                    elif ifi[key][j,k] > rfi[key][j,k]:
                        big += 1
                    elif ifi[key][j,k] < rfi[key][j,k]:
                        small += 1
        print( "SUMMARY FOR KEY %s" % key)
        print( "Entries equal: %d" % equal)
        print("Entries larger: %d" % big)
        print("Entries smaller: %d" % small)            
    ifi.close()
    rfi.close()
        


if __name__ == '__main__': main()

