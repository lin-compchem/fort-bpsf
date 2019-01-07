#!/usr/bin/env python
"""
Reads an HDF5 file and outputs the name of all of the groups
and the shape of the array for the groups. The file can be 
specified by passing -i to the program.

Example:
./read_hdf5.py -i myfile.h5
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
    with h5.File(args.input, 'r') as ifi:
        keys = list(ifi.keys())
        for key in keys:
            key_shape = ifi[key].shape
            print(key, "shape=",key_shape)
        


if __name__ == '__main__': main()

