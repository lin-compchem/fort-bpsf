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
        default=False
    )
    parser.add_argument(
        '-t',
        '--test',
        help='The test you wish to conduct as an integer',
        default=0,
        type=int
    )
    return parser.parse_args(args)

def check2d(ifi, rfi, key, i_shape, verbose=False):
    equal = 0
    big = 0
    small = 0
    assert ifi[key].shape == rfi[key].shape
    for j in range(i_shape[0]):
        for k in range(i_shape[1]):
            fail = False
            if np.isclose(ifi[key][j,k], rfi[key][j,k]):
                equal += 1
                if verbose > 2:
                    print("PASS for key {} and indices {}   {}".format(key, j,
                                                                       k))
                    print("ref: {}        ifi: {}".format(rfi[key][j,k],
                        ifi[key][j,k]))
            elif ifi[key][j,k] > rfi[key][j,k]:
                big += 1
                fail = True
            elif ifi[key][j,k] < rfi[key][j,k]:
                small += 1
                fail = True
            if fail and verbose:
                print("FAIL For key {}, following inds failed: {}     {}".format(
                    key, j, k))
                if verbose > 1:
                    print("ref: {}        ifi: {}".format(rfi[key][j,k],
                        ifi[key][j,k]))
    return equal, big, small

def check4d(ifi, rfi, key, i_shape, verbose=False):
    equal = 0
    big = 0
    small = 0
    for j in range(i_shape[0]):
     for k in range(i_shape[1]):
         fail = False 
         if np.isclose(ifi[key][j,k], rfi[key][j,k]).all():
             equal += 1
             if verbose > 2:
                 print("PASS for key {} and indices {}   {}".format(key, j,
                                                                    k))
                 print("ref: {}        ifi: {}".format(rfi[key][j,k],
                     ifi[key][j,k]))
         elif (ifi[key][j,k] > rfi[key][j,k]).any():
             big += 1
             fail = True  
         elif (ifi[key][j,k] < rfi[key][j,k]).any():
             small += 1
             fail = True  
         if fail and verbose:
             print("FAIL For key {}, following inds failed: {}     {}".format(
                 key, j, k))
             if verbose > 1:
                 for rr , ii in zip(rfi[key][j,k], ifi[key][j,k]):
                     if (rr == np.zeros(3)).all() and (ii == np.zeros(3)).all():
                         continue
                     print("ref: {}        ifi: {}".format(rr, ii))
                 #print("ref: {}        ifi: {}".format(rfi[key][j,k],
                 #    ifi[key][j,k]))
    return equal, big, small
def main():
    arg_vals = None
    args = get_args(arg_vals)
    keys2check = [
            'h_radial_sym_funcs',
            'h_angular_sym_funcs',
            'o_radial_sym_funcs',
            'o_angular_sym_funcs',
            'h_rad_cartesian_gradient',
            'o_rad_cartesian_gradient',
            'h_ang_cartesian_gradient',
            'o_ang_cartesian_gradient'
            ]
    ifi = h5.File(args.input, 'r')
    rfi = h5.File(args.reference, 'r')
    #
    # Check to makesure we are comparing the same variables, etc.
    #
    ilist = list(ifi.keys())
    rlist = list(rfi.keys())
    num_keys = len(keys2check)
    for i in reversed(range(num_keys)):
        if keys2check[i] not in ilist or keys2check[i] not in rlist:
            print("Removing key from list: " + keys2check[i])
            del keys2check[i]
            num_keys -= 1
    #
    # The main loop
    #
    for key in keys2check:
        i_shape = ifi[key].shape
        r_shape = rfi[key].shape
        if args.verbose:
            print("**************************")
            print("Current key: %s" % key)
            print("Key shape: ", i_shape)
            print("Ref shape: ", r_shape)
        if (i_shape != r_shape) and (args.test == 0):
            print("Error for key %s" % key)
            print("Input shape {0} for key {2} does not equal reference shape {1}".format(i_shape, r_shape, key))
            continue
        #
        # The gradient array have different dimensions than the
        # other arrays
        #
        if 'gradient' in key:
            equal, big, small = check4d(ifi, rfi, key, i_shape, args.verbose)
        else:
            equal, big, small = check2d(ifi, rfi, key, i_shape, args.verbose)
        print("SUMMARY FOR KEY %s" % key)
        print("Entries equal: %d" % equal)
        print("Entries larger: %d" % big)
        print("Entries smaller: %d" % small)            
    ifi.close()
    rfi.close()
        


if __name__ == '__main__': main()

