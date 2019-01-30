#!/usr/bin/env python
#%%
"""
Reads xyz files in a folder, puts coordinate information, atomic numbers,
number of atoms, etc. in associated hdf5 files for further analysis

Example:
   ./xyz2hdf5.py -i .xyzs/ --out_folder o --out_prefix tester
will create the file ./o/tester-cartgeom.h5
"""
__author__ = 'Adam Duster'
__copyright__ = ''
__credits__ = ['Adam Duster']
__license__ = 'CC-BY-SA'
__version__ = '0.1'
__email__ = 'adam.duster@ucdenver.edu'
__status__ = 'Development'

import argparse
import glob
import os
import h5py as h5
import numpy as np
import indicator_tools
import time

#%%
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
        required=False)
    parser.add_argument(
        '--in_folder',
        help="Folder with coordinate files to add",
        required=True
    )
    parser.add_argument(
        '--out_folder',
        help="The folder for the h5py files",
        required=False,
        default="./"
    )
    parser.add_argument(
        '--out_prefix',
        help="Prefix to add for h5py files",
        required=False,
        default=''
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
    parser.add_argument(
        '-a',
        '--max_atoms',
        help="Max atoms for size of h5 array",
        type=int,
        default=151
    )
    return parser.parse_args(args)


#%%
def main():
#%%
    ## Vars
    test = False
    arg_vals = None
    if test:
        arg_vals = ['--in_folder', './xyz_a_tests/',
                '--out_folder', './xyz_a_out',
                '--out_prefix', 'small_test']
    if test:
        arg_vals = ['--in_folder', './behler_xyz/',
                '--out_folder', './behler_xyz',
                '--out_prefix', 'behler_3b']

    args = get_args(arg_vals)

    args.out_stem = args.out_folder + '/'
    if args.out_prefix:
        args.out_stem += args.out_prefix + '-'
    args.in_folder += '/'
    out_type = "xyz"

    coordinates = True
    eigenspectra = False

    # Get the file list
    globarg = args.in_folder + "*." + out_type
    flist = glob.glob(globarg)
    flist.sort()
    num_files = len(flist)

    # Make the h5py files
    max_atoms = 151
    of_path = args.out_stem + "cartgeom.h5"
    ofi = h5.File(of_path, "w")
    dcoords = ofi.create_dataset("cartesian_coords", shape=(num_files,
                                 max_atoms, 3), dtype=float)
    datmnum = ofi.create_dataset("atomic_number", shape=(num_files, max_atoms),
                                 dtype=np.int8)
    dnatoms = ofi.create_dataset("num_atoms", shape=(num_files,),
                                 dtype=np.int16)
    dnames = ofi.create_dataset("fpath", shape=(num_files,),
                                dtype=np.dtype("S25"))
    delsym = ofi.create_dataset("elsym", shape=(num_files, max_atoms),
                                dtype=np.dtype("S2"))

    chunk_size = 10000
    chunk = 0
    counted = 0

    tmp_coords = np.empty((chunk_size, max_atoms, 3), dtype=float)
    tmp_natoms = np.empty(chunk_size, dtype=np.int16)
    tmp_names = np.zeros(chunk_size, dtype="S25")
    tmp_elsym = np.zeros((chunk_size, max_atoms), dtype="S2")
    tmp_atmnum = np.empty((chunk_size, max_atoms), dtype=np.int8)
    begin_time = time.time()
    lfi = None
    while counted < num_files:
        loop_start_time = time.time()
        if counted + chunk_size > num_files:
            chunk_size = num_files - counted
        start = counted
        end = counted + chunk_size

        for i in range(chunk_size):
            tmp_natoms[i], tc, te, ta =  \
                 indicator_tools.read_xyz2(flist[start + i])
            tmp_atmnum[i, :tmp_natoms[i]] = ta
            tmp_coords[i, :tmp_natoms[i], :] = tc
            tmp_elsym[i, :tmp_natoms[i]] = te
            tmp_names[i] = flist[start + i]

        datmnum[start:end] = tmp_atmnum[:chunk_size]
        dcoords[start:end] = tmp_coords[:chunk_size]
        dnames[start:end] = tmp_names[:chunk_size]
        delsym[start:end] = tmp_elsym[:chunk_size]
        dnatoms[start:end] = tmp_natoms[:chunk_size]
        tmp_coords[:] = 0.0
        tmp_elsym[:] = ''
        ofi.flush()
        loop_time = time.time() - loop_start_time
        print("loop time: {0}   chunk: {1}".format(loop_time, end), file=lfi)
        counted = end
    print("Total time: {0}    Total counted:{1}".format(time.time() - begin_time, end), file=lfi)








if __name__ == "__main__":
    main()
