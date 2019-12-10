#!/usr/bin/env python
#
"""
Reads xyz files in a folder, puts important spectra, coordinate information
in associated hdf5 files for further analysis
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
import structures
import sys

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
        '--in_folder',
        help="Folder with coordinate files to add",
        required=True
    )
    parser.add_argument(
        '-o',
        '--out_folder',
        help="The folder for the h5py files",
        required=False,
        default="./"
    )
    parser.add_argument(
        '-p',
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
        '-m',
        '--max_atoms',
        help='Maximum number of atoms in a geometry',
        default=151,
        type=int
    )
    parser.add_argument(
        '-c',
        '--code',
        help='Computer code used for the calculations. Can be mndo or mopac',
        required=True,
        default=''
    )
    return parser.parse_args(args)


def rLine(ifi):
    """
    By Adam Duster, Oct 12, 2017
    Read a line. Skip the line if there is a comment sign (#). If the line
    fails to read, raise an exception to be handeled by the calling program
    """
    while True:
        line = ifi.readline()
        if not line:
            raise IOError

        # skip comment
        if line.strip().startswith('#'):
            continue
        else:
            return line


#%%
def main():
#%%
    ## Vars
    test = 0
    arg_vals = None

    if test == 1:
        arg_vals = ['--in_folder', './fort15s',
                '--out_folder', './h5s',
                '--out_prefix', 'mndos',
                '--code', 'mndo']
    if test == 2:
        arg_vals = ['--in_folder', './mop_o',
                    '--out_folder', './h5s',
                    '--out_prefix', 'mops',
                    '--code', 'mopac']
    if test == 3:
        arg_vals = ['--in_folder', './engrad',
                    '--out_folder', './h5s',
                    '--out_prefix', 'pes',
                    '--code', 'orca']
    args = get_args(arg_vals)

    args.out_stem = args.out_folder + '/'
    if args.out_prefix:
        args.out_stem += args.out_prefix + '-'
    args.in_folder += '/'

    #Code specific stuff
    code = args.code.upper()
    if code == 'MNDO':
        out_type = "fort.15"
        of_suffix = "engrad.h5"
        FileClass = structures.Fort15File
    elif code == 'MOPAC':
        out_type = "out"
        FileClass = structures.MopacOutFile
        of_suffix = "engrad.h5"
    elif code == 'ORCA':
        out_type = "engrad"
        FileClass = structures.OrcaEngradFile
        of_suffix = "engrad.h5"
    else:
        sys.exit("Error, code not recognized")



    num2el = {1: 'H',
              8: 'O',
              6: 'C',
              7: 'N'}
    # Get the file list
    globarg = args.in_folder + "*." + out_type
    flist = glob.glob(globarg)
    flist.sort()
    num_files = len(flist)
    assert num_files > 1

    # Make the h5py files
    max_atoms = args.max_atoms

    of_path = args.out_stem + of_suffix
    ofi = h5.File(of_path, "w")
    dcoords = ofi.create_dataset("cartesian_coords", shape=(num_files,
                                 max_atoms, 3), dtype=float)
    datmnum = ofi.create_dataset("atomic_number", shape=(num_files, max_atoms),
                                 dtype=np.int8)
    dnatoms = ofi.create_dataset("num_atoms", shape=(num_files,),
                                 dtype=np.int16)
    dnames = ofi.create_dataset("fpath", shape=(num_files,),
                                dtype=np.dtype("S60"))
    delsym = ofi.create_dataset("elsym", shape=(num_files, max_atoms),
                                dtype=np.dtype("S2"))
    dener = ofi.create_dataset("potential_energies", shape=(num_files,),
                             dtype=float)
    dgrads = ofi.create_dataset("cartesian_grads", shape=(num_files,
                                max_atoms, 3), dtype=float)
    chunk_size = 10000
    chunk = 0
    counted = 0

    tmp_coords = np.full((chunk_size, max_atoms, 3), np.NaN, dtype=float)
    tmp_natoms = np.empty(chunk_size, dtype=np.int16)
    tmp_names = np.zeros(chunk_size, dtype="S60")
    tmp_elsym = np.zeros((chunk_size, max_atoms), dtype="S2")
    tmp_atmnum = np.zeros((chunk_size, max_atoms), dtype=np.int8)
    tmp_ener = np.empty((chunk_size), dtype=float)
    tmp_grad = np.full((chunk_size, max_atoms, 3), np.NaN, dtype=float)
    begin_time = time.time()
    lfi = None
    while counted < num_files:
        loop_start_time = time.time()
        if counted + chunk_size > num_files:
            chunk_size = num_files - counted
        start = counted
        end = counted + chunk_size

        for i in range(chunk_size):
            try:
                ifi = FileClass(flist[start + i])
                geom = ifi.parseAll()
                tmp_natoms[i] = geom.num_atoms
                tmp_atmnum[i, :tmp_natoms[i]] = geom.atomic_numbers
                tmp_coords[i, :tmp_natoms[i], :] = geom.atomic_coords
                tmp_elsym[i, :tmp_natoms[i]] = [num2el[n] for n in geom.atomic_numbers]
                tmp_names[i] = flist[start + i]
                tmp_ener[i] = geom.ener
                tmp_grad[i, :tmp_natoms[i], :] = geom.grad
                tmp_names[i] = flist[start + i]
            except:
                print("Error parsing file " + flist[start + i])
                sys.exit()

        datmnum[start:end] = tmp_atmnum[:chunk_size]
        dcoords[start:end] = tmp_coords[:chunk_size]
        dnames[start:end] = tmp_names[:chunk_size]
        delsym[start:end] = tmp_elsym[:chunk_size]
        dnatoms[start:end] = tmp_natoms[:chunk_size]
        dener[start:end] = tmp_ener[:chunk_size]
        dgrads[start:end] = tmp_grad[:chunk_size]
        tmp_coords[:] = np.NaN
        tmp_grad[:] = np.NaN
        tmp_names[:] = ''
        tmp_elsym[:] = ''
        tmp_atmnum[:] = 0

        ofi.flush()
        loop_time = time.time() - loop_start_time
        print("loop time: {0}   chunk: {1}".format(loop_time, end), file=lfi)
        counted = end
    print("Total time: {0}    Total counted:{1}".format(time.time() - begin_time,
                                                        end), file=lfi)


if __name__ == "__main__":
    main()
