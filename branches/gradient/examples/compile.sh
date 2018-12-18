#!/usr/bin/env bash
gfortran $1 -I /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran -o $2
