#!/usr/bin/env bash
cwd=`pwd`
cd ../../src
make debug
mod_file=`realpath bp_symfuncs.mod`
cd $cwd

FLAGS="-fopenmp -fcheck=all -g -fbacktrace -O0 -fcheck=all -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran"

echo "Compiling"
#gfortran -c ../../src/bp_symfuncs.f08 $FLAGS -o bp_symfuncs.o
gfortran -c caller.f90 $FLAGS -o caller.o
gfortran -o caller $FLAGS caller.o bp_symfuncs.o


echo "Linking"

#gfortran -o caller -fopenmp -fcheck=all -g -fbacktrace -O0 -fcheck=all -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran caller.o $mod_file 
