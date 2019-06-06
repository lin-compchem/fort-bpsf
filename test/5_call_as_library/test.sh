#!/usr/bin/env bash
ofi=./output.h5
rfi=../../test_files/test_3b_1-fortsym.h5
make clean
make
./caller
../../scripts/compare_basis_functions.py -i $ofi -r $rfi
echo "For debugging, type:"
echo "gdb ./caller"
echo "gdb --args $bin $ifi $ofi"

