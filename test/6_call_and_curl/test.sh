#!/usr/bin/env bash
set -x
ofi=./output.h5
rfi=../../test_files/test_3b_1-fortsym.h5
rm output.h5
make clean
make
if [ $? -ne 0 ]; then
    exit
fi
./passer
if [ $? -ne 0]; then
    exit
fi
#../../scripts/compare_basis_functions.py -i $ofi -r $rfi
echo "For debugging, type:"
echo "gdb ./passer"
echo "gdb --args $bin $ifi $ofi"

