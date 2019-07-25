#!/usr/bin/env bash
#
# Run the executable!
#
ifi=../../test_files/2_water_clusters.h5
ofi=./output.h5
rfi=../../test_files/reference.h5
bin=../../bin/gen_symfuncs_debug
cwd=`pwd`
rm $ofi
cd ../../src
make clean
make debug
cd $cwd 
$bin $ifi $ofi 
../../scripts/compare_basis_functions.py -i $ofi -r $rfi -v
echo "For debugging, type:"
echo "gdb --args $bin $ifi $ofi"
