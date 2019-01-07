#!/usr/bin/env bash
#
# Run the executable!
#
ifi=../../test_files/2_water_clusters.h5
ofi=./2_waters-symfuncs.h5
rfi=../../test_files/reference.h5
bin=../../bin/gen_symfuncs_parallel
cwd=`pwd`
cd ../../src
make all 
cd $cwd 
$bin $ifi $ofi 
../../test_files/compare_basis_functions.py -i $ofi -r $rfi 
