#!/usr/bin/env bash
#
# Run the executable!
#
ifi=bp.inp
ofi=./output.h5
rfi=../../test_files/reference.h5
bin=../../bin/gen_symfuncs_parallel
cwd=`pwd`
cd ../../src
make all 
cd $cwd 
$bin $ifi $ofi 
