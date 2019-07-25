#!/usr/bin/env bash
#
# Run the executable!
#
set -x
export OMP_NUM_THREADS=10
v='-vv'
v=''
ifi=bp.inp
ofi=./output.h5
rfi=../../test_files/python_reference.h5
rfi=../../test_files/reference.h5
rfi=../../test_files/test_3b_1-cartgeom.h5
bin=../../bin/gen_symfuncs_debug
cwd=`pwd`
rm $ofi
cd ../../src
make clean
make debug
cd $cwd 
$bin $ifi $ofi 
../../scripts/compare_basis_functions.py -i $ofi -r $rfi $v
echo "For debugging, type:"
echo "gdb --args $bin $ifi $ofi"
