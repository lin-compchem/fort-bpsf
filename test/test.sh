#!/usr/bin/env bash
#
# Run the executable!
#
echo $1
if [ -n "$1" ]; then
    if [ $1 == "debug" ]; then
        t="debug"
        export OMP_NUM_THREADS=1
    elif [ $1 == "serial" ]; then
        t="serial"
        export OMP_NUM_THREADS=1
    elif [ $1 == "regular" ]; then
        t="regular"
    else
	echo "Could not find test"
	exit 1
    fi
else
    t="regular"
fi

tdir="1_${t}"

echo "RUNNING TEST: $tdir"

ifi=../../test_files/2_water_clusters.h5
ofi=./2_waters-symfuncs.h5
rfi=../../test_files/reference.h5
bin=../../bin/gen_symfuncs_serial

cd $tdir
cwd=`pwd`
echo $cwd
rm $bin
cd ../../src
make serial
cd $cwd
rm $ofi
$bin $ifi $ofi
../../scripts/compare_basis_functions.py -i $ofi -r $rfi

