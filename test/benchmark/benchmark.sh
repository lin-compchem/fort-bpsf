#!/usr/bin/env bash
#
# Run the executable!
#
ifi=../../test_files/bench_md_1000-engrad.h5
ofi=./bench_md_1000-symfuncs.h5
bin=../../bin/gen_symfuncs_parallel

cwd=`pwd`
echo $cwd
rm $bin
cd ../../src
make all
cd $cwd
rm $ofi
$bin $ifi $ofi

