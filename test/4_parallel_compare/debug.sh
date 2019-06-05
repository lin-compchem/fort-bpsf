#!/usr/bin/env bash
#
# Run the executable!
#
delta=1e-4
first=1oh.xyz
second=2oh.xyz
pre=h2o
h5=../3_ang_grad_small/${pre}-cartgeom.h5
o5=./${pre}-bpsf.h5
bin=../../bin/gen_symfuncs_debug
cwd=`pwd`
rm -rf *.h5 $bin

cd ../../src
make clean
make debug
cd $cwd 

echo "run the program"
$bin $h5 $o5
