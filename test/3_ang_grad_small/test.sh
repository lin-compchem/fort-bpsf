#!/usr/bin/env bash
#
# Run the executable!
#
delta=1e-4
first=1oh.xyz
second=2oh.xyz
pre=h2o
h5=${pre}-cartgeom.h5
o5=${pre}-bpsf.h5
bin=../../bin/gen_symfuncs_debug
cwd=`pwd`
export OMP_NUM_THREADS=1
rm -rf *.h5 $bin

cd ../../src
make clean
make debug
cd $cwd 

printf "\n###############\nSTEP 1\n###############\n"
echo "make the coordinate files"
#python make_grad_coords.py -i $first -o $second
cp ./ref_xyzs/${delta}/*.xyz ./ref_xyzs/
xyz_nn_analysis.py --in_folder ./ref_xyzs --out_folder ./ --out_prefix ${pre}

printf "\n###############\nSTEP 2\n###############\n"
echo "run the program"
$bin $h5 $o5
