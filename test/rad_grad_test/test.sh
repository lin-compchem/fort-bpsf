#!/usr/bin/env bash
#
# Run the executable!
#
first=1oh.xyz
second=2oh.xyz
h5=oh-cartgeom.h5
o5=oh-bpsf.h5
bin=../../bin/gen_symfuncs_debug
cwd=`pwd`
export OMP_NUM_THREADS=1
rm -rf *.xyz *.h5 $bin

cd ../../src
make clean
make debug
cd $cwd 

printf "\n###############\nSTEP 1\n###############\n"
echo "make the coordinate files"
python make_grad_coords.py -i $first -o $second
xyz_nn_analysis.py --in_folder ./ --out_folder ./ --out_prefix oh

printf "\n###############\nSTEP 2\n###############\n"
echo "run the program"
$bin $h5 $o5

printf "\n###############\nSTEP 3\n###############\n"
echo "analyze the data"
../../scripts/compare_basis_functions.py -i $h5 -r $o5
