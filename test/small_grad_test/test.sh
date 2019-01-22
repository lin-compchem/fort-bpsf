#!/usr/bin/env bash
#
# Run the executable!
#
ifi="./h5/oh-cartgeom.h5"
ifirst="./step1_coords.h5"
ofirst="./step1_bf.h5"
inext="./step2_coords.h5"
onext="./step2_bf.h5"
bin=../../bin/gen_symfuncs_debug
cwd=`pwd`
export OMP_NUM_THREADS=1
rm -rf $ifirst $ofirst $inext $onext $bin

cd ../../src
make clean
make debug
cd $cwd 

printf "\n###############\nSTEP 1\n###############\n"
cp $ifi $ifirst
if [ ! -e $ifirst ]; then echo "could not find $ifirst file"; exit 1; fi;
python make_grad_coords.py -i $ifirst -o $inext -s 1e-6
if [ ! -e $inext ]; then echo "could not find the newly created coordinates at $inext file"; exit 2; fi;

printf "\n###############\nSTEP 2\n###############\n"
echo "Creating step 1 symfuncs"
echo "$ifirst $ofirst"
$bin $ifirst $ofirst
echo "$inext $onext"
echo "Creating step 2 symfuncs"
$bin $inext $onext 

printf "\n###############\nSTEP 3\n###############\n"

../../scripts/compare_basis_functions.py -i $ofirst -r $onext 
