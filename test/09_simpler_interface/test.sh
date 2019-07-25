#!/usr/bin/env bash
set -x
echo "Cleaning Test Directory"
make clean
echo "Make the program"
make
if [ $? -ne 0 ]; then
    exit
fi

echo "Start the tensorflow server and spam generic rest api calls"
source ./start_tfserving
./passer
if [ $? -ne 0]; then
    echo "DId the program crash? Error code .ne. 0"
fi
kill -9 $tfs_pid
echo "For debugging, type:"
echo "gdb ./passer"
echo "gdb --args $bin $ifi $ofi"

echo "ENERGY SHOULD BE: 11.9272804"
