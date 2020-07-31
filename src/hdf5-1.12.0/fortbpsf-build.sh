set -x
./configure --enable-fortran --prefix=`readlink -f ../../` --enable-tools=no --enable-tests=no
make -j16
make install
