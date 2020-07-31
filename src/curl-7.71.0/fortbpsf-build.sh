./configure --prefix=`readlink -f ../../`
make -j8
cd lib
make install
cd ../
cd include
make install
cd ../
