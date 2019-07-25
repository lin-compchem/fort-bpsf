FC=gfortran
FFLAGS= -cpp -Wunused -fopenmp -fcheck=all -g -lhdf5_fortran -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -lrapidjson -lstdc++
CC=g++
CFLAGS=
SRCS=../../src/bp_symfuncs.f08 passer.f08
default: debug
OBJS= str.o bp_symfuncs.o tfserving_utils.o tfserving_utils_mod.o tfs_messages.o tfserving_utils_capi.o curler.o passer.o  reciever.o

LFLAGS=-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
MFLAGS=-fopenmp -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lstdc++ -lcurl

CFLAGS= -Wunused -cpp -lrapidjson -g 

debug: passer

passer: $(OBJS)
	$(FC) $(LFLAGS) $(OBJS)  -o passer $(MFLAGS)

str.o: ../../src/str.f08
	$(FC) $(FFLAGS) -c ../../src/str.f08

bp_symfuncs.o: ../../src/bp_symfuncs.f08
	$(FC) $(FFLAGS) -c ../../src/bp_symfuncs.f08

passer.o: passer.f08
	$(FC) $(FFLAGS) -c passer.f08

clean:
	rm -rf *.mod
	rm *.o
reciever.o: ./reciever.cpp
	    $(CC) $(CFLAGS) -c $< -o $@ $(LIBS)
tfserving_utils.o: ../../src/interface/tfserving_utils.cpp
	    $(CC) $(CFLAGS) -c $< -o $@ $(LIBS)
tfserving_utils_capi.o: ../../src/interface/tfserving_utils_capi.cpp
	    $(CC) $(CFLAGS) -c $< -o $@ $(LIBS)
tfs_messages.o: ../../src/interface/tfs_messages.cpp
	    $(CC) $(CFLAGS) -c $< -o $@ $(LIBS)
curler.o: ../../src/interface/curler.cpp
	    $(CC) $(CFLAGS) -c $< -o $@ $(LIBS)
tfserving_utils_mod.o: ../../src/interface/tfserving_utils_mod.f03 ../../src/interface/tfserving_utils_capi.o
	    $(FC) $(FFLAGS) -c $< -o $@ $(LIBS)
