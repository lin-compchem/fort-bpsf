FC=gfortran
CC=g++
SRCS=../../src/bp_symfuncs.f08 passer.f08
default: debug
OBJS= bp_symfuncs.o reciever.o passer.o  
FFLAGS= -cpp -Wunused -fopenmp -fcheck=all -g -lhdf5_fortran -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -lrapidjson -lstdc++

LFLAGS=-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
MFLAGS=-fopenmp -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lstdc++

CFLAGS= -Wunused -cpp -lrapidjson -g -fbacktrace

debug: passer

passer: $(OBJS)
	$(FC) $(LFLAGS) $(OBJS)  -o passer $(MFLAGS)

bp_symfuncs.o: ../../src/bp_symfuncs.f08
	$(FC) $(FFLAGS) -c ../../src/bp_symfuncs.f08

passer.o: passer.f08
	$(FC) $(FFLAGS) -c passer.f08

reciever.o: reciever.cpp
	$(CC) $(CFLAGS) -o reciever.o -c reciever.cpp 
clean:
	rm -rf *.mod
	rm *.o
