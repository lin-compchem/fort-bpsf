FC=gfortran
SRCS=../../src/bp_symfuncs.f08 caller.f08
default: debug
OBJS= bp_symfuncs.o caller.o 
FFLAGS= -cpp -Wunused -fopenmp -fcheck=all -g -lhdf5_fortran -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include 

LFLAGS=-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
MFLAGS=-fopenmp -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib

debug: caller

caller: $(OBJS)
	$(FC) $(LFLAGS) $(OBJS)  -o caller $(MFLAGS)

bp_symfuncs.o: ../../src/bp_symfuncs.f08
	$(FC) $(FFLAGS) -c ../../src/bp_symfuncs.f08

caller.o: caller.f08
	$(FC) $(FFLAGS) -c caller.f08

clean:
	rm -rf *.mod
	rm *.o
