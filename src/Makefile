COMPILER=gfortran
OBJS=bp_symfuncs.f08 h5_file_info.f08 main.f08
default: all
debug: $(OBJS)
	$(COMPILER)  $(OBJS) -fopenmp -fcheck=all -g -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran -o ../bin/gen_symfuncs_debug -O0 
all: $(OBJS)
	$(COMPILER)  $(OBJS) -fopenmp -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran  -o ../bin/gen_symfuncs_parallel -O3
serial: $(OBJS)
	$(COMPILER)  $(OBJS) -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran -o ../bin/gen_symfuncs_serial -O3
profile: $(OBJS)
	$(COMPILER)  $(OBJS) -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran -g -O5 -pg -o ../bin/gen_symfuncs_profile
clean:
	rm -rf *.mod
	rm -rf ../bin/gen_symfuncs*
