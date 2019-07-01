COMPILER=gfortran
OBJS=str.f08 bp_symfuncs.f08 h5_file_info.f08 main.f08
ALL_FLAGS= -cpp
LIBS= -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran
default: all
debug: $(OBJS)
	$(COMPILER)  $(OBJS) $(ALL_FLAGS) -Wunused -fopenmp -fcheck=all -g $(LIBS) -o ../bin/gen_symfuncs_debug -O0 
all: $(OBJS)
	$(COMPILER)  $(OBJS) $(LIBS) $(ALL_FLAGS) -fopenmp -o ../bin/gen_symfuncs_parallel -O5
serial: $(OBJS)
	$(COMPILER)  $(OBJS) $(LIBS) $(ALL_FLAGS) -o ../bin/gen_symfuncs_serial -O3
profile: $(OBJS)
	$(COMPILER)  $(OBJS) $(LIBS) $(ALL_FLAGS) -g -O5 -pg -o ../bin/gen_symfuncs_profile
clean:
	rm -rf *.mod
	rm -rf ../bin/gen_symfuncs*
