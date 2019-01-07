COMPILER=gfortran
OBJS=bp_symfuncs.f08 main.f08
default: debug
debug: $(OBJS)
	$(COMPILER) -o gen_symfuncs_debug $(OBJS) -I -fopenmp -fcheck=all -g -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran -lblas -o gen_symfuncs_debug -O0
