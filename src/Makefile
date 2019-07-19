FC=gfortran
SRC=str.f08 bp_symfuncs.f08 h5_file_info.f08 main.f08
ALL_FLAGS= -cpp
LIBS= -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include -L /usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran
LIBDIR= ../lib
LIBMOD= $(addprefix $(LIBDIR) , string_functions.mod bp_symfuncs.mod)
LIBMOD= $(addprefix $(LIBDIR) ,  bp_symfuncs.mod)
LIBOBJ= $(addprefix $(LIBDIR) /, bp_symfuncs.o str.o)
default: all
#debug: $(SRC)
#	$(FC)  $(SRC) $(ALL_FLAGS) -Wunused -fopenmp -fcheck=all -g $(LIBS) -o ../bin/gen_symfuncs_debug -O0 
debug: $(SRC)
	$(FC)  $(SRC) $(ALL_FLAGS) -fopenmp -Wunused -fcheck=all -g $(LIBS) -o ../bin/gen_symfuncs_debug -O0 
all: $(SRC)
	$(FC)  $(SRC) $(LIBS) $(ALL_FLAGS) -fopenmp -o ../bin/gen_symfuncs_parallel -O5
serial: $(SRC)
	$(FC)  $(SRC) $(LIBS) $(ALL_FLAGS) -o ../bin/gen_symfuncs_serial -O3
profile: $(SRC)
	$(FC)  $(SRC) $(LIBS) $(ALL_FLAGS) -g -O5 -pg -o ../bin/gen_symfuncs_profile
lib: $(LIBMOD)
	echo "library combiled in ../lib directory"
clean:
	rm -rf *.mod
	rm -rf ../bin/gen_symfuncs*
	rm ../lib/*
string_functions.mod: str.f08
	$(FC) -c str.f08 -o str.o
bp_symfuncs.mod: bp_symfuncs.f08 string_functions.mod
	$(FC) -c bp_symfuncs.f08 $(ALL_FLAGS) -fopenmp -O3 -o $(LIBDIR)/bp_symfuncs.o
	mv bp_symfuncs.mod string_functions.mod str.o $(LIBDIR)
