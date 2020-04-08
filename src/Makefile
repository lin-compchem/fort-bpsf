FC=gfortran
SRC=str.f90 bp_symfuncs.f90 h5_file_info.f90 main.f90
LIBS= -I../obj -lcurl -lstdc++
MODDIR=../obj
ifeq ($(FC), gfortran)
    FFLAGS=-cpp -J$(MODDIR) -g 
    ifneq ($(MAKECMDGOALS), serial)
        FFLAGS+= -fopenmp
    endif
    CFLAGS=-g
    LIBS+=-L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lhdf5_fortran -lgfortran -I  /usr/lib/x86_64-linux-gnu/hdf5/serial/include
    ALL_FLAGS:= $(ALL_FLAGS) -cpp -J $(MODDIR) -g
endif
ifeq ($(FC), ifort)
    FFLAGS=-fpp -module $(MODDIR)
    ifneq ($(MAKECMDGOALS), serial)
        FFLAGS+= -qopenmp
    endif
	ifeq ($(MAKECMDGOALS), debug)
        FFLAGS=-fpp -module $(MODDIR) -g 

    endif
    CFLAGS= -g
    CFLAGS=
	LIBS+=-L/opt/hdf5-intel-1.10.5/lib -lhdf5_fortran -I/opt/hdf5-intel-1.10.5/include
endif
LIBDIR= ../lib
OBJDIR= ../obj
LIBMOD= $(addprefix $(LIBDIR) , string_functions.mod bp_symfuncs.mod)
LIBMOD= $(addprefix $(LIBDIR) ,  bp_symfuncs.mod)
OBJ= $(addprefix $(OBJDIR)/, str.o bp_symfuncs.o h5_file_info.o tfserving_utils_mod.o tfserving_utils_capi.o tfserving_utils.o tfs_messages.o curler.o fdemo1.o demo1.o )
default: all
debug: $(SRC)
	$(FC)  $(SRC) $(LIBS) $(FFLAGS) -Wunused -fcheck=all -g $(LIBS) -o ../bin/gen_symfuncs_debug -O0 
all: $(SRC)
	$(FC)  $(SRC) $(LIBS) $(FFLAGS) -o ../bin/gen_symfuncs_parallel -O5
serial: $(SRC)
	$(FC)  $(SRC) $(LIBS) $(FFLAGS) -o ../bin/gen_symfuncs_serial -O3
profile: $(SRC)
	$(FC)  $(SRC) $(LIBS) $(ALL_FLAGS) -g -O5 -pg -o ../bin/gen_symfuncs_profile
lib: $(LIBMOD)
	echo "library combiled in ../lib directory"
clean:
	@$(RM) *.mod
	@$(RM) $(OBJDIR)/*.a
	@$(RM) $(OBJDIR)/*.o
	@$(RM) $(OBJDIR)/*.mod
	@$(RM) ../bin/gen_symfuncs*
	@$(RM) $(LIBDIR)/*.a
	@$(RM) $(LIBDIR)/*.o
	@$(RM) $(LIBDIR)/*.mod
obj: $(OBJ)
	ar rvs $(OBJDIR)/bplib.a $(OBJ) 
string_functions.mod: str.f90
	$(FC) -c str.f90 -o str.o
bp_symfuncs.mod: bp_symfuncs.f90 string_functions.mod
	$(FC) -c bp_symfuncs.f90 $(ALL_FLAGS) -O3 -o $(LIBDIR)/bp_symfuncs.o
	mv bp_symfuncs.mod string_functions.mod str.o $(LIBDIR)
$(OBJDIR)/str.o: str.f90
	$(FC) -c $< $(FFLAGS) -o $@
$(OBJDIR)/bp_symfuncs.o: bp_symfuncs.f90 string_functions.mod
	$(FC) -c $< $(FFLAGS) -O3 -o $@ 
$(OBJDIR)/h5_file_info.o: h5_file_info.f90 $(OBJDIR)/bp_symfuncs.o
	$(FC) -c $< $(FFLAGS) -O3 -o $@ $(LIBS) 
$(OBJDIR)/tfserving_utils_mod.o: interface/tfserving_utils_mod.f90 $(OBJDIR)/tfserving_utils_capi.o
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)
$(OBJDIR)/tfserving_utils_capi.o: interface/tfserving_utils_capi.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
$(OBJDIR)/tfserving_utils.o: interface/tfserving_utils.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
$(OBJDIR)/tfs_messages.o: interface/tfs_messages.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
$(OBJDIR)/curler.o: interface/curler.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
$(OBJDIR)/demo1.o: interface/demo1.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
$(OBJDIR)/fdemo1.o: interface/fdemo1.f90 $(OBJDIR)/tfserving_utils_mod.o
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)
