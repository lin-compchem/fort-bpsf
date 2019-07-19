CC=g++
CFLAGS=
FC=gfortran
FFLAGS=
SRC=tfserving_utils.cpp tfs_messages.cpp
OBJ=tfserving_utils.o tfs_messages.o tfserving_utils_capi.o curler.o 
LIBS= -lcurl -lgfortran -lstdc++
default: demo1
fdemo1: fdemo1.o tfserving_utils_mod.o $(OBJ)
	$(FC) $(FFLAGS) -o fdemo1.x $< $(OBJ) $(LIBS)
demo1: demo1.o $(OBJ)
	$(CC) $(CFLAGS) -o demo1x $< $(OBJ) $(LIBS)
.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
	rm demo1.x
	rm fdemo1.x
tfserving_utils.o: tfserving_utils.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
tfserving_utils_capi.o: tfserving_utils_capi.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
tfs_messages.o: tfs_messages.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
curler.o: curler.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
tfserving_utils_mod.o: tfserving_utils_mod.f03 tfserving_utils_capi.o
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)
demo1.o: demo1.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) 
fdemo1.o: fdemo1.f03 tfserving_utils_mod.o
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)
