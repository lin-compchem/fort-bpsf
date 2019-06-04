COMPILER=gfortran
OBJS=curl.f90
default: debug
debug: $(OBJS)
	$(COMPILER) -o a.out $(OBJS) -fopenmp -fcheck=all -lcurl
clean:
	rm *.mod
	rm a.out
