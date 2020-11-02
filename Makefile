FC=gfortran
FLAGS=-O3

all:
	$(FC) $(FLAGS) trace.f90 -o trace

