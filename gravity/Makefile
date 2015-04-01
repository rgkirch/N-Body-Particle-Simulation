MPCC = mpicxx
CFLAGS = -O3

TARGETS = mpi-gravity

all:	$(TARGETS)

mpi-gravity: mpi-gravity.o common-gravity.o
	$(MPCC) -o $@ mpi-gravity.o mpi-common.o

common-gravity.o: common-gravity.cpp common-gravity.h
	$
