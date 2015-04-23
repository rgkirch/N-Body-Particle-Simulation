cpp: mpi-gravity.cpp
	mpicxx -O3 mpi-gravity.cpp -o gravity.out -lsfml-graphics -lsfml-window -lsfml-system
	mpiexec -n 4 ./gravity.out
