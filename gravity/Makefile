mpi: mpi-gravity.cpp
	mpicxx -O2 mpi-gravity.cpp -o gravity.out -lsfml-graphics -lsfml-window -lsfml-system
	mpiexec -n 4 ./gravity.out
omp: omp-gravity.cpp
	g++ -o2 -fopenmp omp-gravity.cpp -o omp-gravity.out -lsfml-graphics -lsfml-window -lsfml-system

