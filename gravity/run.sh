#! /bin/bash
mpicxx -O3 mpi-gravity.cpp -lsfml-graphics -lsfml-window -lsfml-system && mpiexec -n 9 ./a.out -n 50
