#! /bin/bash
mpicxx -std=c++11 -O3 mpi-gravity.cpp -lsfml-graphics -lsfml-window -lsfml-system && mpiexec -n 4 ./a.out
