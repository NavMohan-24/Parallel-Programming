#!/bin/bash

# Exit if any command fails
set -e

# create output directory if it doesn't exist
mkdir -p ./bin

# compile Hamiltonian
g++-15 -c -fopenmp ./hamiltonian/hamiltonian.cpp -o ./bin/hamiltonian.o

# compile random density matrix
g++-15 -c -fopenmp ./density_matrix/density_matrix.cpp -o ./bin/density_matrix.o

# compile jump operator
g++-15 -c -fopenmp ./jump_operator/jump_operator.cpp -o ./bin/jump_operator.o 

# compile the Linbladian constructor
g++-15 -c -fopenmp ./linbladian_constructor/linbladian_constructor.cpp -o ./bin/linbladian.o 

# compile diagonalizer
g++-15 -c ./arnoldi_diagonalizer/linbladian_diagonalizer.cpp -o ./bin/diagonalizer.o 

# run function to find liouvillian spectrum of TFIM
g++-15 -fopenmp tfim_linbladian_spectrum.cpp ./bin/diagonalizer.o ./bin/linbladian.o ./bin/density_matrix.o ./bin/jump_operator.o ./bin/hamiltonian.o -o res_spectrum.out -lblas -llapacke -std=c++20