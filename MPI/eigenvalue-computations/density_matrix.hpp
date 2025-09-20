#ifndef RANDOM_RHO
#define RANDOM_RHO

#include <vector>
#include <iostream>
#include <complex>
#include <random>
#include <cblas.h>
#include <omp.h>

//alias
using Complex = std::complex<double>;
//using Vector = std::vector<Complex>;


Complex genMatrixElements(std::mt19937 &gen, std::normal_distribution<double> & dist);
std::vector<Complex> constructRandomRho(int num_states);
std::vector<Complex> constructRhoFromStatevector(std::vector<Complex>& psi);



#endif