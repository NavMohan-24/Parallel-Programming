#ifndef HAM_HPP
#define HAM_HPP

#include <vector>
#include <iostream>
#include <complex>
#include <boost/dynamic_bitset.hpp>


using Complex = std::complex<double>;

std::vector<Complex> constructHamiltonian(int N, int num_states, double J, double h);


#endif