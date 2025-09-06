#ifndef ISING_HPP
#define ISING_HPP

#include <vector>
#include <iostream>

//type alias complex data type
using Complex = std::complex<double>;

// forward declaration - telling compiler that the function exist somewhere.
std::vector<Complex> tfimHamiltonian(int N, double J, double h);

#endif