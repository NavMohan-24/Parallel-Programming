#ifndef TFIM_HPP
#define TFIM_HPP

#include <vector>
#include <iostream>

//type alias complex data type
using Complex = std::complex<double>;

// forward declaration - telling compiler that the function exist somewhere.
std::vector<Complex> tfimHamiltonian(int num_states, int N, double J, double h);

#endif