#include <cmath>
#include <iostream>
#include <cblas.h>


#include <ezarpack/arpack_solver.hpp>
#include <ezarpack/storages/eigen.hpp>
#include <ezarpack/version.hpp>

using namespace ezarpack; // need to check it
using namespace Eigen;

// Size of the matrix
const int N = 10000;

// The number of low-lying eigenvalues we want to compute
const int N_ev = 10;

void applying_commutator(hx, hz, hzz, rho, N){

    




}


int main(){

    std::cout << "ezARPACK verison :" << EZARPACK_VERSION << std::endl;

    using solver_t = arpack_solver<ezarpack::Complex, eigen_storage>; // type aliasing 
    // arpack_solver is a template class that can accept real and complex matrices.
    
    solver_t solver(N); // declares a variable called solver of the type solver_t

    using params_t = solver_t::params_t; // type aliasing params_t attribute of solver

    params_t params(N_ev, paramst::Smallest, true);










    
    




}