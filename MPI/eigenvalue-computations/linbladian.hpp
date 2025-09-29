#ifndef LINBLADIAN_HPP
#define LINBLADIAN_HPP

#include "jump_operator.hpp"
#include<cblas.h>
#include<omp.h>
#include <vector>
#include <complex>

using Complex = std::complex<double>;

class LinbladianSolver
{

    private:
        std::vector<Complex> hamiltonian;
        int N;
        int num_states;
        double rate = 1.0;
        DecayType decay = DecayType::Damping;
        Scope scope = Scope::Local;

    public:

        LinbladianSolver(const std::vector<Complex>& hamiltonian_,
                        int N_, double rate_ = 1.0, 
                        DecayType decay_ = DecayType::Damping,  
                        Scope scope_ = Scope::Local);

        std::vector<Complex> solve(const std::vector<Complex>& rho);

    protected:
        
        std::vector<Complex> applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int num_states);
        std::vector<Complex> applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int num_states);
        std::vector<Complex> constructDissipator(const std::vector<Complex>& mat,int N,int num_states, double rate, DecayType decay, Scope scope);
    
};

#endif