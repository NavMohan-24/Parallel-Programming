#ifndef LINBLADIAN_HPP
#define LINBLADIAN_HPP

#include "jump_operator/jump_operator.hpp"
#include <cblas.h>
#include <omp.h>
#include <vector>
#include <complex>
#include <stdexcept>

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

        std::vector<Complex>applyLinbladian(const std::vector<Complex>& rho);

        std::vector<double> constructHessenbergMatrix(const std::vector<Complex>& rho, int k, double tol = 1e-12);

    protected:
        
        std::vector<Complex> applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);
        std::vector<Complex> applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);
        std::vector<Complex> constructDissipator(const std::vector<Complex>& mat,int N, int num_states, double rate, DecayType decay, Scope scope);
        double computeInnerProduct(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M, double tol = 1e-12);
        void normalizeMatrix(std::vector<Complex>& mat, int M);
    
};

#endif