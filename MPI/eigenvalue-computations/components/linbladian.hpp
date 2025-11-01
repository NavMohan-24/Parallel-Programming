#ifndef LINBLADIAN_HPP
#define LINBLADIAN_HPP

#include "jump_operator/jump_operator.hpp"
#include <cblas.h>
#include <omp.h>
#include <vector>
#include <complex>
#include <stdexcept>
#include <lapacke.h>
#include <algorithm>

using Complex = std::complex<double>;

struct EigenResult {
    std::vector<Complex> eigenvalues;
    std::vector<Complex> eigenvectors;
};


class LinbladianDiagonalizer
{

    private:
        // Hamiltonian of the system
        std::vector<Complex> hamiltonian;
        // number of quantum systems.
        int N;
        // dimension of the quantum systems 
        int d;
        // hilbert space dimension
        int num_states;
        // decay rate
        double rate = 1.0;
        // decay type
        DecayType decay = DecayType::Damping;
        // scope of decay - local implies spin wise decay
        Scope scope = Scope::Local;
        // container to store the krylov vector in flattened manner
        std::vector<Complex>  kyrlov_basis_flatten;
        
        
    public:
         /**
         * @brief Constructs the diagonalizer from a Hamiltonian and specifications of dissipations.
         * 
         * @param H System Hamiltonian matrix.
         * @param N Number of quantum systems.
         * @param d Dimension of quantum systme.
         * @param rate Decay rate.
         * @param scope Scope of the decay. Whether decay happens locally or globally.
         * 
         */
        LinbladianDiagonalizer(const std::vector<Complex>& hamiltonian_,
                        int N_, int d_, double rate_ = 1.0, 
                        DecayType decay_ = DecayType::Damping,  
                        Scope scope_ = Scope::Local);
        /**
         * @brief Applies Linbladian on density matrix.
         * 
         * Computes:
         * \f[
         * \mathcal{L}(\rho) = -i[H, \rho] + \sum_{k}
         * \left( L_k \rho L_k^{\dagger} - \frac{1}{2} \{ L_k^{\dagger} L_k, \rho \} \right)
         * \f]
         * 
         * @param rho Input density matrix.
         * 
         *  */    
        std::vector<Complex> applyLinbladian(std::vector<Complex>& rho);
        EigenResult diagonalize(const std::vector<Complex>& initial_rho, int k, bool descending = false, double tol = 1e-12);
        std::vector<double> constructHessenbergMatrix(const std::vector<Complex>& initial_rho, int k, double tol = 1e-12); // need to make it protected at end


    protected:
        /**
         * @brief Constructs Hessenberg Matrix from the Linbladian through Arnoldi iteration.
         * 
         * Constructs Krylov Subspace from initial density matrix iteratively. Each iteration
         * calls the `applyLinbladian` method to generate Krylov vectors.
         * 
         * @param rho Input density matrix.
         * 
         * @return H Vectorized Hessenberg matrix.
         */
        std::vector<double> constructHessenbergMatrix(const std::vector<Complex>& initial_rho, int k, double tol = 1e-12); // need to make it protected at end
        std::vector<Complex> applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);
        std::vector<Complex> applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);
        std::vector<Complex> constructDissipator(const std::vector<Complex>& mat,int N, int num_states, double rate, DecayType decay, Scope scope);
        double computeInnerProduct(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M, double tol = 1e-12);
        void normalizeMatrix(std::vector<Complex>& mat, int M); // equivalent to setting norm = None
        void sortEigenPairs(std::vector<Complex>& eigenvalues,std::vector<Complex>& eigenvectors,int row, int col, bool descending = false);
        
    
};

#endif