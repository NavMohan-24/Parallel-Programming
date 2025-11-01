#ifndef SPIN_LINBLADIAN_CONSTRUCTOR_HPP
#define SPIN_LINBLADIAN_CONSTRUCTOR_HPP

#include "jump_operator/jump_operator.hpp"
#include <cblas.h>
#include <omp.h>
#include <vector>
#include <complex>


using Complex = std::complex<double>;


/**
 * @class LinbladianConstructor
 * @brief Constructs and Diagonalizes Linbladian Superoperator via Arnoldi Iteration.
 * 
 * Construction of Linbladian superoperator is done implicitly with out full matrix reconstruction.
 * Arnoldi iteration is performed to iteratively construct a Krylov subspace, where the 
 * Linbladian is converted to equivalent Hessenberg matrix and diagonalized.
 * 
 */

class LinbladianConstructor
{

    private:
        // Hamiltonian of the system
        std::vector<Complex> hamiltonian;
        // number of quantum systems.
        int N;
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
         * @brief Constructs Linbladian for spin systems from a Hamiltonian and specifications of dissipations.
         * 
         * @param H System Hamiltonian matrix.
         * @param N Number of spins.
         * @param rate Decay rate.
         * @param decay Type of Decay. 
         *        Decay types supported::
         *          - Damping - Uses $\sigma_{\minus}$ as jump operator.
         *          - Pumping - Uses $\sigma_{\plus}$ as jump operator
         *          - Dephasing - Uses $\sigma_{z}$ as a jump operator
         * @param scope The scope of decay, whether decay happens acting Locally or Globally.
         * 
         */

        LinbladianConstructor(const std::vector<Complex>& hamiltonian_,
                        int N_, double rate_ = 1.0, 
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
         **/

        std::vector<Complex> applyLinbladian(std::vector<Complex>& rho);


    protected:
       
        /**
         * @brief Compute commutator between two matrix.
         * 
         * Takes two square matrix of same dimension as a input and computes 
         * the commutator between them.
         * 
         * @param matA Input vectorized square matrix.
         * @param matB Input vectorized square matrix.
         * @param M Dimension of the square matrix
         * 
         * @return Commutator between input matrices in vectorized form.
         */

        std::vector<Complex> applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);

         /**
         * @brief Compute anti-commutator between two matrix.
         * 
         * Takes two square matrix of same dimension as a input and computes 
         * the aniti-commutator between them.
         * 
         * @param matA Input vectorized square matrix.
         * @param matB Input vectorized square matrix.
         * @param M Dimension of the square matrix
         * 
         * @return Anti-Commutator between input matrices in vectorized form.
         */

        std::vector<Complex> applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);

        /**
         * @brief Construction of dissipator part of Linbladian.
         * 
         * Disipator part of Linbladian is defined as $\sum_{k} \gamma_{k} \left(L_{k}\rho L_{k}^{\dagger} - \frac{1}{2}\{L_{k}^{\dagger}L_{k}, \rho\}\right)$.
         * 
         * @param rho Input vectorized density matrix.
         * @param N Number of quantum systems/particles.
         * @param num_states Dimension of Hilbert Space.
         * @param rate Decay rate.
         * @param decay Type of Decay. It can be one of Damping, Pumping or Dephasing.
         * @param scope The scope of decay, whether it is acting Locally or Globally.
         * 
         * @return Anti-Commutator between input matrices in vectorized form.
         */

        std::vector<Complex> constructDissipator(const std::vector<Complex>& rho,int N, int num_states, double rate, DecayType decay, Scope scope);

        // helper function to enforce Hermiticity
        /**
         * @brief Helper Function to enforce Hermicity of a matrix.
         * 
         * Floating point errors in computation could affect the hermicity 
         * of the matrix. Perform (A+A^H)/2 on an input matrix A.
         * 
         * @param mat Input square matrix.
         * @param M  Dimension of the matrix.
         * @return void
         */

        void enforceHermiticity(std::vector<Complex>& mat, int M);
            
        };

#endif