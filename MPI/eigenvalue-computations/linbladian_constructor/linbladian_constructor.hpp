#ifndef SPIN_LINBLADIAN_CONSTRUCTOR_HPP
#define SPIN_LINBLADIAN_CONSTRUCTOR_HPP

#include "../jump_operator/jump_operator.hpp"
#include <cblas.h>
#include <omp.h>
#include <vector>
#include <complex>


using Complex = std::complex<double>;


/**
 * @brief Constructs and applies Lindbladian superoperators for open quantum systems.
 * 
 * This class provides an implicit construction of the Lindbladian superoperator without
 * full matrix reconstruction. The Lindbladian governs the dynamics of open quantum systems
 * including both unitary evolution (via the Hamiltonian) and non-unitary processes
 * (via dissipation operators).
 * 
 * The implementation supports various decay types (damping, pumping, dephasing) and can
 * apply dissipation either locally (per spin) or globally (collective decay).
 */

class LinbladianConstructor
{

    private:
        std::vector<Complex> hamiltonian;           // Hamiltonian of the quantum system
        int N;                                      // Number of quantum systems (e.g., spins)
        int num_states;                             // Hilbert space dimension (d^N)
        double rate;                                // Decay rate (gamma)
        DecayType decay;                            // Type of decay process
        Scope scope;                                // Scope of decay (local or global)
        std::vector<Complex> kyrlov_basis_flatten;  // Flattened storage for Krylov basis vectors
        
        
    public:

        /**
         * @brief Constructs a Lindbladian for spin systems with specified dissipation.
         * 
         * Initializes the Lindbladian superoperator with a given Hamiltonian and
         * dissipation parameters. The Lindbladian describes both coherent evolution
         * and incoherent decay processes in the open quantum system.
         * 
         * @param hamiltonian_ System Hamiltonian matrix (flattened, row-major order)
         * @param N_ Number of spins or quantum subsystems
         * @param rate_ Decay rate (gamma), default is 1.0
         * @param decay_ Type of decay process, default is Damping
         *        - Damping: Uses σ₋ as jump operator (energy decay)
         *        - Pumping: Uses σ₊ as jump operator (energy gain)
         *        - Dephasing: Uses σz as jump operator (phase decay)
         * @param scope_ Scope of decay, default is Local
         *        - Local: Decay acts independently on each spin
         *        - Global: Collective decay across all spins (!! TO BE IMPLEMENTED !!)
         */

        LinbladianConstructor(const std::vector<Complex>& hamiltonian_,
                        int N_, double rate_ = 1.0, 
                        DecayType decay_ = DecayType::Damping,  
                        Scope scope_ = Scope::Local);
        /**
         * @brief Applies the Lindbladian superoperator to a density matrix.
         * 
         * Computes the action of the Lindbladian on a density matrix:
         * \f[
         * \mathcal{L}(\rho) = -i[H, \rho] + \sum_{k}
         * \left( L_k \rho L_k^{\dagger} - \frac{1}{2} \{ L_k^{\dagger} L_k, \rho \} \right)
         * \f]
         * 
         * The first term represents unitary evolution under the Hamiltonian H,
         * while the sum represents dissipative processes characterized by jump
         * operators Lₖ.
         * 
         * @param rho Input density matrix (flattened, row-major order)
         * 
         * @return Resulting density matrix after applying the Lindbladian (flattened, row-major order)
         */

        std::vector<Complex> applyLinbladian(std::vector<Complex>& rho);


    protected:
       
        /**
         * @brief Computes the commutator between two matrices.
         * 
         * Calculates [A, B] = AB - BA for two square matrices of the same dimension.
         * Used to compute the unitary part of the Lindbladian evolution.
         * 
         * @param matA First input matrix (flattened, row-major order)
         * @param matB Second input matrix (flattened, row-major order)
         * @param M Dimension of the square matrices (M x M)
         * 
         * @return Commutator [A, B] (flattened, row-major order)
         */

        std::vector<Complex> applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);

        /**
         * @brief Computes the anti-commutator between two matrices.
         * 
         * Calculates {A, B} = AB + BA for two square matrices of the same dimension.
         * Used in the dissipative part of the Lindbladian to compute the anti-commutator
         * term {L†L, ρ}.
         * 
         * @param matA First input matrix (flattened, row-major order)
         * @param matB Second input matrix (flattened, row-major order)
         * @param M Dimension of the square matrices (M x M)
         * 
         * @return Anti-commutator {A, B} (flattened, row-major order)
         */

        std::vector<Complex> applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M);

        /**
         * @brief Constructs the dissipator part of the Lindbladian.
         * 
         * Computes the dissipative contribution to the Lindbladian:
         * \f[
         * \sum_{k} \gamma_{k} \left(L_{k}\rho L_{k}^{\dagger} - \frac{1}{2}\{L_{k}^{\dagger}L_{k}, \rho\}\right)
         * \f]
         * 
         * This term describes non-unitary evolution due to coupling with the environment.
         * The specific form depends on the decay type and scope parameters.
         * 
         * @param rho Input density matrix (flattened, row-major order)
         * @param N Number of quantum systems/spins
         * @param num_states Dimension of Hilbert space (d^N)
         * @param rate Decay rate (gamma)
         * @param decay Type of decay (Damping, Pumping, or Dephasing)
         * @param scope Scope of decay (Local or Global)
         * 
         * @return Dissipator contribution to dρ/dt (flattened, row-major order)
         */

        std::vector<Complex> constructDissipator(const std::vector<Complex>& rho,int N, int num_states, double rate, DecayType decay, Scope scope);

        /**
         * @brief Enforces Hermiticity of a matrix.
         * 
         * Corrects numerical errors that may violate Hermiticity by computing
         * (A + A†)/2 for an input matrix A. This ensures that density matrices
         * remain properly Hermitian throughout computations.
         * 
         * @param mat Input matrix to be symmetrized (flattened, row-major order), modified in-place
         * @param M Dimension of the square matrix (M x M)
         * 
         * @note This operation is performed in-place to minimize memory overhead
         */

        void enforceHermiticity(std::vector<Complex>& mat, int M);
            
        };

#endif