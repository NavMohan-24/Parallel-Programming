#ifndef LINBLADIAN_DIAGONALIZER_HPP
#define LINBLADIAN_DIAGONALIZER_HPP

#include "./linbladian_constructor/linbladian_constructor.hpp"

#include <cblas.h>
#include <omp.h>
#include <vector>
#include <complex>
#include <stdexcept>
#include <lapacke.h>
#include <algorithm>


using Complex = std::complex<double>;

/**
 * @brief Structure to hold eigenvalue decomposition results.
 * 
 * Contains the eigenvalues and corresponding eigenvectors from
 * diagonalization of the Lindbladian operator.
 */

struct EigenResult {
    std::vector<Complex> eigenvalues;   // Eigenvalues of the Lindbladian
    std::vector<Complex> eigenvectors;  // Eigenvectors stored in row-major order
};

/**
 * @brief Diagonalizes Lindbladian operators using the Arnoldi iteration method.
 * 
 * This class implements the Arnoldi algorithm to compute eigenvalues and eigenvectors
 * of Lindbladian superoperators that govern the dynamics of open quantum systems.
 * The Arnoldi method constructs a Krylov subspace and projects the Lindbladian onto
 * this subspace to obtain a Hessenberg matrix, which is then diagonalized.
 * 
 * The implementation is optimized for large-scale quantum systems where full
 * diagonalization would be computationally prohibitive.
 */

class ArnoldiLinbladianDiagonalizer
{

    private:
        private:
        LinbladianConstructor linbladianConstructor;  // Constructor for the Lindbladian operator
        int N;                                         // Number of quantum systems in the composite system
        int d;                                         // Dimension of each individual quantum system
        int num_states;                                // Total Hilbert space dimension (d^N)
        std::vector<Complex> kyrlov_basis_flatten;     // Flattened storage for Krylov basis vectors
        
        
    public:

        /**
         * @brief Constructs the Arnoldi diagonalizer for a given Lindbladian.
         * 
         * Initializes the diagonalizer with the structure of the open quantum system,
         * including the Hamiltonian and dissipation operators through the provided
         * LinbladianConstructor.
         * 
         * @param linbladianConstructor_ Constructor object defining the Lindbladian operator
         * @param N_ Number of quantum systems in the composite system
         * @param d_ Dimension of each individual quantum system (e.g., 2 for qubits)
         */

        ArnoldiLinbladianDiagonalizer(LinbladianConstructor& linbladianConstructor_,int N_, int d_);

        /**
         * 
         * @brief Computes eigenvalues and eigenvectors of the Lindbladian operator.
         * 
         * Uses the Arnoldi iteration to build a Krylov subspace of dimension k starting
         * from an initial density matrix, then diagonalizes the resulting Hessenberg matrix
         * to obtain approximations to the Lindbladian's eigenspectrum.
         * 
         * @param initial_rho Initial density matrix (flattened) used to seed the Krylov subspace
         * @param k Dimension of the Krylov subspace to construct
         * @param descending If true, sort eigenvalues in descending order; otherwise ascending
         * @param tol Tolerance for numerical operations (orthogonalization, normalization)
         * 
         * @return EigenResult Structure containing eigenvalues and corresponding eigenvectors
         * 
         * @throws std::runtime_error LAPACK routines encounter errors.
         * 
         * @note The Arnoldi iteration may terminate early (at dimension < k) if the norm of a 
         * newly generated basis vector falls below the tolerance, indicating the Krylov subspace 
         * has reached its maximum dimension. 
         * 
         * When Arnoldi iteration is terminated early
         * 
         */

        EigenResult diagonalize(const std::vector<Complex>& initial_rho, int k, bool descending = false, double tol = 1e-12);


    protected:
        /**
         * @brief Constructs a Hessenberg matrix from the Lindbladian via Arnoldi iteration.
         * 
         * Builds a Krylov subspace {rho, L(rho), L²(rho), ..., L^(k-1)(rho)} where L is the
         * Lindbladian superoperator. At each iteration, applies the Lindbladian to the most
         * recent Krylov vector, orthogonalizes against all previous vectors using modified
         * Gram-Schmidt, and stores the resulting coefficients in the Hessenberg matrix.
         * 
         * The Hessenberg matrix is an upper Hessenberg form that represents the projection
         * of the Lindbladian onto the Krylov subspace.
         * 
         * @param initial_rho Initial density matrix (flattened) to seed the Krylov subspace
         * @param k Dimension of the Krylov subspace
         * @param tol Tolerance for determining linear dependence in orthogonalization
         * 
         * @return Vectorized Hessenberg matrix stored in row-major order
         * 
         */

        std::vector<double> constructHessenbergMatrix(const std::vector<Complex>& initial_rho, int k, double tol = 1e-12); 

        /**
         * @brief Computes the Frobenius inner product of two matrices.
         * 
         * Calculates the trace of the conjugate transpose product: Tr(A† * B).
         * This is used for orthogonalization in the Arnoldi iteration.
         * 
         * @param matA First matrix (flattened, row-major order)
         * @param matB Second matrix (flattened, row-major order)
         * @param M Dimension of the square matrices (M x M)
         * @param tol Tolerance for considering values as zero
         * 
         * @return Real part of the inner product Tr(A† * B)
         */

        double computeInnerProduct(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M, double tol = 1e-12);
    
        /**
         * @brief Sorts eigenvalue-eigenvector pairs.
         * 
         * Reorders eigenvalues and their corresponding eigenvectors either in ascending
         * or descending order based on the magnitude of the eigenvalues. Used to identify
         * dominant decay modes or steady states of the open quantum system.
         * 
         * @param eigenvalues Vector of eigenvalues, reordered in-place
         * @param eigenvectors Matrix of eigenvectors in row-major order, reordered in-place
         * @param row Number of rows in the eigenvector matrix
         * @param col Number of columns (number of eigenvectors)
         * @param descending If true, sort by decreasing magnitude; otherwise by increasing magnitude
         */
        
        void sortEigenPairs(std::vector<Complex>& eigenvalues,std::vector<Complex>& eigenvectors,int row, int col, bool descending = false);
    
};

#endif