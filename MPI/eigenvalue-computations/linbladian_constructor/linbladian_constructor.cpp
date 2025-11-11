#include "linbladian_constructor.hpp"
#include "utils.hpp"


//=======================================================================================
// Helper Functions
//=======================================================================================

/**
 * @brief Cleans matrix elements by setting values below tolerance to zero.
 * 
 * @param z Complex number to clean
 * @param tol Tolerance threshold for considering values as zero
 * @return Complex number with small components zeroed out
 */
inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    double re = (std::abs(z.real()) < tol) ? 0: z.real();
    double img = (std::abs(z.imag()) < tol) ? 0: z.imag();
    return {re, img};
}


//=======================================================================================
// Lindbladian Constructor Implementation
//=======================================================================================

LinbladianConstructor::LinbladianConstructor(const std::vector<Complex>& hamiltonian_,
                int N_, double rate_, DecayType decay_, Scope scope_)
                : hamiltonian(hamiltonian_), N(N_), num_states(1 << N_),
                rate(rate_), decay(decay_), scope(scope_){}


std::vector<Complex> LinbladianConstructor::applyLinbladian(std::vector<Complex>& rho)
{
    // Enforce Hermiticity to correct numerical errors
    enforceHermiticity(rho, num_states);


    // Compute unitary part: -i[H, ρ]
    std::vector<Complex> comm = applyCommutator(hamiltonian, rho, num_states);
    
    // Compute dissipative part: Σₖ(Lₖ ρ Lₖ† - ½{Lₖ†Lₖ, ρ})
    std::vector<Complex> diss = constructDissipator(rho, N, num_states, rate, decay, scope);
    
    
    // Combine unitary and dissipative contributions
    std::vector<Complex> linbladian(num_states*num_states);
    Complex i_unit(0.0, -1.0);

    #pragma omp parallel for if (num_states > std::pow(2,7))
    for (int i = 0; i < num_states*num_states; i++){
        linbladian[i] = cleanMatrixElements(i_unit*comm[i] + diss[i]);
    }

    return linbladian;
}

std::vector<Complex> LinbladianConstructor::applyCommutator(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M)
{
    Complex alpha = {1.0, 0.0};
    Complex beta = {0.0, 0.0};

    std::vector<Complex> AB(M*M);
    std::vector<Complex> BA(M*M);

    // Compute AB = matA * matB
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, 
        matA.data(), M, 
        matB.data(), M, 
        &beta, AB.data(), M
    );

    // Compute BA = matB * matA
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, 
        matB.data(), M, 
        matA.data(), M, 
        &beta, BA.data(), M
    );

    // Compute commutator: [A, B] = AB - BA
    #pragma omp parallel for if (M > std::pow(2,7))
    for (int i = 0; i < M*M; i++){
        AB[i] -= BA[i];
    }

    return AB;
}

std::vector<Complex> LinbladianConstructor::applyAntiCommutator(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M)
{
    Complex alpha = {1.0, 0.0};
    Complex beta = {0.0, 0.0};

    std::vector<Complex> AB(M*M);
    std::vector<Complex> BA(M*M);

    // Compute AB = matA * matB
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, 
        matA.data(), M, 
        matB.data(), M, 
        &beta, AB.data(), M
    );

    // Compute BA = matB * matA
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, 
        matB.data(), M, 
        matA.data(), M, 
        &beta, BA.data(), M
    );

    // Compute anti-commutator: {A, B} = AB + BA
    #pragma omp parallel for if (M > std::pow(2,7)) 
    for (int i = 0; i < M*M; i++){
        AB[i] += BA[i];
    }

    return AB;
}

std::vector<Complex> LinbladianConstructor::constructDissipator(const std::vector<Complex>& rho, int N, 
    int num_states, double rate, DecayType decay, Scope scope)
{
    // Compute dissipator: D(ρ) = Σₖ(Jₖ ρ Jₖ† - ½{Jₖ†Jₖ, ρ})
    
    std::vector<Complex> Dissipator(num_states*num_states, Complex(0, 0));
    Complex beta = {0.0, 0.0};

    // Sum over all jump operators
    #pragma omp parallel for schedule(dynamic) 
    for (int k = 1; k <= N; k++){

        // Construct jump operator Jₖ for site k
        std::vector<Complex> Jk = constructJumpOperator(N, num_states, k, rate, decay, scope);

        // Compute Jₖ†Jₖ with factor of 1/2 for anti-commutator term
        Complex alpha = {0.5, 0.0};
        std::vector<Complex> K(num_states*num_states, Complex(0, 0));
        cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    Jk.data(), num_states, 
                    &beta, K.data(), num_states
                );

        // Compute anti-commutator term: ½{Jₖ†Jₖ, ρ}
        std::vector<Complex> SP = applyAntiCommutator(K, rho, num_states);

        // Compute first part: Jₖ ρ
        alpha = {1.0, 0.0};
        std::vector<Complex> tempA(num_states*num_states, Complex(0, 0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    rho.data(), num_states, 
                    &beta, tempA.data(), num_states
                );

        // Compute full positive term: Jₖ ρ Jₖ†
        std::vector<Complex> FP(num_states*num_states, Complex(0, 0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                    num_states, num_states, num_states,
                    &alpha, tempA.data(), num_states,
                    Jk.data(), num_states, 
                    &beta, FP.data(), num_states
                );

        // Accumulate contribution: (Jₖ ρ Jₖ† - ½{Jₖ†Jₖ, ρ})
        #pragma omp critical{
        for (int i = 0; i < num_states*num_states; i++){
            Dissipator[i] += FP[i] - SP[i];
        }
    }
    }
    
    return Dissipator;
}

void LinbladianConstructor::enforceHermiticity(std::vector<Complex>& mat, int M) {
    // Symmetrize matrix to enforce Hermiticity: mat = (mat + mat†) / 2
    
    #pragma omp parallel for if (M > std::pow(2,7))
    for (int i = 0; i < M; i++) {
        // Diagonal elements must be real
        mat[i*M + i] = Complex(mat[i*M + i].real(), 0.0);
        
        // Symmetrize off-diagonal elements
        for (int j = i+1; j < M; j++) {
            Complex m_ij = mat[i*M + j];
            Complex m_ji = mat[j*M + i];
            
            // Average to enforce Hermiticity: (m_ij + conj(m_ji)) / 2
            Complex sym_val = 0.5 * (m_ij + std::conj(m_ji));
            
            mat[i*M + j] = sym_val;
            mat[j*M + i] = std::conj(sym_val);
        }
    }
}