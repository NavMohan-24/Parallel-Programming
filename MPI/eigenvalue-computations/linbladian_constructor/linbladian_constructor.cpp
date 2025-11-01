#include "linbladian_constructor.hpp"
#include "utils.hpp"


//_______________helper-functions_______________________________________________________

// inline function to convert the matrix elements less than tolerance to zero.
inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    double re = (std::abs(z.real()) < tol) ? 0: z.real();
    double img = (std::abs(z.imag()) < tol) ? 0: z.imag();
    return {re, img};
}


//____________________Linbladian-Constructor-Implementation________________________________________________

LinbladianConstructor::LinbladianConstructor(const std::vector<Complex>& hamiltonian_,
                int N_, double rate_, DecayType decay_, Scope scope_)
                : hamiltonian(hamiltonian_), N(N_), num_states(1 << N),
                rate(rate_), decay(decay_), scope(scope_){}
                // validate the hamiltonian in the constructor


std::vector<Complex> LinbladianConstructor::applyLinbladian(std::vector<Complex>& rho)
{
    // validate the rho
    enforceHermiticity(rho,num_states);

    std::vector<Complex> comm = applyCommutator(hamiltonian,rho, num_states);
    std::vector<Complex> diss = constructDissipator(rho, N, num_states, rate, decay, scope);
    std::vector<Complex> linbladian(num_states*num_states);

    Complex i_unit(0.0, -1.0);

    for (int i = 0; i < num_states*num_states; i++){
        linbladian[i] =  cleanMatrixElements(i_unit*comm[i]+diss[i]);
    };

    return linbladian;
};  

std::vector<Complex> LinbladianConstructor::applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M)
{

    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};

    std::vector<Complex> AB(M*M);
    std::vector<Complex> BA(M*M);
    // std::vector<Complex> C(num_states*num_states);

    // AB = matA*matB (H * rho)
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matA.data(),
        M, matB.data(), M, &beta, AB.data(), M
    );

    // BA = matB*matA (rho * H) 
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matB.data(),
        M, matA.data(), M, &beta, BA.data(), M
    );

    // Commutator = AB - BA
    for (int i = 0; i < M*M; i++){
        AB[i] -= BA[i];
    };

    return AB;

};

std::vector<Complex> LinbladianConstructor::applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M)
{

    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};

    std::vector<Complex> AB(M*M);
    std::vector<Complex> BA(M*M);
    

    // A = matA * matB
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matA.data(),
        M, matB.data(), M, &beta, AB.data(), M
    );

    // A = matB * matA
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matB.data(),
        M, matA.data(), M, &beta, BA.data(), M
    );

    // Commutator = AB + BA
    for (int i = 0; i < M*M; i++){
        AB[i] += BA[i];
    };

    return AB;
};

std::vector<Complex> LinbladianConstructor::constructDissipator(const std::vector<Complex>& rho,int N, 
    int num_states, double rate, DecayType decay, Scope scope)
    {
    // D(ρ)=∑_k​(J_k ​ρ J_k^†​ − 1/2{J_k^†​ J_k​,ρ})

    // TODO need to optimise Dissipator construction remove additional matrices

    std::vector<Complex> Dissipator(num_states*num_states, Complex(0,0));
    Complex beta =  {0.0,0.0};

    for (int k = 1; k <= N; k++){

        std::vector<Complex> Jk = constructJumpOperator(N,num_states,k,rate,decay,scope);

        //compute Jk^†​Jk
        Complex alpha = {0.5,0.0}; // 1/2 infront of anti commutator

        std::vector<Complex> K(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    Jk.data(), num_states, 
                    &beta, K.data(), num_states
                );

        // compute  1/2{Jk^†​ Jk​,ρ})
        std::vector<Complex> SP = applyAntiCommutator(K,rho,num_states);

        //compute L_k ​ρ
        alpha = {1.0,0.0};
        std::vector<Complex> tempA(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    rho.data(), num_states, 
                    &beta, tempA.data(), num_states
                );

        //compute L_k ​ρ L_k^†
        std::vector<Complex> FP(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
                    num_states, num_states, num_states,
                    &alpha, tempA.data(), num_states,
                    Jk.data(), num_states, 
                    &beta, FP.data(), num_states
                );

        for (int i = 0; i < num_states*num_states; i++){
            Dissipator[i] += FP[i] - SP[i];
        };
            
    };
    return Dissipator;
};

void LinbladianConstructor::enforceHermiticity(std::vector<Complex>& mat, int M) {
    for (int i = 0; i < M; i++) {
        // Diagonal must be real
        mat[i*M + i] = Complex(mat[i*M + i].real(), 0.0);
        
        // Off-diagonal: symmetrize
        for (int j = i+1; j < M; j++) {
            Complex m_ij = mat[i*M + j];
            Complex m_ji = mat[j*M + i];
            
            // Average: (m_ij + conj(m_ji))/2
            Complex sym_val = 0.5 * (m_ij + std::conj(m_ji));
            
            mat[i*M + j] = sym_val;
            mat[j*M + i] = std::conj(sym_val);
        }
    }
}
