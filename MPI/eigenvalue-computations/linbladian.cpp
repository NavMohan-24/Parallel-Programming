#include "linbladian.hpp"



LinbladianSolver::LinbladianSolver(const std::vector<Complex>& hamiltonian_,
                int N_, double rate_, DecayType decay_, Scope scope_)
                : hamiltonian(hamiltonian_), N(N_), num_states(1<< N_),
                rate(rate_), decay(decay_), scope(scope_){}
                // validate the hamiltonian in the constructor


std::vector<Complex> LinbladianSolver::solve(const std::vector<Complex>& rho)
{
    // validate the rho

    std::vector<Complex> comm = applyCommutator(hamiltonian,rho, num_states);
    std::vector<Complex> diss = constructDissipator(rho, N, num_states, rate, decay, scope);
    std::vector<Complex> linbladian(num_states*num_states);

    Complex i_unit(0.0, 1.0);

    for (int i = 0; i < num_states*num_states; i++){
        linbladian[i] =  -i_unit*comm[i]- diss[i];
    };

    return linbladian;
};        


std::vector<Complex> LinbladianSolver::applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int num_states)
{

    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};

    std::vector<Complex> AB(num_states*num_states);
    std::vector<Complex> BA(num_states*num_states);
    // std::vector<Complex> C(num_states*num_states);

    // AB = matA*matB (H * rho)
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matA.data(),
        num_states, matB.data(), num_states, &beta, AB.data(), 
        num_states
    );

    // BA = matB*matA (rho * H) 
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matB.data(),
        num_states, matA.data(), num_states, &beta, BA.data(),
        num_states
    );

    // Commutator = AB - BA
    for (int i = 0; i < num_states*num_states; i++){
        AB[i] -= BA[i];
    };

    return AB;

};

std::vector<Complex> LinbladianSolver::applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int num_states)
{

    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};

    std::vector<Complex> AB(num_states*num_states);
    std::vector<Complex> BA(num_states*num_states);
    

    // A = matA * matB
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matA.data(),
        num_states, matB.data(), num_states, &beta, AB.data(), 
        num_states
    );
    // A = matB * matA
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matB.data(),
        num_states, matA.data(), num_states, &beta, BA.data(),
        num_states
    );
    // Commutator = AB + BA
    for (int i = 0; i < num_states*num_states; i++){
        AB[i] += BA[i];
    };

    return AB;
};

std::vector<Complex> LinbladianSolver::constructDissipator(const std::vector<Complex>& mat,int N, 
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
        cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    Jk.data(), num_states, 
                    &beta, K.data(), num_states
                );

        // compute  1/2{Jk^†​ Jk​,ρ})
        std::vector<Complex> SP = applyAntiCommutator(K,mat,num_states);

        //compute L_k ​ρ
        alpha = {1.0,0.0};
        std::vector<Complex> tempA(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    mat.data(), num_states, 
                    &beta, tempA.data(), num_states
                );

        //compute L_k ​ρ L_k^†
        std::vector<Complex> FP(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
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

