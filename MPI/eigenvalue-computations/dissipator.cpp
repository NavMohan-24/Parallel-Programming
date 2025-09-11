#include <stdio.h>
#include "tfim.hpp"
#include "random_rho.hpp"
#include <vector>
#include <cblas.h>


using Complex = std::complex<double>;

std::vector<Complex> applyAntiCommutator(std::vector<Complex>& matA, std::vector<Complex>& matB, int num_states){
    
    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};
    // int num_states = 1 << N;

    std::vector<Complex> A(num_states*num_states);
    std::vector<Complex> B(num_states*num_states);
    std::vector<Complex> C(num_states*num_states);

    // A = matA * matB
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matA.data(),
        num_states, matB.data(), num_states, &beta, A.data(), 
        num_states
    );
    // A = matB * matA
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matB.data(),
        num_states, matA.data(), num_states, &beta, B.data(),
        num_states
    );
    // C = A + B
    for (int i = 0; i < num_states*num_states; i++){
        C[i] = A[i] + B[i];
    };

    return C;

}

std::vector<Complex> constructDissipator(std::vector<Complex>&, rho, std::vector<Complex>&, J, Complex r){
    
}



