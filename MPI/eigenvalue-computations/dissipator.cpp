#include <stdio.h>
#include <vector>
#include <cblas.h>
#include <iomanip>


#include "tfim.hpp"
#include "random_rho.hpp"
#include "jump_op.hpp"

using Complex = std::complex<double>;

std::vector<Complex> applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int num_states){
    
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

std::vector<Complex> constructDissipator(const std::vector<Complex>& rho, int N, int num_states, double rate = 1.0){
    // D(ρ)=∑_k​(L_k ​ρ L_k^†​ − 1/2{L_k^†​ L_k​,ρ})
   
    std::vector<Complex> Dissipator(num_states*num_states, Complex(0,0));
    Complex beta =  {0.0,0.0};

    for (int k = 1; k <= N; k++){

        std::vector<Complex> Jk = constructJumpOperator(N,num_states,k,rate,DecayType::SigmaMinus,Scope::Local);

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
}


#ifdef BUILD_MAIN
bool complexEqual(const Complex& x, const Complex& y, double tol = 1e-12){
    return std::abs(x-y) < tol;
}

bool vectorEqual(const std::vector<Complex>& A, const std::vector<Complex>& B,double eps = 1e-12){
    if (A.size() != B.size()){return false;}
    for (int i = 0; i < A.size(); i++){
        if (!complexEqual(A[i], B[i], eps)){return false;}
    }
    return true;
}

Complex computeTrace(const std::vector<Complex>& A, int N){

    //assumes a square matrix
    Complex trace = (0.0,0.0);
    for (int i = 0; i < N; i++){
        trace += A[i*N+i];
    }
    return trace;
}

bool checkHermicity(const std::vector<Complex>& A, int N, double tol = 1e-12){

    for (int i = 0; i < N; i++){
       for (int j = 0; j < N; j++){
           if (std::abs(A[i*N+j] - std::conj(A[j*N+i])) > tol) {return false;} 
       } 
    }
    return true;
}   


int main(){

    // coherent evolution test
    int N = 3;
    int num_states = 1 << N;
    double k = 0.0;

    std::vector<Complex> rho = createRandomRho(num_states);
    std::vector<Complex> Dissipator = constructDissipator(rho, N, num_states,k);
    std::vector<Complex> ref(num_states*num_states, Complex(0,0));

    if (vectorEqual(Dissipator, ref)){
        std::cout<<"\n Dissipator is Zero Matrix for rate = 0..✅" << "\n" <<std::endl;
    }
    else{
        std::cout<<"\n Dissipator is not Zero Matrix for rate = 0..❌" << "\n" <<std::endl;
    }
    
    //unit test

    N = 1;
    num_states = 1 << N;
    k = 1.0;

    std::vector<Complex> psi = {
        Complex (1.0 / std::sqrt(2.0)),
        Complex (1.0 / std::sqrt(2.0)) 
    };
 
    rho = createRhoFromStatevector(psi);
    Dissipator = constructDissipator(rho, N, num_states,k);

    ref = {0.5, -0.25, -0.25, -0.5};

    if (vectorEqual(Dissipator,ref)){
        std::cout<<"\n Passed Unit test for Dissipator..✅" << "\n" <<std::endl;
    }
    else{
        std::cout<<"\n Failed Unit test for Dissipator..❌" << "\n" <<std::endl;
    }

    // testing properties of dissipator
    N = 3;
    num_states = 1 << N;
    k = 1.0;

    rho = createRandomRho(num_states);
    Dissipator = constructDissipator(rho, N, num_states, k);

    Complex trace = computeTrace(Dissipator,num_states);
    double tol = 1e-12;

    if (trace.imag() == 0.0 && std::abs(trace)< tol){
        std::cout<<"\n Trace of the Disspator is Zero..✅" << "\n" <<std::endl;
    }
    else{
        std::cout<<"\n Trace of the Disspator is NOT Zero..❌" << "\n" <<std::endl;
    }

    if (checkHermicity(Dissipator, num_states)){
        std::cout<<"\n Dissipator is Hermitian..✅" << "\n" <<std::endl;
    }
    else{
        std::cout<<"\n Dissipator is NOT Hermitian..❌" << "\n" <<std::endl;
    }
}
#endif


