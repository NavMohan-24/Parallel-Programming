#include "ising.hpp"
#include <vector>
#include <cblas.h>


std::vector<double> applyCommutator(std::vector<double> H, std::vector<double> rho, int N){
    
    double alpha = 1.0, beta = 0.0;
    int num_states = 1 << N;

    std::vector<double> A(num_states*num_states);
    std::vector<double> B(num_states*num_states);
    std::vector<double> C(num_states*num_states);

    // A = H * rho
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, alpha, H.data(),
        num_states, rho.data(), num_states, beta, A.data(), num_states
    );

    // B = rho * H
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, alpha, rho.data(),
        num_states, H.data(), num_states, beta, B.data(), num_states
    );

    // C = A - B
    for (int i = 0; i < num_states*num_states; i++){
        C[i] = A[i] - B[i];
    };

    return C;

}

#ifdef BUILD_MAIN
int main(){

    int N = 3;
    int num_states = 1 << N;
    double J = 1.0; 
    double h = 1.0;

    std::vector<double> rho(num_states*num_states);

    std::cout << "Printing rho matrix.." << std::endl;
    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            rho[i*num_states+j] = (i==j) ? 1:0;
            std::cout << rho[i*num_states+j] << " ";
        }
        std::cout << std::endl;
    }

    std::vector<double> H = tfimHamiltonian(N,J,h);
    // std::vector<double> rho = constructRandomRho(N);
    std::vector<double> matrix = applyCommutator(H,rho,N);

    
    
    std::cout << "Printing commutator matrix.." << std::endl;
    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << matrix[i*num_states+j] << " "; 
        }
         std::cout << std::endl;
    }

}
#endif