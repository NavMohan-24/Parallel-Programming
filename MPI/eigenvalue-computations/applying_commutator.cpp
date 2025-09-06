#include "ising.hpp"
#include "random_rho.hpp"
#include <vector>
#include <cblas.h>
#include <iomanip>

using Complex = std::complex<double>;


std::vector<Complex> applyCommutator(std::vector<Complex> H, std::vector<Complex> rho, int N){
    
    Complex alpha = {1.0,0.0};
    Complex beta =  {1.0,0.0};
    int num_states = 1 << N;

    std::vector<Complex> A(num_states*num_states);
    std::vector<Complex> B(num_states*num_states);
    std::vector<Complex> C(num_states*num_states);

    // A = H * rho
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
        num_states, num_states, num_states, &alpha, H.data(),
        num_states, rho.data(), num_states, &beta, A.data(), 
        num_states
    );

    // B = rho * H
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
        num_states, num_states, num_states, &alpha, rho.data(),
        num_states, H.data(), num_states, &beta, B.data(),
        num_states
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

    // std::vector<Complex> rho(num_states*num_states);

    // std::cout << "Printing rho matrix.." << std::endl;
    // for (int i = 0; i < num_states; i++){
    //     for (int j = 0; j < num_states; j++){
    //         rho[i*num_states+j] = Complex((i==j) ? 1:0,0);
    //         std::cout << rho[i*num_states+j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    std::vector<Complex> H = tfimHamiltonian(N,J,h);
    std::vector<Complex> rho = createRandomRho(N);
    std::vector<Complex> matrix = applyCommutator(H,rho,N);
    
    std::cout << "Printing commutator matrix.." << std::endl;
    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << matrix[i*num_states+j] << " "; 
        }
         std::cout << std::endl;
    }

}
#endif