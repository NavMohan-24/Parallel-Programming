#include "tfim.hpp"
#include "random_rho.hpp"
#include <vector>
#include <cblas.h>


using Complex = std::complex<double>;


std::vector<Complex> applyCommutator(std::vector<Complex>& matA, std::vector<Complex>& matB, int num_states){
    
    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};
    // int num_states = 1 << N;

    std::vector<Complex> A(num_states*num_states);
    std::vector<Complex> B(num_states*num_states);
    std::vector<Complex> C(num_states*num_states);

    // A = H * rho
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matA.data(),
        num_states, matB.data(), num_states, &beta, A.data(), 
        num_states
    );

    // B = rho * H 
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states, num_states, num_states, &alpha, matB.data(),
        num_states, matA.data(), num_states, &beta, B.data(),
        num_states
    );

    // C = A - B
    for (int i = 0; i < num_states*num_states; i++){
        C[i] = A[i] - B[i];
    };

    return C;

}

#ifdef BUILD_MAIN
#include <iomanip>
int main(){

    int N = 3;
    int num_states = 1 << N;
    double J = 1.0; 
    double h = 1.0;
    double tol = 1e-12;

    std::vector<Complex> psi = {
        Complex (1.0 / std::sqrt(2.0)),
        Complex (0.0),
        Complex (0.0), 
        Complex (0.0), 
        Complex (0.0), 
        Complex (0.0), 
        Complex (0.0),
        Complex (1.0 / std::sqrt(2.0)) 
    };


    std::vector<Complex> H = tfimHamiltonian(num_states, N, J,h);
    // std::vector<Complex> rho = createRandomRho(num_states);
    std::vector<Complex> rho = createRhoFromStatevector(psi);
    std::vector<Complex> matrix = applyCommutator(H,rho,num_states);
    std::vector<Complex> anti_matrix = applyCommutator(rho,H,num_states);
    
    bool isThereAntisymmetry = true;
    Complex trace = {0.0, 0.0};

    std::cout << "Printing commutator matrix.." << std::endl;
    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            if(std::abs(matrix[i*num_states+j]+anti_matrix[i*num_states+j])>tol){
                isThereAntisymmetry = false;
            }
            std::cout << std::fixed << std::setprecision(3) << matrix[i*num_states+j] << "\t"; 
        }
        trace += matrix[i*num_states+i];
        std::cout << std::endl;
    }

    std::cout << "\n" <<"Checking the Anti-Hermicity of the Commutator..🔬" << "\n" << std::endl;
    if (isThereAntisymmetry){
        std::cout << "Commutator satisfy Anti-Hermiticity ✅" << "\n" << std::endl;
    }
    else{
        std::cout << "Commutator do not satisfy Anti-Heremicity ❌" << "\n" << std::endl;
    }

    std::cout << "Checking the Trace of the Commutator..🔬" << "\n" <<std::endl;
    std::cout << std::fixed << std::setprecision(3) << "Trace of the Commutator : " << trace  << "\n"  << std::endl;

    if (std::abs(trace) < tol){
        std::cout << "Trace is close to zero..✅" << "\n" << std::endl;
    }
    else{
        std::cout << "Trace is not close to zero..❌" << "\n" << std::endl;
    }
}
#endif