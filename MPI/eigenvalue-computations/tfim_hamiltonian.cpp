#include "ising.hpp"
#include <boost/dynamic_bitset.hpp>
#include <iomanip>

std::vector<Complex>  tfimHamiltonian(int num_states, int N, double J, double h){

    // int num_states = (int) std::pow(2, N); - risky and can have floating imprecision

    // int num_states = 1 << N;

    std::vector<Complex> hamiltonian(num_states*num_states);

    for (int k = 0 ; k < num_states; k++){
        for (int l = 0; l < num_states; l++){
            Complex ham(0.0,0.0);
            if (l == k){
                boost::dynamic_bitset<> b(N, k);
                for (int i = 0; i < N - 1; i++){
                    ham += Complex((b[i]==b[i+1]) ? -J:J,0);
                }
            }
            else{
                boost::dynamic_bitset<> bl(N, l); 
                boost::dynamic_bitset<> bk(N, k);
                boost::dynamic_bitset<> c = bl ^ bk;
                ham = Complex((c.count()==1) ? -h:0,0); 
            }
            hamiltonian[k*num_states+l] = ham;
        }
    }
    return hamiltonian;

}

#ifdef BUILD_MAIN
#include <iomanip>
#include <Dense>
int main(){

    int N = 3;
    int num_states = 1 << N;
    double J = 1.0;
    double h = 1.0;

    std::vector<Complex> matrix = tfimHamiltonian(num_states, N, J, h);
    Eigen::MatrixXcd eigen_matrix(num_states, num_states);

    Complex trace = {0.0,0.0};

    std::cout << "Printing density matrix.." << "\n" << std::endl;
    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << std::fixed << std::setprecision(3) << matrix[i*num_states+j] << "\t"; 
            eigen_matrix(i,j) = matrix[i*num_states+j];
        }
        std::cout << std::endl;
        trace += matrix[i*num_states+i];
    }

    std::cout<< "\n" << "Checking the Trace of the dentsity matrix.. 🔬" << "\n" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "Trace of the matrix : " << trace  << "\n"  << std::endl;

    std::cout<< "Checking the Hermicity of the density matrix.. 🔬" << "\n" << std::endl;

    bool isHermitian = true;
    double tol = 1e-12;

    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            if (std::abs(matrix[i*num_states + j] - std::conj(matrix[j*num_states + i])) > tol) {
                isHermitian = false;
                break;
            }
        }
        if (!isHermitian) break;
    }

    if (isHermitian){
        std::cout << "Matrix is Hermitian ✅" << "\n" << std::endl;
    }
    else{
        std::cout << "Matrix is NOT Hermitian ❌" << "\n" << std::endl;
    }

    std::cout<< "Checking the sum of eigenvalues... 🔬" << "\n" << std::endl;
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(eigen_matrix);
    Eigen::VectorXd eigenvalues = es.eigenvalues();
    Complex sum_eigs = es.eigenvalues().sum();
    if (std::abs(sum_eigs.real()) < 1e-12) sum_eigs.real(0.0);
    if (std::abs(sum_eigs.imag()) < 1e-12) sum_eigs.imag(0.0);
    if (trace==sum_eigs){
        std::cout << "Sum of Eigenvalues = Sum of Trace ✅" << std::endl;
    }
    else{
        std::cout << "Sum of Eigenvalues = Sum of Trace ❌" << std::endl;
    }

}
#endif