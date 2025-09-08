#include "random_rho.hpp"


Complex genMatrixElements(std::mt19937 &gen, std::normal_distribution<double> &dist){
    return {dist(gen), dist(gen)};
}

inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    double re = (std::abs(z.real()) < tol) ? 0: z.real();
    double img = (std::abs(z.imag()) < tol) ? 0: z.imag();
    return {re, img};
}

std::vector<Complex> createRandomRho(int num_states){

    // int num_states = 1 << N;

    std::vector<Complex> A(num_states*num_states);
    std::vector<Complex> rho(num_states*num_states);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0,1.0);

    for (int i = 0; i < num_states ; i++){
        for (int j = 0; j < num_states ; j++){
            A[i*num_states+j] = genMatrixElements(gen,dist);
        }
    }
    Complex alpha = {1.0, 0.0}; // make alpha and beta complex and pass by reference
    Complex beta  = {0.0, 0.0};

    // multiply AA^{H} 
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasConjTrans,
        num_states,num_states,num_states,&alpha,A.data(), 
        num_states,A.data(), num_states,&beta,rho.data(),
        num_states
    ); 
    
    
    double trace = 0.0;
    for (int i = 0; i < num_states; i++){
            trace += std::real(rho[i*num_states+i]); // diagonal values are guarnteed to be real for hermitian.
        
    }
    for (int i =0; i < num_states*num_states; i++){
        rho[i] = cleanMatrixElements(rho[i]/trace);
    }


    return rho;
}


#ifdef BUILD_MAIN
#include <iomanip>
#include <Dense>
int main() {

    int N = 3;
    int num_states = 1 << N;

    std::vector<Complex> matrix = createRandomRho(num_states);
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

    std::cout<< "Checking the positive semi-definity of the dentsity matrix.. 🔬" << "\n" << std::endl;
    
    bool isPositiveSemidefinite = true;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(eigen_matrix);
    Eigen::VectorXd eigenvalues = es.eigenvalues();

    for (int i = 0; i < eigenvalues.size(); i++){
        if (eigenvalues[i] < -tol){
            isPositiveSemidefinite = false;
            break;
        }
    }
    if (isPositiveSemidefinite){
        std::cout << "Matrix is Positive Semi-Definite ✅" << "\n" << std::endl;
    }
    else{
        std::cout << "Matrix is NOT Positive Semi-Definite ❌" << "\n" << std::endl;
    }
}
#endif