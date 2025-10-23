#include "utils.hpp"


bool complexEqual(const Complex& x, const Complex& y, double tol) {
    return std::abs(x-y) < tol;
    };

bool vectorEqual(const std::vector<Complex>& A, const std::vector<Complex>& B,double tol){
    if (A.size() != B.size()){return false;}
    for (int i = 0; i < A.size(); i++){
        if (!complexEqual(A[i], B[i], tol)){return false;}
    }
    return true;
}

//matrix multiplication
std::vector<Complex> matrixMultiplication(const std::vector<Complex>& A,
                                          const std::vector<Complex>& B,
                                          int M, int N, int K) 
{
    std::vector<Complex> mat(M * N, Complex(0,0));

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < K; k++) {
                mat[i*N + j] += A[i*K + k] * B[k*N + j];
            }
        }
    }

    return mat;
}

// Kronecker product for square matrices (N×N and M×M)
std::vector<Complex> kroneckerProduct(const std::vector<Complex>& A, int N,
                                      const std::vector<Complex>& B, int M) {
    int size = N * M;  // new dimension
    std::vector<Complex> result(size * size, Complex(0,0));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Complex a = A[i*N + j];  // element A(i,j)

            for (int k = 0; k < M; k++) {
                for (int l = 0; l < M; l++) {
                    // map to new indices
                    int row = i*M + k;
                    int col = j*M + l;
                    result[row*size + col] = a * B[k*M + l];
                }
            }
        }
    }
    return result;
}

std::vector<Complex> addMatrix(const std::vector<Complex>& A,
                               const std::vector<Complex>& B,
                               int N) {
    std::vector<Complex> C(N * N, Complex(0, 0));

    for (int i = 0; i < N * N; i++) {
        C[i] = A[i] + B[i];
    }

    return C;
}

// template<typename T>
// void printMatrix(const std::vector<T>& M, int N) {
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             std::cout << std::setw(10) << std::fixed << std::setprecision(2) 
//                       << M[i * N + j] << " ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }

Complex computeTrace(const std::vector<Complex>& A, int N){

    //assumes a square matrix
    Complex trace = (0.0,0.0);
    for (int i = 0; i < N; i++){
        trace += A[i*N+i];
    }
    //std::cout << std::fixed << std::setprecision(3) << "Trace of the matrix : " << trace  << "\n"  << std::endl;
    return trace;
}

bool checkHermicity(const std::vector<Complex>& A, int N, double tol){

    for (int i = 0; i < N; i++){
       for (int j = 0; j < N; j++){
           if (std::abs(A[i*N+j] - std::conj(A[j*N+i])) > tol) {return false;} 
       } 
    }
    return true;
}   

// bool checkPositiveSemiDefinitivity(const std::vector<Complex>& A, int N, double tol){
//     Eigen::MatrixXcd eigen_matrix(N,N);
//     // Copy flat vector into Eigen matrix
//     for (int i = 0; i < N; i++){
//         for (int j = 0; j < N; j++){
//             eigen_matrix(i,j) = A[i*N + j];  // row-major
//         }
//     }
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(eigen_matrix);
//     Eigen::VectorXd eigenvalues = es.eigenvalues();
//     std::cout << "Eigenvalues: " << eigenvalues.transpose() << "\n";
//     for (int i = 0; i < eigenvalues.size(); i++){
//         if (eigenvalues[i] < -tol) return false; // allow small negative tol due to FP
//     }
//     return true;
// }

Eigen::VectorXd findEigenvalues(const std::vector<Complex>& A, int N){
    Eigen::MatrixXcd eigen_matrix(N,N);
   
    // Copy flat vector into Eigen matrix
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            eigen_matrix(i,j) = A[i*N + j];  // row-major
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(eigen_matrix);
    Eigen::VectorXd eigenvalues = es.eigenvalues();

    return eigenvalues;

}

std::vector<Complex> hermitian(const std::vector<Complex>& A, int N){
    std::vector<Complex> B(N * N, Complex(0, 0));

    for (int i = 0; i < N; i++){
        for (int j=0; j < N; j++){
            B[j*N+i] = std::conj(A[i*N+j]);
        }
    }
    return B;
}

void printResult(bool condition, const std::string& testName) {
        
    if (condition)
        std::cout << testName << " ✅" << std::endl;
    else
        std::cout << testName << " ❌" << std::endl;
    
    
    }
