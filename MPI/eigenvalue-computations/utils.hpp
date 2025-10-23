#ifndef UTILS_HPP
#define UTILS_HPP

#include <iomanip>
#include <vector>
#include <iostream>
#include <complex>
#include <Dense>
#include <eigen>
#include <cmath>

using Complex = std::complex<double>;

bool complexEqual(const Complex& x, const Complex& y, double tol = 1e-12);
bool vectorEqual(const std::vector<Complex>& A, const std::vector<Complex>& B,double tol = 1e-12);

std::vector<Complex> matrixMultiplication(const std::vector<Complex>& A,
                                          const std::vector<Complex>& B,
                                          int M, int N, int K);

std::vector<Complex> kroneckerProduct(const std::vector<Complex>& A, int N,
                                      const std::vector<Complex>& B, int M);

std::vector<Complex> addMatrix(const std::vector<Complex>& A,
                               const std::vector<Complex>& B,
                               int N);

template<typename T>
void printMatrix(const std::vector<T>& M, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            std::cout << M[i * col + j] << " ";
        }
        std::cout << std::endl;
    }
}
Complex computeTrace(const std::vector<Complex>& A, int N);

bool checkHermicity(const std::vector<Complex>& A, int N, double tol = 1e-10);

// bool checkPositiveSemiDefinitivity(const std::vector<Complex>& A, int N, double tol = 1e-12);

Eigen::VectorXd findEigenvalues(const std::vector<Complex>& A, int N);

std::vector<Complex> hermitian(const std::vector<Complex>& A, int N);

void printResult(bool condition, const std::string& testName);

#endif