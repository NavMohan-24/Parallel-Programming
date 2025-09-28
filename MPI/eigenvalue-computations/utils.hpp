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

void printMatrix(const std::vector<Complex>& M, int N);

Complex computeTrace(const std::vector<Complex>& A, int N);

bool checkHermicity(const std::vector<Complex>& A, int N, double tol = 1e-12);

// bool checkPositiveSemiDefinitivity(const std::vector<Complex>& A, int N, double tol = 1e-12);

Eigen::VectorXd findEigenvalues(const std::vector<Complex>& A, int N);

std::vector<Complex> hermitian(const std::vector<Complex>& A, int N);

void printResult(bool condition, const std::string& testName);

#endif