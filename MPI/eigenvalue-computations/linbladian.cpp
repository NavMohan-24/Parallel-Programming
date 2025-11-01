#include "linbladian.hpp"
#include "utils.hpp"

//_______________helper-functions_______________________________________________________

// inline function to convert the matrix elements less than tolerance to zero.
inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    double re = (std::abs(z.real()) < tol) ? 0: z.real();
    double img = (std::abs(z.imag()) < tol) ? 0: z.imag();
    return {re, img};
}


// helper function to enforce Hermiticity
/**
 * @brief Helper Function to enforce Hermicity of a matrix.
 * 
 * Floating point errors in computation could affect the hermicity 
 * of the matrix. Perform (A+A^H)/2 on an input matrix A.
 * 
 * @param mat Input square matrix.
 * @param M  Dimension of the matrix.
 * @return void
 */

void enforceHermiticity(std::vector<Complex>& mat, int M) {
    for (int i = 0; i < M; i++) {
        // Diagonal must be real
        mat[i*M + i] = Complex(mat[i*M + i].real(), 0.0);
        
        // Off-diagonal: symmetrize
        for (int j = i+1; j < M; j++) {
            Complex m_ij = mat[i*M + j];
            Complex m_ji = mat[j*M + i];
            
            // Average: (m_ij + conj(m_ji))/2
            Complex sym_val = 0.5 * (m_ij + std::conj(m_ji));
            
            mat[i*M + j] = sym_val;
            mat[j*M + i] = std::conj(sym_val);
        }
    }
}

// inline function to compute trace of a matrix.
inline Complex computeTrace(const std::vector<Complex>& mat, int M){

    // assumes a square matrix
    Complex trace(0.0, 0.0);
    for (int i = 0; i < M; i++){
        trace += mat[i*M + i];
    }
    return trace;
    
}

//_______________________Linbladian-Solver-Implementation________________________________________________

LinbladianDiagonalizer::LinbladianDiagonalizer(const std::vector<Complex>& hamiltonian_,
                int N_, int d_, double rate_, DecayType decay_, Scope scope_)
                : hamiltonian(hamiltonian_), N(N_), d(d_), num_states(std::pow(d_,N_)),
                rate(rate_), decay(decay_), scope(scope_){}
                // validate the hamiltonian in the constructor


std::vector<Complex> LinbladianDiagonalizer::applyLinbladian(std::vector<Complex>& rho)
{
    // validate the rho
    enforceHermiticity(rho,num_states);

    std::vector<Complex> comm = applyCommutator(hamiltonian,rho, num_states);
    std::vector<Complex> diss = constructDissipator(rho, N, num_states, rate, decay, scope);
    std::vector<Complex> linbladian(num_states*num_states);

    Complex i_unit(0.0, -1.0);

    for (int i = 0; i < num_states*num_states; i++){
        linbladian[i] =  cleanMatrixElements(i_unit*comm[i]+diss[i]);
    };

    return linbladian;
};        



std::vector<double> LinbladianDiagonalizer::constructHessenbergMatrix(const std::vector<Complex>& initial_rho, int k, double tol)
{     
    //lambda function to calculate index
    auto idx = [k](int i, int j){return i*(k+1)+j;};

    // Hessenberg Matrix Initialization. 
    // Elements are expected to be double since 
    // it is constructed from inner product between two Hermitian matrices
    std::vector<double> H((k+1)*(k+1),0.0); 

    // container to store kyrlov basis vectors
    std::vector<std::vector<Complex>> kyrlov_basis(k+2);
  
    // normalize the input matrix
    kyrlov_basis[0] = initial_rho;
    normalizeMatrix(kyrlov_basis[0],num_states);
    
    for (int j = 0; j <= k; j++){

        // apply the Linbladian superoperator
        kyrlov_basis[j+1] = applyLinbladian(kyrlov_basis[j]);
        
        // orthogonalize density matrix
        for (int i = 0; i <= j; i++)
        {   
            // compute Hessenberg matrix elements
            H[idx(i,j)] = computeInnerProduct(kyrlov_basis[i],kyrlov_basis[j+1], num_states);

            // substract the projection on basis[i]
            for (int m = 0; m < num_states*num_states; m++){
                kyrlov_basis[j+1][m] -= H[idx(i,j)]*kyrlov_basis[i][m];
            }

        }
        if (j < k){

            // compute norm of new basis vector
            double nf = computeInnerProduct(kyrlov_basis[j+1], kyrlov_basis[j+1], num_states);

            // graceful exit when norm is zero
            if (nf < 1e-12){
                
                // resizing Hessenberg matrix 
                std::vector<double> H_new((j+1)*(j+1)); // size of the H matrix depends upon current iteration value.

                // allocating matrix elements from H_new to H
                for (int m=0; m < j+1; m++){
                    for (int n=0; n < j+1; n++){
                        H_new[m*(j+1)+n] = H[m*(k+1)+n]; // j+1^th row and column should be zero 
                    }
                }
                
                // resizing and flatten kyrlov basis vector
                kyrlov_basis_flatten.resize(num_states*num_states*(j+1));
                for (int i = 0; i < j+1; i++){
                    for (int m = 0; m < num_states*num_states; m++){
                        kyrlov_basis_flatten[m*(j+1)+ i] = kyrlov_basis[i][m]; 
                    }
                }

                return H_new;
                // throw std::runtime_error("Arnoldi Method Breakdown : Norm became zero");

            }
            // norm of j+1 th kryrlov vector becomes sub-diagonal elements of H matrix
            H[idx(j+1,j)] = std::sqrt(nf);
            
            

            // normalzing the krylov basis
            for (int i = 0; i < num_states*num_states; i++){
                kyrlov_basis[j+1][i] = kyrlov_basis[j+1][i]/std::sqrt(nf); 
            };
            
            // calculate (A+A^H)/2 to ensure hermicity of density matrix after floating point errors
            // enforceHermiticity(kyrlov_basis[j+1],num_states); // (can move it back to apply linbladian rho)


        }    
       
    }
    //Flattening the kyrlov basis 
    kyrlov_basis_flatten.resize(num_states*num_states*(k+1));
    for (int i = 0; i < k+1; i++){
        for (int j = 0; j < num_states*num_states; j++){
            kyrlov_basis_flatten[j*(k+1) + i] = kyrlov_basis[i][j]; // flattening and storing it in transposed manner
        }
    }

    return H;
}

EigenResult LinbladianDiagonalizer::diagonalize(const std::vector<Complex>& init_rho, int k, bool descending, double tol){

    int max_k = num_states*num_states - 1;
    if (k > max_k){
        throw std::runtime_error("Number of requested eigen vals/vecs must be <= N.");
    }

    std::vector<double> H = constructHessenbergMatrix(init_rho, k, tol);
    EigenResult result;
    
    int size = static_cast<int>(H.size()); 
    int n = std::sqrt(size); // dimension of Hessenberg matrix (k+1) x (k+1)

    int ilo = 1, ihi = n; // lapack allows working with submatrices - here we have chosen full matrix

    std::vector<double> wr(n), wi(n); // containers to store real and imaginary parts of eigenvalues computed.
    std::vector<double> Z(n*n, 0.0); // container to store Schur vectors.
    // Initialize Z as identity matrix
    for (int i = 0; i < n; i++) {
        Z[i*n + i] = 1.0;
    }

    // converting H to Schur form + Schur vectors using QR algorithm
    int info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR,
        'S', 'V',  
        n, ilo, ihi,
        H.data(), n, 
        wr.data(), wi.data(),
        Z.data(), n
    ); 
    
    // checking whether Schur form construction is success
    if (info != 0) {
        throw std::runtime_error("DHSEQR failed in diagonalizing the matrix");
    }

    // adding eigenvalues to output struct - eigen values of Hessenberg matrix approximates that of Linbladian
    result.eigenvalues.resize(n);
    for (int i=0; i < n; i++){
        result.eigenvalues[i] = std::complex<double>(wr[i], wi[i]);
    };
    
    // diagonalizing schur form of H
    std::vector<double> VR(n*n); // container to store right eigenvector of schur form.
    int mm = n; // leading dimension matrix to store right eigen vector (row of VR)
    int m; // number of eigenvectors to be found (coloum of the VR )
    info = LAPACKE_dtrevc(
        LAPACK_ROW_MAJOR,
        'R', 'A',
        nullptr, n,
        H.data(), n,
        nullptr, n,
        VR.data(), n,
        mm, &m
    );

    if (info != 0) {
        throw std::runtime_error("DTREVC failed in diagonalizing the matrix");
    }

    // finding eigenvector of Hessenberg Matrix
    std::vector<double> eigenvecs_H(n*n);
    cblas_dgemm(
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        n, n, n, 1.0, Z.data(), n,
        VR.data(), n, 0.0,
        eigenvecs_H.data(), n 
    );

    // converting eigenvecs_H to complex to do further transformation
    std::vector<Complex> eigenvecs_H_complex(n*n);
    for (int i = 0; i < n*n; i++) {
        eigenvecs_H_complex[i] = Complex(eigenvecs_H[i], 0.0);
    }


    // transforming eigenvectors of Hessenberg matrix to density matrix basis
    result.eigenvectors.resize(num_states*num_states*n);
    Complex alpha = {1.0, 0.0};
    Complex beta = {0.0, 0.0};
    cblas_zgemm(
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states*num_states, n, n, &alpha, 
        kyrlov_basis_flatten.data(), n,
        eigenvecs_H_complex.data(), n, &beta,
        result.eigenvectors.data(), n 
    );

    sortEigenPairs(result.eigenvalues, result.eigenvectors, num_states*num_states, n, descending);
    return result;

}

double LinbladianDiagonalizer::computeInnerProduct(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M, double tol)
{
    std::vector<Complex> matC(M*M);

    Complex alpha = {1.0,0.0};
    Complex beta = {0.0,0.0};

    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
            M, M, M, &alpha, matA.data(),
            M, matB.data(), M, &beta, matC.data(), M);
    
    Complex trace = computeTrace(matC, M);
    Complex ctrace = cleanMatrixElements(trace);

    if (std::abs(ctrace.imag()) < tol){
        return ctrace.real();
    }
    else {
        // std::cout << "\nimaginary part of norm causing the problem:"<< "\n";
        // std::cout << std::abs(ctrace.imag()) << "\n";
        throw std::runtime_error("Norm is not real, Matrix is not Hermitian");
    };
}

// while doing inplace modification matrix need to be passed by reference.
void LinbladianDiagonalizer::normalizeMatrix(std::vector<Complex>& mat, int M){

    // set norm if norm argument is not specified
    double n = computeInnerProduct(mat, mat, M); 

    //inplace modification of matrix elements
    for (int i = 0; i < M*M; i++){
        mat[i] /= std::sqrt(n); 
    };

}

void LinbladianDiagonalizer::sortEigenPairs(std::vector<Complex>& eigenvalues,
                    std::vector<Complex>& eigenvectors,int row, int col,
                    bool descending) {
    // Create index vector [0, 1, 2, ..., n-1]
    std::vector<int> indices(eigenvalues.size());
    for (int i = 0; i < indices.size(); ++i)
        indices[i] = i;

    // Sort indices based on corresponding eigenvalue's real part
    std::sort(indices.begin(), indices.end(),
              [descending, &eigenvalues](int i, int j) {
                  return (descending)
                      ? (std::abs(eigenvalues[i]) > std::abs(eigenvalues[j]))
                      : (std::abs(eigenvalues[i]) < std::abs(eigenvalues[j]));
        });

    // Apply sorted order

    std::vector<Complex> sortedVals(eigenvalues.size());
    std::vector<Complex> sortedVecs(row * col);

    for (int new_j = 0; new_j < col; ++new_j) {
        int old_j = indices[new_j];
        sortedVals[new_j] = eigenvalues[old_j];

        // Copy column `old_j` → `new_j`
        for (int r = 0; r < row; ++r) {
            sortedVecs[r * col + new_j] = eigenvectors[r * col + old_j];
        }
    }

    eigenvalues = std::move(sortedVals);
    eigenvectors = std::move(sortedVecs);
}


std::vector<Complex> LinbladianDiagonalizer::applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M)
{

    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};

    std::vector<Complex> AB(M*M);
    std::vector<Complex> BA(M*M);
    // std::vector<Complex> C(num_states*num_states);

    // AB = matA*matB (H * rho)
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matA.data(),
        M, matB.data(), M, &beta, AB.data(), M
    );

    // BA = matB*matA (rho * H) 
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matB.data(),
        M, matA.data(), M, &beta, BA.data(), M
    );

    // Commutator = AB - BA
    for (int i = 0; i < M*M; i++){
        AB[i] -= BA[i];
    };

    return AB;

};

std::vector<Complex> LinbladianDiagonalizer::applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M)
{

    Complex alpha = {1.0,0.0};
    Complex beta =  {0.0,0.0};

    std::vector<Complex> AB(M*M);
    std::vector<Complex> BA(M*M);
    

    // A = matA * matB
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matA.data(),
        M, matB.data(), M, &beta, AB.data(), M
    );

    // A = matB * matA
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, M, M, &alpha, matB.data(),
        M, matA.data(), M, &beta, BA.data(), M
    );

    // Commutator = AB + BA
    for (int i = 0; i < M*M; i++){
        AB[i] += BA[i];
    };

    return AB;
};


std::vector<Complex> LinbladianDiagonalizer::constructDissipator(const std::vector<Complex>& mat,int N, 
    int num_states, double rate, DecayType decay, Scope scope)
    {
    // D(ρ)=∑_k​(J_k ​ρ J_k^†​ − 1/2{J_k^†​ J_k​,ρ})

    // TODO need to optimise Dissipator construction remove additional matrices

    std::vector<Complex> Dissipator(num_states*num_states, Complex(0,0));
    Complex beta =  {0.0,0.0};

    for (int k = 1; k <= N; k++){

        std::vector<Complex> Jk = constructJumpOperator(N,num_states,k,rate,decay,scope);

        //compute Jk^†​Jk
        Complex alpha = {0.5,0.0}; // 1/2 infront of anti commutator

        std::vector<Complex> K(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    Jk.data(), num_states, 
                    &beta, K.data(), num_states
                );

        // compute  1/2{Jk^†​ Jk​,ρ})
        std::vector<Complex> SP = applyAntiCommutator(K,mat,num_states);

        //compute L_k ​ρ
        alpha = {1.0,0.0};
        std::vector<Complex> tempA(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    num_states, num_states, num_states,
                    &alpha, Jk.data(), num_states,
                    mat.data(), num_states, 
                    &beta, tempA.data(), num_states
                );

        //compute L_k ​ρ L_k^†
        std::vector<Complex> FP(num_states*num_states, Complex(0,0));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
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
};

#ifdef BUILD_MAIN
#include "density_matrix/density_matrix.hpp"
#include "hamiltonian/hamiltonian.hpp"
#include "jump_operator/jump_operator.hpp"
#include "utils.hpp"

int main(){

    int N = 2;
    int num_states = 1 << N;
    int K = 15;
    double J = -1.0;
    double h = 0.1;

    std::vector<Complex> rho = constructRandomRho(num_states);
    
    // std::vector<Complex> ghz = {
    //     Complex (1.0 / std::sqrt(2.0)),
    //     Complex (0.0),
    //     Complex (0.0), 
    //     Complex (0.0), 
    //     Complex (0.0), 
    //     Complex (0.0), 
    //     Complex (0.0),
    //     Complex (1.0 / std::sqrt(2.0)) 
    // };

    //std::vector<Complex> rho = constructRhoFromStatevector(ghz);

    // std::vector<Complex> bell = {
    //     Complex (1.0 / std::sqrt(2.0)),
    //     Complex (0.0),
    //     Complex (0.0),
    //     Complex (1.0 / std::sqrt(2.0)) 
    // };

    // std::vector<Complex> rho = constructRhoFromStatevector(bell);

    // std::cout<<"Printing Random Matrix" << "\n";
    // printMatrix(rho, num_states, num_states);
    // std::cout<<"\n";

    // bool isHermitian = checkHermicity(rho, num_states);
    // printResult(isHermitian, "Hermicity of rho");

    // std::cout<<"Printing Hamiltonian" << "\n";
    std::vector<Complex> hamiltonian = constructHamiltonian(N,num_states,J,h);
    // printMatrix(hamiltonian, num_states, num_states);
    // std::cout<<"\n";

    LinbladianDiagonalizer LSolver(hamiltonian, N);
    // std::cout<<"Printing dRho" << "\n";
    // std::vector<Complex> rho_new = LSolver.applyLinbladian(rho);
    // printMatrix(rho_new, num_states, num_states);
    // std::cout<<"\n";
    // isHermitian = checkHermicity(rho_new, num_states);
    // printResult(isHermitian, "Hermicity of d_rho");

    // std::cout<<"Printing Hessenberg Matrix" << "\n";    std::vector<double> H = LSolver.constructHessenbergMatrix(rho, K);
    // printMatrix(H, K+1, K+1);
    // std::cout<<"\n";

    EigenResult result = LSolver.diagonalize(rho, K);
    std::cout << "\nRitzvectors\n" << "\n";
    int size = static_cast<int> (result.eigenvalues.size());
    // std::cout<< size<< "\n";
    //printMatrix(result.eigenvectors,num_states*num_states,size);
    //std::cout << "Steady State\n" << "\n";
    for (int i = 0; i <  size ; i++){
        std::cout << "eigenvector :" << i+1 << "\n";
        for (int j = 0; j <  num_states*num_states; j++){
            std::cout<<result.eigenvectors[j*(size)+i]<< "\n"; //transposed
        }

    }

    // for (int i = 0; i <  size ; i++){
    //     std::cout << "eigenvector :" << i+1 << "\n";
    //     for (int j = 0; j <  num_states; j++){
    //         int m = j*(size)+i;
    //         for (int k =0; k < num_states; k++){
    //             std::cout<<result.eigenvectors[m*num_states+k]<< " "; //transposed
    //         }
    //         std::cout << std::endl;
            
    //     }

    // }
    //printMatrix(result.eigenvectors, num_states*num_states, size);

    std::cout << "\nRitzvalues\n" << "\n";
    for(int k = 0; k < size; k++){
        std::cout << "eigenvalue " << k+1 << ":" << result.eigenvalues[k] << "\n";
    }
}
#endif