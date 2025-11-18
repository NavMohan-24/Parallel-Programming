#include "linbladian_diagonalizer.hpp"

//=======================================================================================
// Helper Functions
//=======================================================================================

/**
 * @brief Cleans matrix elements by setting values below tolerance to zero.
 * 
 * @param z Complex number to clean
 * @param tol Tolerance threshold for considering values as zero
 * @return Complex number with small components zeroed out
 */
inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    double re = (std::abs(z.real()) < tol) ? 0: z.real();
    double img = (std::abs(z.imag()) < tol) ? 0: z.imag();
    return {re, img};
}

/**
 * @brief Computes the trace of a square matrix.
 * 
 * @param mat Matrix in flattened row-major format
 * @param M Dimension of the square matrix (M x M)
 * @return Complex trace value
 */
inline Complex computeTrace(const std::vector<Complex>& mat, int M){
    Complex trace(0.0, 0.0);
    for (int i = 0; i < M; i++){
        trace += mat[i*M + i];
    }
    return trace;
}

void EigenResult::sortEigenPairs(int row, int col, bool descending) {
    // Create index vector for indirect sorting
    std::vector<int> indices(eigenvalues.size());
    for (size_t i = 0; i < indices.size(); ++i)
        indices[i] = i;

    // Sort indices by magnitude of corresponding eigenvalues
    std::sort(indices.begin(), indices.end(),
              [descending, this](int i, int j) {
                  return (descending)
                      ? (std::abs(eigenvalues[i].real()) > std::abs(eigenvalues[j].real()))
                      : (std::abs(eigenvalues[i].real()) < std::abs(eigenvalues[j].real()));
        });

    // Apply sorted order to eigenvalues and eigenvectors
    std::vector<Complex> sortedVals(eigenvalues.size());
    std::vector<Complex> sortedVecs(row * col);

    for (int new_j = 0; new_j < col; ++new_j) {
        int old_j = indices[new_j];
        sortedVals[new_j] = eigenvalues[old_j];

        // Copy column old_j to column new_j
        for (int r = 0; r < row; ++r) {
            sortedVecs[r * col + new_j] = eigenvectors[r * col + old_j];
        }
    }

    eigenvalues = std::move(sortedVals);
    eigenvectors = std::move(sortedVecs);
}

//=======================================================================================
// Lindbladian Diagonalizer Implementation
//=======================================================================================

ArnoldiLinbladianDiagonalizer::ArnoldiLinbladianDiagonalizer(LinbladianConstructor& linbladianConstructor_, int N_, int d_)
                : linbladianConstructor(linbladianConstructor_), N(N_), d(d_), num_states(std::pow(d_,N_)){}

EigenResult ArnoldiLinbladianDiagonalizer::diagonalize(const std::vector<Complex>& init_rho, int k, bool descending, double tol){

    // Validate requested Krylov subspace dimension
    int max_k = num_states*num_states - 1;
    if (k > max_k){
        throw std::runtime_error("Number of requested eigen vals/vecs must be <= N.");
    }

    // Construct Hessenberg matrix via Arnoldi iteration
    std::vector<double> H = constructHessenbergMatrix(init_rho, k, tol);
    EigenResult result;
    
    int size = static_cast<int>(H.size()); 
    int n = std::sqrt(size); // Dimension of Hessenberg matrix (actual size may be < k+1 if early termination occurred)

    // LAPACK parameters for submatrix specification (using full matrix here)
    int ilo = 1, ihi = n;

    // Containers for eigenvalue decomposition
    std::vector<double> wr(n), wi(n); // Real and imaginary parts of eigenvalues
    std::vector<double> Z(n*n, 0.0);  // Schur vectors
    
    // Initialize Z as identity matrix
    for (int i = 0; i < n; i++) {
        Z[i*n + i] = 1.0;
    }

    // Convert H to Schur form and compute Schur vectors using QR algorithm
    int info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR,
        'S', 'V',  // 'S' = compute Schur form, 'V' = compute Schur vectors
        n, ilo, ihi,
        H.data(), n, 
        wr.data(), wi.data(),
        Z.data(), n
    ); 
    
    if (info != 0) {
        throw std::runtime_error("DHSEQR failed in diagonalizing the matrix");
    }

    // Store eigenvalues (eigenvalues of Hessenberg matrix approximate those of Lindbladian)
    result.eigenvalues.resize(n);
    for (int i=0; i < n; i++){
        result.eigenvalues[i] = std::complex<double>(wr[i], wi[i]);
    }
    
    // Compute eigenvectors of Schur form
    std::vector<double> VR(n*n); // Right eigenvectors of Schur form
    int mm = n; // Leading dimension for VR
    int m;      // Number of eigenvectors found (output parameter)
    
    info = LAPACKE_dtrevc(
        LAPACK_ROW_MAJOR,
        'R', 'A',  // 'R' = right eigenvectors, 'A' = compute all eigenvectors
        nullptr, n,
        H.data(), n,
        nullptr, n,
        VR.data(), n,
        mm, &m
    );

    if (info != 0) {
        throw std::runtime_error("DTREVC failed in diagonalizing the matrix");
    }

    // Transform eigenvectors back: eigenvectors of H = Z * VR
    std::vector<double> eigenvecs_H(n*n);
    cblas_dgemm(
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        n, n, n, 1.0, 
        Z.data(), n,
        VR.data(), n, 
        0.0, eigenvecs_H.data(), n 
    );

    // Convert to complex for final transformation
    std::vector<Complex> eigenvecs_H_complex(n*n);
    for (int i = 0; i < n*n; i++) {
        eigenvecs_H_complex[i] = Complex(eigenvecs_H[i], 0.0);
    }

    // Transform eigenvectors from Krylov subspace back to full density matrix basis
    // result_eigenvectors = kyrlov_basis_flatten * eigenvecs_H_complex
    result.eigenvectors.resize(num_states*num_states*n);
    Complex alpha = {1.0, 0.0};
    Complex beta = {0.0, 0.0};
    cblas_zgemm(
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        num_states*num_states, n, n, &alpha, 
        kyrlov_basis_flatten.data(), n,
        eigenvecs_H_complex.data(), n, 
        &beta, result.eigenvectors.data(), n 
    );

    // Sort eigenvalue-eigenvector pairs
    result.sortEigenPairs(num_states*num_states, n , descending);
    return result;

}


std::vector<double> ArnoldiLinbladianDiagonalizer::constructHessenbergMatrix(const std::vector<Complex>& initial_rho, int k, double tol)
{     
    // Lambda function to calculate matrix index in row-major order
    auto idx = [k](int i, int j){return i*(k+1)+j;};

    // Initialize Hessenberg matrix (elements are real since constructed from 
    // inner products of Hermitian matrices)
    std::vector<double> H((k+1)*(k+1), 0.0); 

    // Container to store Krylov basis vectors (size k+2 to accommodate j+1 indexing)
    std::vector<std::vector<Complex>> kyrlov_basis(k+2);

    double nf; // Norm variable
  
    // Normalize the initial density matrix to form first Krylov basis vector
    kyrlov_basis[0] = initial_rho;
    nf = computeInnerProduct(kyrlov_basis[0], kyrlov_basis[0], num_states);
    
    if (nf < tol){
        throw std::runtime_error("Input density matrix has zero norm. Arnoldi iteration could not be performed.");
    }
    
    for (int i = 0; i < num_states*num_states; i++){
        kyrlov_basis[0][i] = kyrlov_basis[0][i]/std::sqrt(nf); 
    }
        
    // Arnoldi iteration to build Krylov subspace and Hessenberg matrix
    for (int j = 0; j <= k; j++){

        // Apply Lindbladian superoperator to generate next Krylov vector
        kyrlov_basis[j+1] = linbladianConstructor.applyLinbladian(kyrlov_basis[j]);
        
        // Orthogonalize against all previous basis vectors (Modified Gram-Schmidt)
        for (int i = 0; i <= j; i++)
        {   
            // Compute Hessenberg matrix element: H[i,j] = <basis[i], L(basis[j])>
            H[idx(i,j)] = computeInnerProduct(kyrlov_basis[i], kyrlov_basis[j+1], num_states);

            // Subtract projection onto basis[i]
            for (int m = 0; m < num_states*num_states; m++){
                kyrlov_basis[j+1][m] -= H[idx(i,j)]*kyrlov_basis[i][m];
            }
        }
        
        if (j < k){
            // Compute norm of orthogonalized vector
            nf = computeInnerProduct(kyrlov_basis[j+1], kyrlov_basis[j+1], num_states);

            // Early termination if Krylov subspace is complete (basis vector has zero norm)
            if (nf < 1e-12){
                
                // Resize Hessenberg matrix to actual dimension reached
                std::vector<double> H_new((j+1)*(j+1));

                // Copy elements from full-sized H to resized H_new
                for (int m=0; m < j+1; m++){
                    for (int n=0; n < j+1; n++){
                        H_new[m*(j+1)+n] = H[m*(k+1)+n];
                    }
                }
                
                // Flatten and resize Krylov basis vectors to actual dimension
                kyrlov_basis_flatten.resize(num_states*num_states*(j+1));
                for (int i = 0; i < j+1; i++){
                    for (int m = 0; m < num_states*num_states; m++){
                        kyrlov_basis_flatten[m*(j+1)+ i] = kyrlov_basis[i][m]; 
                    }
                }
                return H_new;
            }
            
            // Store norm as subdiagonal element of Hessenberg matrix
            H[idx(j+1,j)] = std::sqrt(nf);
        
            // Normalize the new Krylov basis vector
            for (int i = 0; i < num_states*num_states; i++){
                kyrlov_basis[j+1][i] = kyrlov_basis[j+1][i]/std::sqrt(nf); 
            }
        }          
    }

    // Flatten Krylov basis vectors into column-major format for matrix operations
    kyrlov_basis_flatten.resize(num_states*num_states*(k+1));
    for (int i = 0; i < k+1; i++){
        for (int j = 0; j < num_states*num_states; j++){
            kyrlov_basis_flatten[j*(k+1) + i] = kyrlov_basis[i][j];
        }
    }

    return H;
}

double ArnoldiLinbladianDiagonalizer::computeInnerProduct(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M, double tol)
{
    // Compute Frobenius inner product: Tr(A† * B)
    std::vector<Complex> matC(M*M);

    Complex alpha = {1.0, 0.0};
    Complex beta = {0.0, 0.0};

    // Compute matC = A† * B
    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
            M, M, M, &alpha, 
            matA.data(), M, 
            matB.data(), M, 
            &beta, matC.data(), M);
    
    // Compute trace of result
    Complex trace = computeTrace(matC, M);
    Complex ctrace = cleanMatrixElements(trace);

    // Inner product must be real for Hermitian matrices
    if (std::abs(ctrace.imag()) < tol){
        return ctrace.real();
    }
    else {
        throw std::runtime_error("Norm is not real, Matrix is not Hermitian");
    }
}
