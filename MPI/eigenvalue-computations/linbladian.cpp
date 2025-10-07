#include "linbladian.hpp"


inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    double re = (std::abs(z.real()) < tol) ? 0: z.real();
    double img = (std::abs(z.imag()) < tol) ? 0: z.imag();
    return {re, img};
}

inline Complex computeTrace(const std::vector<Complex>& mat, int M){

    // assumes a square matrix
    Complex trace(0.0, 0.0);
    for (int i = 0; i < M; i++){
        trace += mat[i*M + i];
    }
    return trace;
    
}

//_________________________________________________________________________

LinbladianSolver::LinbladianSolver(const std::vector<Complex>& hamiltonian_,
                int N_, double rate_, DecayType decay_, Scope scope_)
                : hamiltonian(hamiltonian_), N(N_), num_states(1<< N_),
                rate(rate_), decay(decay_), scope(scope_){}
                // validate the hamiltonian in the constructor


std::vector<Complex> LinbladianSolver::applyLinbladian(const std::vector<Complex>& rho)
{
    // validate the rho

    std::vector<Complex> comm = applyCommutator(hamiltonian,rho, num_states);
    std::vector<Complex> diss = constructDissipator(rho, N, num_states, rate, decay, scope);
    std::vector<Complex> linbladian(num_states*num_states);

    Complex i_unit(0.0, 1.0);

    for (int i = 0; i < num_states*num_states; i++){
        linbladian[i] =  cleanMatrixElements(i_unit*comm[i]- diss[i]);
    };

    return linbladian;
};        



std::vector<double> LinbladianSolver::constructHessenbergMatrix(const std::vector<Complex>& initial_rho, int k, double tol)
{   
    //lambda function to calculate index
    auto idx = [k](int i, int j){return i*(k+1)+j;};

    // Hessenberg Matrix Initialization. 
    // Elements are expected to be double since it is constructed from inner product between two Hermitian matrices
    std::vector<double> H((k+1)*(k+1),0.0); 

    // container to store kyrlov basis vectors
    // std::vector<std::vector<Complex>> basis(k+2);
    kyrlov_basis.resize(k+2);
    
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
            // compute Hessenberg matrix elements
            double nf = computeInnerProduct(kyrlov_basis[j+1], kyrlov_basis[j+1], num_states);

            if (nf < 1e-12){
                throw std::runtime_error("Arnoldi Method Breakdown : Norm became zero");
            }
            H[idx(j+1,j)] = std::sqrt(nf);
            
            //normalize basis[j+1]
            //normalizeMatrix(kyrlov_basis[j+1],num_states);
            for (int i = 0; i < num_states*num_states; i++){
                kyrlov_basis[j+1][i] /= std::sqrt(nf); 
            };
        }             
    }
    return H;
}

double LinbladianSolver::computeInnerProduct(const std::vector<Complex>& matA, const std::vector<Complex>& matB, int M, double tol)
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
        throw std::runtime_error("Norm is not real, Matrix is not Hermitian");
    };
}

// while doing inplace modification matrix need to be passed by reference.
void LinbladianSolver::normalizeMatrix(std::vector<Complex>& mat, int M){

    // set norm if norm argument is not specified
    double n = computeInnerProduct(mat, mat, M); 

    //inplace modification of matrix elements
    for (int i = 0; i < M*M; i++){
        mat[i] /= std::sqrt(n); // * to extract double from std::optional<double>
    };

}

std::vector<Complex> LinbladianSolver::applyCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M)
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

std::vector<Complex> LinbladianSolver::applyAntiCommutator(const std::vector<Complex>& matA,const std::vector<Complex>& matB, int M)
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

std::vector<Complex> LinbladianSolver::constructDissipator(const std::vector<Complex>& mat,int N, 
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
        cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
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
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
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

    int N = 3;
    int num_states = 1 << N;
    int K = 3;
    double J = -1.0;
    double h = 2.0;

    std::vector<Complex> rho = constructRandomRho(num_states);
    std::cout<<"Printing Random Matrix" << "\n";
    printMatrix(rho, num_states);
    std::cout<<"\n";

    bool isHermitian = checkHermicity(rho, num_states);
    printResult(isHermitian, "Hermicity of rho");

    std::cout<<"Printing Hamiltonian" << "\n";
    std::vector<Complex> hamiltonian = constructHamiltonian(N,num_states,J,h);
    printMatrix(hamiltonian, num_states);
    std::cout<<"\n";

    LinbladianSolver LSolver(hamiltonian, N);
    std::cout<<"Printing dRho" << "\n";
    std::vector<Complex> rho_new = LSolver.applyLinbladian(rho);
    printMatrix(rho_new, num_states);
    std::cout<<"\n";
    isHermitian = checkHermicity(rho_new, num_states);
    printResult(isHermitian, "Hermicity of d_rho");

    std::cout<<"Printing Hessenberg Matrix" << "\n";
    std::vector<double> H = LSolver.constructHessenbergMatrix(rho, K);
    printMatrix(H, K+1);
    std::cout<<"\n";


}
#endif