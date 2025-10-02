#include "density_matrix.hpp"



Complex genMatrixElements(std::mt19937 &gen, std::normal_distribution<double> &dist){
    return {dist(gen), dist(gen)};
}

inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    double re = (std::abs(z.real()) < tol) ? 0: z.real();
    double img = (std::abs(z.imag()) < tol) ? 0: z.imag();
    return {re, img};
}

std::vector<Complex> constructRandomRho(int num_states){
    

    std::vector<Complex> A(num_states*num_states);
    std::vector<Complex> rho(num_states*num_states);

    #pragma omp parallel shared(A)
    {

        std::normal_distribution<double> dist(0.0,1.0);
        std::random_device rd;
        std::seed_seq seq{rd(),rd(),rd(),static_cast<unsigned>(omp_get_thread_num())}; //seeding class that takes few seed values and expand them into many high-quality seed values.
        std::mt19937 gen(seq);

        #pragma omp for collapse(2) schedule(static)
        for (int i = 0; i < num_states ; i++){
            for (int j = 0; j < num_states ; j++){ 
                A[i*num_states+j] = genMatrixElements(gen,dist);
            }
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
    #pragma omp parallel for reduction(+:trace)
    for (int i = 0; i < num_states; i++){
            trace += std::real(rho[i*num_states+i]); // diagonal values are guarnteed to be real for hermitian.
        
    }
    #pragma omp parallel for shared(rho)
    for (int i = 0; i < num_states*num_states; i++){
        rho[i] = cleanMatrixElements(rho[i]/trace);
    }


    return rho;
}


std::vector<Complex> constructRhoFromStatevector(std::vector<Complex>& psi){

    int N = (int) psi.size();

    std::vector<Complex> rho(N*N);

    #pragma omp parallel for collapse(2) schedule(static) shared(rho)
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            rho[i*N+j] = psi[i]*std::conj(psi[j]); // outer product
        }
    }

    return rho;

}

