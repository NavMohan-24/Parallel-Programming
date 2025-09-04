#include "random_rho.hpp"



Complex genMatrixElements(std::mt19937 &gen, std::normal_distribution<double> &dist){
    return {dist(gen), dist(gen)};
}
inline Complex cleanMatrixElements(const Complex& z, double tol=1e-12){
    // need to define it 
    // why we need const?
    
}


Vector createRandomRho(int N){

    int num_states = 1 << N;

    Vector A(num_states*num_states);
    Vector rho(num_states*num_states); // declares a vector with all elements zero.

    std::random_device rd;
    std::mt19937 gen(rd);
    std::normal_distribution<double> dist(0.0,1.0);

    for (int i = 0; i < num_states ; i++){
        for (int j = 0; j < num_states ; j++){
            A[i*num_states+j] = genMatrixElements(gen,dist);
        }
    }
    double alpha = 1.0, beta = 0.0;

    // multiply AA^{H} 
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasConjTrans,
        num_states,num_states,num_states,
        alpha, A.data(), num_states, A.data(), num_states,
        beta, rho.data(), num_states); 
    
    
    double trace = 0.0;
    for (int i = 0; i < num_states; i++){
            trace += std::real(i*num_states+i); // diagonal values are guarnteed to be real for hermitian.
        
    }

    



    


    





    return rho
}


#ifdef BUILD_MAIN
int main() {
    std::complex<double> z(3.0, 4.0);   // 3 + 4i
    auto z_conj = std::conj(z);         // 3 - 4i

    std::cout << "z = " << z << "\n";
    std::cout << "conj(z) = " << z_conj << "\n";
}
#endif