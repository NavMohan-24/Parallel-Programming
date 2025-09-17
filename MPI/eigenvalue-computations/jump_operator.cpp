#include "jump_op.hpp"
#include <boost/dynamic_bitset.hpp>

std::vector<Complex> constructJumpOperator(int N, int num_states, int k, double rate = 1.0, DecayType type = DecayType::Damping, Scope scope = Scope::Local){

    std::vector<Complex> jump_k(num_states*num_states, Complex(0,0)); 

    double r = std::pow(rate,0.5);

    switch(scope)
    {
    case Scope::Local:
        if (type == DecayType::Damping){
            for (int m = 0; m < num_states; m++){
                int n = m + (1 << (k-1));
                if (n < num_states){
                    boost::dynamic_bitset<> bn(N, n);
                    jump_k[m*num_states+n] = Complex((bn.test(k-1)==true) ? r*1:0, 0);
                };
            };
            return jump_k;
        }
        if (type == DecayType::Pumping){
            for (int m = 0; m < num_states; m++){
                int n = m - (1 << (k-1));
                if (n >= 0){
                    boost::dynamic_bitset<> bn(N, n);
                    jump_k[m*num_states+n] = Complex((bn.test(k-1)==false) ? r*1:0, 0);
                };
            };
            return jump_k;
        }
        if (type == DecayType::Dephasing){
            for (int m =0; m < num_states; m++){
                boost::dynamic_bitset<> bn(N, m);
                if ((k-1) < N){
                    jump_k[m*num_states+m] = Complex((bn.test(k-1)==true) ? -r*1:r*1, 0);
                }
            };
            return jump_k;
        }
        break;

    case Scope::Collective:
        //TODO: Collective Decay's
        return jump_k;
    };
   return jump_k;
};

#ifdef BUILD_MAIN
#include <iomanip>
#include "random_rho.hpp"

bool complexEqual(const Complex& x, const Complex& y, double eps = 1e-12) {
    return std::abs(x - y) < eps;
    };

bool vectorEqual(const std::vector<Complex>& a,
                const std::vector<Complex>& b,
                double eps = 1e-12) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (!complexEqual(a[i], b[i], eps)) return false;
    }
    return true;
};

std::vector<Complex> matrixMultiplication(std::vector<Complex>& A,std::vector<Complex>&  B, int N){

    std::vector<Complex> mat(N*N);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){

                mat[i*N+j] += A[i*N+k]*B[k*N+j];
                 
            }
        }
    }
    return mat;
};

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

std::vector<Complex> hermitian(const std::vector<Complex>& A, int N){
    std::vector<Complex> B(N * N, Complex(0, 0));

    for (int i = 0; i < N; i++){
        for (int j=0; j < N; j++){
            B[j*N+i] = std::conj(A[i*N+j]);
        }
    }
    return B;
}

void printMatrix(const std::vector<Complex>& M, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(2) 
                      << M[i * N + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

Complex trace(const std::vector<Complex>& M, int N){ 
    Complex tr(0.0,0.0);
    for (int i = 0; i < N; i++) {
        tr += M[i*N+i];
    }
    return tr;
}

int main(){

    std::vector<Complex> Identity = {Complex(1,0), Complex(0,0), Complex(0,0), Complex(1,0)};
    std::vector<Complex> SigmaMinus = {Complex(0,0), Complex(1,0), Complex(0,0), Complex(0,0)};
    std::vector<Complex> SigmaPlus = {Complex(0,0), Complex(0,0), Complex(1,0), Complex(0,0)};
    std::vector<Complex> SigmaZ = {Complex(1,0), Complex(0,0), Complex(0,0), Complex(-1,0)};

    DecayType damping = DecayType::Damping;
    DecayType pumping = DecayType::Pumping;
    DecayType dephasing = DecayType::Dephasing;

    std::cout << "Running Unit tests..🔬" << "\n" << std::endl;
    
    int N = 1;
    int num_states = 1 << N;
    int k = 1;

    // constructing jump operators for N=1 and k=1
    std::vector<Complex> mat11m = constructJumpOperator(N,num_states,k,1.0,damping);
    std::vector<Complex> mat11p = constructJumpOperator(N,num_states,k,1.0,pumping);
    std::vector<Complex> mat11z = constructJumpOperator(N,num_states,k,1.0,dephasing);
    
    // Test
    if (vectorEqual(mat11m,SigmaMinus)){
        std::cout << "Damping decay type pass the unit test for N=1 ✅" << std::endl;
    }
    else{
        std::cout << "Damping decay type fails the unit test for N=1 ❌" << std::endl;
    }

    if (vectorEqual(mat11p,SigmaPlus)){
        std::cout << "Pumping decay type pass the unit test for N=1 ✅" << std::endl;
    }
    else{
        std::cout << "Pumping decay type fails the unit test for N=1 ❌" << std::endl;
    }

    if (vectorEqual(mat11z,SigmaZ)){
        std::cout << "SigmaZ decay type pass the unit test for N=1 ✅" << std::endl;
    }
    else{
        std::cout << "SigmaZ decay type fails the unit test for N=1 ❌" << std::endl;
    }

    std::cout << "\n" << std::endl;

    N = 2;
    num_states = 1 << N;
    
    // constructing jump operators for N=1 and k=1
    std::vector<Complex> mat21m = constructJumpOperator(N,num_states,1,1.0,damping);
    std::vector<Complex> mat22m = constructJumpOperator(N,num_states,2,1.0,damping);


    std::vector<Complex> m1 = kroneckerProduct(Identity, N, SigmaMinus, N);
    std::vector<Complex> m2 = kroneckerProduct(SigmaMinus, N, Identity, N);

    if (vectorEqual(mat21m,m1) & vectorEqual(mat22m, m2)){
        std::cout << "Damping decay type pass the unit test for N=2 ✅" << std::endl;
    }
    else{
        std::cout << "Damping decay type fails the unit test for N=2 ❌" << std::endl;
    }

    std::vector<Complex> mat21p = constructJumpOperator(N,num_states,1,1.0,pumping);
    std::vector<Complex> mat22p = constructJumpOperator(N,num_states,2,1.0,pumping);

    std::vector<Complex> p1 = kroneckerProduct(Identity, N, SigmaPlus, N);
    std::vector<Complex> p2 = kroneckerProduct(SigmaPlus, N, Identity, N);

    // printMatrix(mat21p,num_states);
    // printMatrix(p1, num_states);


    if (vectorEqual(mat21p,p1) && vectorEqual(mat22p, p2)){
        std::cout << "Pumping decay type pass the unit test for N=2 ✅" << std::endl;
    }
    else{
        std::cout << "Pumping decay type fails the unit test for N=2 ❌" << std::endl;
    }

    std::vector<Complex> mat21z = constructJumpOperator(N,num_states,1,1.0,dephasing);
    std::vector<Complex> mat22z = constructJumpOperator(N,num_states,2,1.0,dephasing);

    std::vector<Complex> z1 = kroneckerProduct(Identity, N, SigmaZ, N);
    std::vector<Complex> z2 = kroneckerProduct(SigmaZ, N, Identity, N);

    if (vectorEqual(mat21z,z1) && vectorEqual(mat22z, z2)){
        std::cout << "Dephasing decay type pass the unit test for N=2 ✅" << std::endl;
    }
    else{
        std::cout << "Dephasing decay type fails the unit test for N=2 ❌" << std::endl;
    }

    // Testing Properties of Jump Operators
   
    N = 3;
    num_states = 1 << N;
    std::vector<Complex> rho = createRandomRho(num_states);
    Complex trace1(0.0,0.0);
    Complex trace2(0.0,0.0);

    std::cout << "\nTesting Cyclic Properties of Damping Jump Operators..🔬" << std::endl;

    std::vector<Complex> mat32m = constructJumpOperator(N,num_states,2,1.0,damping);
    std::vector<Complex> mat32mt = hermitian(mat32m, num_states);

    // Calculate Tr(J*rho*J_dagger)
    std::vector<Complex> temp1m = matrixMultiplication(mat32m, rho, num_states);
    std::vector<Complex> M1 = matrixMultiplication(temp1m, mat32mt, num_states);
    trace1 = trace(M1, num_states);

    // Calculate Tr(J_dagger*J*rho)
    std::vector<Complex> temp2m = matrixMultiplication(mat32mt, mat32m, num_states);
    std::vector<Complex> M2 = matrixMultiplication(temp2m, rho, num_states);
    trace2 = trace(M2, num_states);
    
    if (complexEqual(trace1, trace2)){
        std::cout << "Traces are Equal ✅"<<std::endl;
    }
    else{
        std::cout << "Traces are Not Equal ❌"<<std::endl;
        std::cout << "Trace 1: " << trace1 << std::endl;
        std::cout << "Trace 2: " << trace2 << std::endl;
    }

    std::cout << "\nTesting Cyclic Properties of Pumping Jump Operators..🔬" << std::endl;

    // Change from damping to pumping, and k from 2 to 1
    std::vector<Complex> mat31p = constructJumpOperator(N, num_states, 1, 1.0, pumping);
    std::vector<Complex> mat31pt = hermitian(mat31p, num_states); // Using mat31p for hermitian
    // std::vector<Complex> rho = createRandomRho(num_states);

    // Calculate Tr(J*rho*J_dagger)
    std::vector<Complex> temp1p = matrixMultiplication(mat31p, rho, num_states);
    std::vector<Complex> P1 = matrixMultiplication(temp1p, mat31pt, num_states);
    trace1 = trace(P1, num_states);

    // Calculate Tr(J_dagger*J*rho)
    std::vector<Complex> temp2p = matrixMultiplication(mat31pt, mat31p, num_states);
    std::vector<Complex> P2 = matrixMultiplication(temp2p, rho, num_states);
    trace2 = trace(P2, num_states);
        
    if (complexEqual(trace1, trace2)){
        std::cout << "Traces are Equal ✅"<<std::endl;
    }
    else{
        std::cout << "Traces are Not Equal ❌"<<std::endl;
        std::cout << "Trace 1: " << trace1 << std::endl;
        std::cout << "Trace 2: " << trace2 << std::endl;
    }

    std::cout << "\nTesting Cyclic Properties of Dephasing Jump Operators..🔬"  << std::endl;
    

    // Change from pumping to dephasing, and k from 1 to 3
    std::vector<Complex> mat33z = constructJumpOperator(N, num_states, 3, 1.0, dephasing);
    std::vector<Complex> mat33zt = hermitian(mat33z, num_states); // Using mat33z for hermitian
    // std::vector<Complex> rho = createRandomRho(num_states);

    // Calculate Tr(J*rho*J_dagger)
    std::vector<Complex> temp1z = matrixMultiplication(mat33z, rho, num_states);
    std::vector<Complex> Z1 = matrixMultiplication(temp1z, mat33zt, num_states);
    trace1 = trace(Z1, num_states);

    // Calculate Tr(J_dagger*J*rho)
    std::vector<Complex> temp2z = matrixMultiplication(mat33zt, mat33z, num_states);
    std::vector<Complex> Z2 = matrixMultiplication(temp2z, rho, num_states);
    trace2 = trace(Z2, num_states);
        
    if (complexEqual(trace1, trace2)){
        std::cout << "Traces are Equal ✅"<<std::endl;
    }
    else{
        std::cout << "Traces are Not Equal ❌"<<std::endl;
        std::cout << "Trace 1: " << trace1 << std::endl;
        std::cout << "Trace 2: " << trace2 << std::endl;
    }
}
#endif

