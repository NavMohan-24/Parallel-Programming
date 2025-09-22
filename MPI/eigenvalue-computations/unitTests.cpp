#include "density_matrix.hpp"
#include "hamiltonian.hpp"
#include "jump_operator.hpp"
#include "utils.hpp"


class TestDensityMatrix{

    private:
        int N;
        std::vector<Complex> psi_matrix;
        std::vector<Complex> matrix;
        int num_states;
        double tol;
    
    public:
        TestDensityMatrix(int N_, std::vector<Complex> psi_) : N(N_),tol(1e-12){
            num_states = 1 << N;
            matrix = constructRandomRho(num_states);
            psi_matrix = constructRhoFromStatevector(psi_);
        } 

        void testTraceofRandomMatrix(){
            Complex trace = computeTrace(matrix, num_states);
            auto val = std::abs(trace.real() - 1.0) < tol;
            auto val2 = std::abs(trace.imag()) < tol;

            // std::cout << val << '\t' << val2 << std::endl;
            if (std::abs(trace.real() - 1.0) < tol && std::abs(trace.imag()) < tol){
                std::cout << "Trace of the random Density Matrix constructed is correct...✅ " << '\n';
            }
            else{
                 std::cout << "Trace of the random Density Matrix constructed is NOT correct...❌" << '\n';
            }
        }
        
        void testHermicityofRandomMatrix(){

            bool isHermitian = checkHermicity(matrix,num_states);

            if (isHermitian){
                std::cout << "Random Density Matrix is Hermitian ✅" << "\n" << std::endl;
            }
            else{
                std::cout << "Random Density Matrix is NOT Hermitian ❌" << "\n" << std::endl;
            }

        }

        void testPositiveSemiDefinitivityofRandomMatrix(){
           
            // bool isPositiveSemiDefinite = checkPositiveSemiDefinitivity(psi_matrix,num_states);
            bool isPositiveSemiDefinite = true;

            auto es = findEigenvalues(matrix, num_states);

            for (int i = 0; i < es.size(); i++){
                if (es[i] < -tol) {isPositiveSemiDefinite= false;} // allow small negative tol due to FP
            }

            if (isPositiveSemiDefinite){
                std::cout << "Random Density Matrix is Positive Semi-Definite ✅" << "\n" << std::endl;
            }
            else{
                std::cout << "Random Density Matrix is NOT Positive Semi-Definite ❌" << "\n" << std::endl;
            }
            

        }

        void testTraceofPsiMatrix(){
            Complex trace = computeTrace(psi_matrix, num_states);
            if (std::abs(trace.real() - 1.0) < tol && std::abs(trace.imag()) < tol){
                std::cout << "Trace of the Density Matrix from statevector constructed is correct...✅ " << '\n';
            }
            else{
                std::cout << "Trace of the Density Matrix from statevector constructed is NOT correct...❌" << '\n';
            }
        }
        
        void testHermicityofPsiMatrix(){

            bool isHermitian = checkHermicity(psi_matrix,num_states);
            if (isHermitian){
                std::cout << "Density Matrix from statevector is Hermitian ✅" << "\n" << std::endl;
            }
            else{
                std::cout << "Density Matrix from statevector is NOT Hermitian ❌" << "\n" << std::endl;
            }
        }

        void testPositiveSemiDefinitivityofPsiMatrix(){
            // bool isPositiveSemiDefinite = checkPositiveSemiDefinitivity(psi_matrix,num_states);
            bool isPositiveSemiDefinite = true;
            auto es = findEigenvalues(psi_matrix, num_states);
            for (int i = 0; i < es.size(); i++){
                if (es[i] < -tol) {isPositiveSemiDefinite= false;} // allow small negative tol due to FP
            }
            if (isPositiveSemiDefinite){
                std::cout << "Density Matrix from statevector is Positive Semi-Definite ✅" << "\n" << std::endl;
            }
            else{
                std::cout << "Density Matrix from statevector is NOT Positive Semi-Definite ❌" << "\n" << std::endl;
            }
        }
};

class TestHamiltonian{

    private:
        const std::vector<Complex>& hamiltonian;
        int num_states;

    public:

        TestHamiltonian(const std::vector<Complex>& hamiltonian_, int num_states_) : hamiltonian(hamiltonian_), num_states(num_states_){}

        void testHermicity(){
            bool isHermitian = checkHermicity(hamiltonian, num_states);
            if (isHermitian){
                std::cout << "Hamiltonian is Hermitian ✅" << "\n" << std::endl;
            }
            else{
                std::cout << "Hamitonian is NOT Hermitian ❌" << "\n" << std::endl;
            }
        }

        void testSumofEigenvalues(){        
            Eigen::VectorXd es = findEigenvalues(hamiltonian, num_states);
            Complex trace = computeTrace(hamiltonian,num_states);
            Complex sum_eigs = es.sum();

            if (std::abs(sum_eigs.real()) < 1e-12) sum_eigs.real(0.0);
            if (std::abs(sum_eigs.imag()) < 1e-12) sum_eigs.imag(0.0);
            if (trace==sum_eigs){
                std::cout << "Sum of Eigenvalues = Sum of Trace ✅" << std::endl;
            }
            else{
                std::cout << "Sum of Eigenvalues = Sum of Trace ❌" << std::endl;
            }

        }

        void test1QHamiltonian(const std::vector<Complex>& ref){

            //std::vector<Complex> ref = {Complex(0,0), Complex(1,0), Complex(0,0), Complex(1,0)};
            std::vector<Complex> ham = constructHamiltonian(1,2,1.0,1.0);
            if (vectorEqual(ham,ref)){
                std::cout << "Single Qubit TFIM Hamiltonian is Correct ✅" << std::endl;
            }
            else{
                std::cout << "Single Qubit TFIM Hamiltonian is NOT Correct ❌" << std::endl;
            }

        }

        void test2QHamiltonian(const std::vector<Complex>& ref){

            // std::vector<Complex> ref = {Complex(1,0), Complex(1,0), Complex(1,0),Complex(0,0),
            //                             Complex(1,0), Complex(-1,0), Complex(0,0), Complex(1,0),
            //                             Complex(1,0), Complex(0,0), Complex(-1,0), Complex(1,0),
            //                             Complex(0,0), Complex(1,0), Complex(1,0), Complex(1,1)
            //                         };
            
            std::vector<Complex> ham = constructHamiltonian(2,4,1.0,1.0);
            if (vectorEqual(ham,ref)){
                std::cout << "Two Qubit TFIM Hamiltonian is Correct ✅" << std::endl;
            }
            else{
                std::cout << "Two Qubit TFIM Hamiltonian is NOT Correct ❌" << std::endl;
            }
             
        }

};

int main(){

    // unit tests for 
    int N = 3;
    int num_states = 1 << N;

    std::vector<Complex> psi = {
        Complex (1.0 / std::sqrt(2.0)),
        Complex (0.0),
        Complex (0.0), 
        Complex (0.0), 
        Complex (0.0), 
        Complex (0.0), 
        Complex (0.0),
        Complex (1.0 / std::sqrt(2.0)) 
    };

    TestDensityMatrix DensitMatrixTester(N, psi);
    DensitMatrixTester.testTraceofRandomMatrix();
    DensitMatrixTester.testHermicityofRandomMatrix();
    DensitMatrixTester.testPositiveSemiDefinitivityofRandomMatrix();
    DensitMatrixTester.testTraceofPsiMatrix();
    DensitMatrixTester.testHermicityofPsiMatrix();
    DensitMatrixTester.testPositiveSemiDefinitivityofPsiMatrix();

    // Hamiltonian Testing 

    double J = 1.0;
    double h = 1.0;

    std::vector<Complex> test_hamiltonian = constructHamiltonian(N, num_states, J, h);
    // Eigen::VectorXd evs = findEigenvalues(test_hamiltonian,num_states);
    // std::cout << evs.transpose()<<std::endl;
    TestHamiltonian HamiltonianTester(test_hamiltonian, num_states);
    HamiltonianTester.testHermicity();
    HamiltonianTester.testSumofEigenvalues();

    //TFIM unit tests
    std::vector<Complex> ref1q = {Complex(0,0), Complex(-1,0), Complex(-1,0), Complex(0,0)};
    // std::vector<Complex> ham1q = constructHamiltonian(1,2,1.0,1.0);
    std::vector<Complex> ref2q = {Complex(-1,0), Complex(-1,0), Complex(-1,0),Complex(0,0),
                                Complex(-1,0), Complex(1,0), Complex(0,0), Complex(-1,0),
                                Complex(-1,0), Complex(0,0), Complex(1,0), Complex(-1,0),
                                Complex(0,0), Complex(-1,0), Complex(-1,0), Complex(-1,0)
                                };

    // std::vector<Complex> ham2q = constructHamiltonian(2,4,1.0,1.0);
    // printMatrix(ref2q, 4);
    // printMatrix(ham2q, 4);

    HamiltonianTester.test1QHamiltonian(ref1q);
    HamiltonianTester.test2QHamiltonian(ref2q);
    
}
