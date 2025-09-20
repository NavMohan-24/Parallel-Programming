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
           
            bool isPositiveSemiDefinite = checkPositiveSemiDefinitivity(matrix, num_states);

            if (isPositiveSemiDefinite){
                std::cout << "Random Matrix is Positive Semi-Definite ✅" << "\n" << std::endl;
            }
            else{
                std::cout << "Random Matrix is NOT Positive Semi-Definite ❌" << "\n" << std::endl;
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
           
            bool isPositiveSemiDefinite = checkPositiveSemiDefinitivity(psi_matrix,num_states);

            if (isPositiveSemiDefinite){
                std::cout << "Density Matrix from statevector is Positive Semi-Definite ✅" << "\n" << std::endl;
            }
            else{
                std::cout << "Density Matrix from statevector is NOT Positive Semi-Definite ❌" << "\n" << std::endl;
            }
            

        }

};

int main(){

    // unit tests for 
    int N = 3;

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
    
}
