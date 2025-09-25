#include "density_matrix.hpp"
#include "hamiltonian.hpp"
#include "jump_operator.hpp"
#include "utils.hpp"


class TestDensityMatrix{  
    public:
        void testTrace(const std::vector<Complex>& matrix, int num_states, double tol = 1e-12){
            Complex trace = computeTrace(matrix,  num_states);
            bool condition = std::abs(trace.real() - 1.0) < tol && std::abs(trace.imag()) < tol;
            printResult(condition,"Trace of the Density Matrix");
        }
        
        void testHermicity(const std::vector<Complex>& matrix, int num_states){

            bool isHermitian = checkHermicity(matrix,num_states);
            printResult(isHermitian,"Hermicity of the Density Matrix");

        }

        void testPositiveSemiDefinitivity(const std::vector<Complex>& matrix, int num_states, double tol = 1e-12){
           
            bool isPositiveSemiDefinite = true;
            auto es = findEigenvalues(matrix, num_states);
            for (int i = 0; i < es.size(); i++){
                if (es[i] < -tol) {isPositiveSemiDefinite= false;} // allow small negative tol due to FP
            }
            printResult(isPositiveSemiDefinite,"Positive SemiDefinitivity of the Density Matrix");
        }

    };

class TestHamiltonian{

    public:
        void testHermicity(const std::vector<Complex>& hamiltonian, int num_states){
            bool isHermitian = checkHermicity(hamiltonian, num_states);
            printResult(isHermitian, "Hermicity of Hamiltonian" );
        }

        void testSumofEigenvalues(const std::vector<Complex>& hamiltonian, int num_states){        
            Eigen::VectorXd es = findEigenvalues(hamiltonian, num_states);
            Complex trace = computeTrace(hamiltonian,num_states);
            Complex sum_eigs = es.sum();

            if (std::abs(sum_eigs.real()) < 1e-12) sum_eigs.real(0.0);
            if (std::abs(sum_eigs.imag()) < 1e-12) sum_eigs.imag(0.0);
            printResult(trace==sum_eigs, "Sum of Eigenvalues");

        }

        void testHamiltonian(const std::vector<Complex>& hamiltonian, const std::vector<Complex>& ref){

            printResult(vectorEqual(hamiltonian,ref), "Hamiltonian Unit Test");
        }


};

class TestJumpOperator{


    public:

        void testJumpOperator(const std::vector<Complex>& mat, const std::vector<Complex>& ref){

            printResult(vectorEqual(mat,ref), "Jump Operator Unit Test");
        }

        void testCyclicProperties(const std::vector<Complex>& mat, const std::vector<Complex>& rho, int num_states, double tol = 1e-12){

            std::vector<Complex> mat_trans = hermitian(mat, num_states);

            // Calculate Tr(J*rho*J_dagger)
            std::vector<Complex> temp1m = matrixMultiplication(mat, rho, num_states, num_states, num_states);
            std::vector<Complex> M1 = matrixMultiplication(temp1m, mat_trans, num_states, num_states, num_states);
            Complex trace1 = computeTrace(M1, num_states);

            // Calculate Tr(J_dagger*J*rho)
            std::vector<Complex> temp2m = matrixMultiplication(mat_trans, mat, num_states, num_states, num_states);
            std::vector<Complex> M2 = matrixMultiplication(temp2m, rho, num_states, num_states, num_states);
            Complex trace2 = computeTrace(M2, num_states);

            bool isCyclic = std::abs(trace1 - trace2) < tol;

            printResult(isCyclic, "Cyclic Property of Jump Operator");

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

    std::vector<Complex> matrix = constructRandomRho(num_states);
    std::vector<Complex> psi_matrix = constructRhoFromStatevector(psi);
    
    std::cout << "\nRunning Tests for Random Density Matrix ..🔬" << "\n";
    TestDensityMatrix DensitMatrixTester;
    DensitMatrixTester.testTrace(matrix,num_states);
    DensitMatrixTester.testHermicity(matrix, num_states);
    DensitMatrixTester.testPositiveSemiDefinitivity(matrix, num_states);

    std::cout << "\nRunning Tests for GHZ Density Matrix ..🔬" << "\n";
    DensitMatrixTester.testTrace(psi_matrix,num_states);
    DensitMatrixTester.testHermicity(psi_matrix, num_states);
    DensitMatrixTester.testPositiveSemiDefinitivity(psi_matrix, num_states);
 

    // Hamiltonian Testing 
    double J = 1.0;
    double h = 1.0;
    std::vector<Complex> test_hamiltonian = constructHamiltonian(N, num_states, J, h);

    TestHamiltonian HamiltonianTester;
    std::cout << "\nRunning Tests for Hamilitonian ..🔬" << "\n";
    HamiltonianTester.testHermicity(test_hamiltonian, num_states);
    HamiltonianTester.testSumofEigenvalues(test_hamiltonian, num_states);

    //TFIM unit tests
    std::cout << "\nRunning Tests for Single Qubit Hamilitonian ..🔬" << "\n";
    std::vector<Complex> ref1q = {Complex(0,0), Complex(-1,0), Complex(-1,0), Complex(0,0)};
    std::vector<Complex> ham1q = constructHamiltonian(1,2,1.0,1.0);
    HamiltonianTester.testHamiltonian(ham1q, ref1q);
    std::cout << "\nRunning Tests for Two Qubit Hamilitonian ..🔬" << "\n";    
    std::vector<Complex> ref2q = {Complex(-1,0), Complex(-1,0), Complex(-1,0),Complex(0,0),
                                Complex(-1,0), Complex(1,0), Complex(0,0), Complex(-1,0),
                                Complex(-1,0), Complex(0,0), Complex(1,0), Complex(-1,0),
                                Complex(0,0), Complex(-1,0), Complex(-1,0), Complex(-1,0)
                                };
    std::vector<Complex> ham2q = constructHamiltonian(2,4,1.0,1.0);
    HamiltonianTester.testHamiltonian(ham2q, ref2q);


    //Jump Operator testing
    std::vector<Complex> Identity = {Complex(1,0), Complex(0,0), Complex(0,0), Complex(1,0)};
    std::vector<Complex> SigmaMinus = {Complex(0,0), Complex(1,0), Complex(0,0), Complex(0,0)};
    std::vector<Complex> SigmaPlus = {Complex(0,0), Complex(0,0), Complex(1,0), Complex(0,0)};
    std::vector<Complex> SigmaZ = {Complex(1,0), Complex(0,0), Complex(0,0), Complex(-1,0)};

    DecayType damping = DecayType::Damping;
    DecayType pumping = DecayType::Pumping;
    DecayType dephasing = DecayType::Dephasing;

    // constructing jump operators for N=1 and k=1
    std::vector<Complex> mat11m = constructJumpOperator(1,2,1,1.0,damping);
    std::vector<Complex> mat11p = constructJumpOperator(1,2,1,1.0,pumping);
    std::vector<Complex> mat11z = constructJumpOperator(1,2,1,1.0,dephasing);

    std::cout << "\nRunning Unit Tests for Single Qubit Jump Operators ..🔬" << "\n";
    TestJumpOperator JumpOperatorTester;
    JumpOperatorTester.testJumpOperator(mat11m,SigmaMinus);
    JumpOperatorTester.testJumpOperator(mat11p,SigmaPlus);
    JumpOperatorTester.testJumpOperator(mat11z,SigmaZ);

     // constructing jump operators for N=2 
    std::vector<Complex> mat21m = constructJumpOperator(2,4,1,1.0,damping);
    std::vector<Complex> mat22m = constructJumpOperator(2,4,2,1.0,damping);
    std::vector<Complex> m1 = kroneckerProduct(Identity, 2, SigmaMinus, 2);
    std::vector<Complex> m2 = kroneckerProduct(SigmaMinus, 2, Identity, 2);

    std::cout << "\nRunning Unit Tests for Two Qubit Damping Jump Operators ..🔬" << "\n";
    JumpOperatorTester.testJumpOperator(mat21m,m1);
    JumpOperatorTester.testJumpOperator(mat22m,m2);

    std::vector<Complex> mat21p = constructJumpOperator(2,4,1,1.0,pumping);
    std::vector<Complex> mat22p = constructJumpOperator(2,4,2,1.0,pumping);
    std::vector<Complex> p1 = kroneckerProduct(Identity, 2, SigmaPlus, 2);
    std::vector<Complex> p2 = kroneckerProduct(SigmaPlus, 2, Identity, 2);

    std::cout << "\nRunning Unit Tests for Two Qubit Pumping Jump Operators ..🔬" << "\n";
    JumpOperatorTester.testJumpOperator(mat21p,p1);
    JumpOperatorTester.testJumpOperator(mat22p,p2);

    std::vector<Complex> mat21z = constructJumpOperator(2,4,1,1.0,dephasing);
    std::vector<Complex> mat22z = constructJumpOperator(2,4,2,1.0,dephasing);
    std::vector<Complex> z1 = kroneckerProduct(Identity, 2, SigmaZ, 2);
    std::vector<Complex> z2 = kroneckerProduct(SigmaZ, 2, Identity, 2);

    std::cout << "\nRunning Unit Tests for Two Qubit Dephasing Jump Operators ..🔬" << "\n";
    JumpOperatorTester.testJumpOperator(mat21z,z1);
    JumpOperatorTester.testJumpOperator(mat22z,z2);

    std::cout << "\nTesting Cyclic Properties of Pumping Jump Operators..🔬" << std::endl;


    std::vector<Complex> mat32m = constructJumpOperator(N,num_states,2,1.0,damping);
    JumpOperatorTester.testCyclicProperties(mat32m, matrix, num_states);

    std::vector<Complex> mat31p = constructJumpOperator(N, num_states, 1, 1.0, pumping);
    JumpOperatorTester.testCyclicProperties(mat31p, matrix, num_states);

    std::vector<Complex> mat33z = constructJumpOperator(N, num_states, 3, 1.0, dephasing);
    JumpOperatorTester.testCyclicProperties(mat33z, matrix, num_states);

    
}
