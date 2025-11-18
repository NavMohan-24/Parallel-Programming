//=======================================================================================
// TFIM  Linbladian Spectrum
//=======================================================================================

#include "density_matrix/density_matrix.hpp"
#include "hamiltonian/hamiltonian.hpp"
#include "jump_operator/jump_operator.hpp"
#include "linbladian_diagonalizer.hpp"
#include <fstream>
#include <format>

int main(){

    // System parameters
    int N = 3;                    // Number of qubits
    int num_states = 1 << N;      // Hilbert space dimension (2^N)
    int K = 40;                   // Krylov subspace dimension
    double J = -1.0;              // Coupling strength
    double h = 0.1;               // Magnetic field strength

    // Initialize with random density matrix
    std::vector<Complex> rho = constructRandomRho(num_states);


    // Construct Lindbladian and diagonalize
    std::vector<Complex> hamiltonian = constructHamiltonian(N, num_states, J, h);
    LinbladianConstructor linb(hamiltonian, N);
    ArnoldiLinbladianDiagonalizer LSolver(linb, N, 2);
    EigenResult result = LSolver.diagonalize(rho, K);

    // Display results
    std::cout << "\n=== Ritz Vectors (Eigenvectors) ===\n\n";
    int size = static_cast<int>(result.eigenvalues.size());
    for (int i = 0; i < size; i++){
        std::cout << "Eigenvector " << i+1 << ":\n";
        for (int j = 0; j < num_states*num_states; j++){
            std::cout << result.eigenvectors[j*size + i] << "\n";
        }
        std::cout << "\n";
    }
    std::string filename = std::format("eigVals{}.csv", K);
    std::ofstream file(filename);
    std::cout << "=== Ritz Values (Eigenvalues) ===\n\n";
    for(int k = 0; k < size; k++){
        std::cout << "Eigenvalue " << k+1 << ": " << result.eigenvalues[k] << "\n";
        file << result.eigenvalues[k].real() << "," << result.eigenvalues[k].imag() << "\n";

    }
    file.close();


}