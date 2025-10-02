#include "hamiltonian.hpp"


std::vector<Complex>  constructHamiltonian(int N, int num_states, double J, double h){
    /**
    * @brief Construct the Hamiltonian for the TFIM with open boundary condition.
    *
    * @param N Number of spins.
    * @param num_states Hilbert space dimension (2^N).
    * @param J Interaction strength.
    * @param h Strength of transverse magnetic field.
    * @return Vectorized TFIM Hamiltonian.
    * 
    *
    */

    
    std::vector<Complex> hamiltonian(num_states*num_states);

    #pragma omp parallel for collapse(2) shared(hamiltonian)
    for (int k = 0 ; k < num_states; k++){
        for (int l = 0; l < num_states; l++){
            Complex ham(0.0,0.0);
            // bitset to represent basis state
            boost::dynamic_bitset<> bl(N, l); // bit operations without bitset (for eg; 1 << N) is not scalable

     
            if (l == k){
                // Generate diagonal elements
                for (int i = 0; i < N - 1; i++){
                    ham += Complex((bl[i]==bl[i+1]) ? -J:J,0);
                }
            }
            else{
                //Generate off-diagonal elements
                boost::dynamic_bitset<> bk(N, k);
                //bitwise XOR
                boost::dynamic_bitset<> c = bl ^ bk;
                //count number of non zero elements
                ham += Complex((c.count()==1) ? -h:0,0); 
            }
            hamiltonian[k*num_states+l] = ham;
        }
    }
    return hamiltonian;

}