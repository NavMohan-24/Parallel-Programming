#include "ising.hpp"
#include <boost/dynamic_bitset.hpp>
#include <iomanip>

std::vector<Complex>  tfimHamiltonian(int N, double J, double h){

    // int num_states = (int) std::pow(2, N); - risky and can have floating imprecision

    int num_states = 1 << N;

    std::vector<Complex> hamiltonian(num_states*num_states);

    for (int k = 0 ; k < num_states; k++){
        for (int l = 0; l < num_states; l++){
            Complex ham(0.0,0.0);
            if (l == k){
                boost::dynamic_bitset<> b(N, k);
                for (int i = 0; i < N - 1; i++){
                    ham += Complex((b[i]==b[i+1]) ? -J:J,0);
                }
            }
            else{
                boost::dynamic_bitset<> bl(N, l); 
                boost::dynamic_bitset<> bk(N, k);
                boost::dynamic_bitset<> c = bl ^ bk;
                ham = Complex((c.count()==1) ? -h:0,0); 
            }
            hamiltonian[k*num_states+l] = ham;
        }
    }
    return hamiltonian;

}

#ifdef BUILD_MAIN
int main(){

    int N = 3;
    int num_states = 1 << N;
    double J = 1.0;
    double h = 1.0;

    std::vector<Complex> matrix = tfimHamiltonian(N, J, h);

    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << std::fixed << std::setprecision(2) 
          << matrix[i*num_states + j] << " "; 
        }
         std::cout << std::endl;
    }

}
#endif