#include <iostream>
#include <vector>
#include <boost/dynamic_bitset.hpp>

// using namespace std -- can be thought analogus to "from __ import*" in python.


std::vector<std::vector<double>>  isingHamiltonian(int N, double J){

    // int num_states = (int) std::pow(2, N); - risky and can have floating imprecision

    int num_states = 1 << N;

    std::vector<std::vector<double>> hamiltonian(num_states, std::vector<double>(num_states,0));

    for (int k = 0 ; k < num_states; k++){
        for (int l = 0; l < num_states; l++){
            if (l == k){
                double hzz = 0.0;
                boost::dynamic_bitset<> b(N, k);
                // for (int i = 0; i <= b.size()-1; i++){
                //     hzz -= J*std::pow(-1,b[i]+b[i+1]);
                // }
                for (int i = 0; i < N - 1; i++){
                    if (b[i] == b[i+1]) {
                        hzz -= J;   // same spin → σᶻσᶻ = +1
                    } 
                    else {
                    hzz += J;   // opposite spins → σᶻσᶻ = -1
                    }
                }
                hamiltonian[l][k] = hzz;
            }

        }
    }
    return hamiltonian;

}

int main(){

    int N = 4;
    int num_states = 1 << N;
    double J = 1.0;

    std::vector<std::vector<double>> matrix = isingHamiltonian(N, J);

    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << matrix[i][j] << " "; 
        }
         std::cout << std::endl;
    }

}