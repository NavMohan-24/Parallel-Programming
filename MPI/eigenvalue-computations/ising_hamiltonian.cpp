#include <iostream>
#include <vector>
#include <boost/dynamic_bitset.hpp>

// using namespace std -- can be thought analogus to "from __ import*" in python.


std::vector<std::vector<double>>  isingHamiltonian(int N, double J, double h){

    // int num_states = (int) std::pow(2, N); - risky and can have floating imprecision

    int num_states = 1 << N;

    std::vector<std::vector<double>> hamiltonian(num_states, std::vector<double>(num_states,0));

    for (int k = 0 ; k < num_states; k++){
        for (int l = 0; l < num_states; l++){
            double ham = 0.0;
            if (l == k){
                boost::dynamic_bitset<> b(N, k);
                for (int i = 0; i < N - 1; i++){
                    ham += (b[i]==b[i+1]) ? -J:J;
                }
            }
            else{
                boost::dynamic_bitset<> bl(N, l);
                boost::dynamic_bitset<> bk(N, k);
                boost::dynamic_bitset<> c = bl ^ bk;
                ham = (c.count()==1) ? -h:0; 
            }
            hamiltonian[l][k] = ham;
        }
    }
    return hamiltonian;

}

int main(){

    int N = 3;
    int num_states = 1 << N;
    double J = 1.0;
    double h = 1.0;

    std::vector<std::vector<double>> matrix = isingHamiltonian(N, J, h);

    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << matrix[i][j] << " "; 
        }
         std::cout << std::endl;
    }

}