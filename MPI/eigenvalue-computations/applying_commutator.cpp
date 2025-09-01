#include "ising.hpp"
#include <vector>
#include <cblas.h>


std::vector<std::vector<double>> applyCommutator(std::vector<std::vector<double>> H, std::vector<std::vector<double>> rho){

    

    



}

#ifdef BUILD_MAIN
int main(){

    int N = 3;
    double J = 1.0; 
    double h = 1.0;

    std::vector<std::vector<double>> H = tfimHamiltonian(N,J,h);
    std::vector<std::vector<double>> rho = constructRandomRho(N);
    std::vector<std::vector<double>> matrix = applyCommutator(H,rho);

    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << matrix[i][j] << " "; 
        }
         std::cout << std::endl;
    }

}