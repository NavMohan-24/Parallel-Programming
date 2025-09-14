#include "jump_op.hpp"
#include <boost/dynamic_bitset.hpp>

std::vector<Complex> constructJumpOperator(int N, int num_states, int k, DecayType type = DecayType::SigmaMinus, Scope scope = Scope::Local){

    std::vector<Complex> jump_k(num_states*num_states, Complex(0,0)); 

    switch(scope)
    {
    case Scope::Local:
        if (type == DecayType::SigmaMinus){
            for (int m = 0; m < num_states; m++){
                int n = m + (1 << (k-1));
                if (n < num_states){
                    boost::dynamic_bitset<> bn(N, n);
                    jump_k[m*num_states+n] = Complex((bn.test(k-1)==true) ? 1:0, 0);
                };
            };
            return jump_k;
        }
        if (type == DecayType::SigmaPlus){
            for (int m = 0; m < num_states; m++){
                int n = m - (1 << (k-1));
                if (n >= 0){
                    boost::dynamic_bitset<> bn(N, n);
                    jump_k[m*num_states+n] = Complex((bn.test(k-1)==false) ? 1:0, 0);
                };
            };
            return jump_k;
        }
        if (type == DecayType::Dephasing){
            for (int m =0; m < num_states; m++){
                boost::dynamic_bitset<> bn(N, m);
                if ((k-1) < N){
                    jump_k[m*num_states+m] = Complex((bn.test(k-1)==true) ? -1:1, 0);
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
int main(){

    int N = 3;
    int num_states = 1 << N;
    int k = 2;
    DecayType dephasing = DecayType::Dephasing;

    std::vector<Complex> matrix = constructJumpOperator(N,num_states,k,dephasing);

    std::cout << "Printing Jump matrix.." << "\n" << std::endl;
    for (int i = 0; i < num_states; i++){
        for (int j = 0; j < num_states; j++){
            std::cout << std::fixed << std::setprecision(3) << matrix[i*num_states+j] << "\t"; 
        }
        std::cout << std::endl;
    };

}
#endif

