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

bool complex_equal(const Complex& x, const Complex& y, double eps = 1e-12) {
    return std::abs(x - y) < eps;
    };

bool vector_equal(const std::vector<Complex>& a,
                const std::vector<Complex>& b,
                double eps = 1e-12) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (!complex_equal(a[i], b[i], eps)) return false;
    }
    return true;
};
int main(){

    std::cout << "Running Unit tests..🔬" << "\n" << std::endl;
    
    int N = 1;
    int num_states = 1 << N;
    int k = 1;
    DecayType sigmaminus = DecayType::SigmaMinus;
    DecayType sigmaplus = DecayType::SigmaPlus;
    DecayType dephasing = DecayType::Dephasing;

    std::vector<Complex> matrix_1 = constructJumpOperator(N,num_states,k,sigmaminus);
    
    std::vector<Complex> matrix_2 = constructJumpOperator(N,num_states,k,sigmaplus);
    std::vector<Complex> matrix_3 = constructJumpOperator(N,num_states,k,dephasing);

    std::vector<Complex> SigmaMinus = {Complex(0,0), Complex(1,0), Complex(0,0), Complex(0,0)};
    std::vector<Complex> SigmaPlus = {Complex(0,0), Complex(0,0), Complex(1,0), Complex(0,0)};
    std::vector<Complex> SigmaZ = {Complex(1,0), Complex(0,0), Complex(0,0), Complex(-1,0)};

    if (vector_equal(matrix_1,SigmaMinus)){
        std::cout << "SigmaMinus decay type pass the unit test ✅" << std::endl;
    }
    else{
        std::cout << "SigmaMinus decay type fails the unit test ❌" << std::endl;
    }

    if (vector_equal(matrix_2,SigmaPlus)){
        std::cout << "SigmaPlus decay type pass the unit test ✅" << std::endl;
    }
    else{
        std::cout << "SigmaPlus decay type fails the unit test ❌" << std::endl;
    }

    if (vector_equal(matrix_3,SigmaZ)){
        std::cout << "SigmaZ decay type pass the unit test ✅" << std::endl;
    }
    else{
        std::cout << "SigmaZ decay type fails the unit test ❌" << std::endl;
    }

    

    






}
#endif

