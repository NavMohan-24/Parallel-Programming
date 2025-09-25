#include "jump_operator.hpp"


std::vector<Complex> constructJumpOperator(int N, int num_states, int k, double rate , DecayType type , Scope scope){

    //initialise jump operator with zero
    std::vector<Complex> jump_k(num_states*num_states, Complex(0,0)); 

    //square root of rate of decay
    double r = std::pow(rate,0.5);

    switch(scope)
    {
    case Scope::Local:
       // local damping via sigmaminus
        if (type == DecayType::Damping){

            #pragma omp parallel for shared(jump_k)
            for (int m = 0; m < num_states; m++){
                // finding states for which <m|L_k|n> is not zero.
                int n = m + (1 << (k-1)); // for instance let k=1, then n = m +1, such that L_1|n> = |m>
                if (n < num_states){
                    boost::dynamic_bitset<> bn(N, n);
                    jump_k[m*num_states+n] = Complex((bn.test(k-1)==true) ? r*1:0, 0); //the bit at k-1 th position should be one to apply jump operator.
                };
            };
            return jump_k;
        }
        if (type == DecayType::Pumping){
            #pragma omp parallel for shared(jump_k)
            for (int m = 0; m < num_states; m++){
                // finding states for which <m|L_k|n> is not zero.
                int n = m - (1 << (k-1)); // for instance let k=1, then n = m - 1, such that L_1|n> = |m>
                if (n >= 0){
                    boost::dynamic_bitset<> bn(N, n);
                    jump_k[m*num_states+n] = Complex((bn.test(k-1)==false) ? r*1:0, 0); //the bit at k-1 th position should be one to apply jump operator.
                };
            };
            return jump_k;
        }
        if (type == DecayType::Dephasing){
            #pragma omp parallel for shared(jump_k)
            for (int m =0; m < num_states; m++){
                // finding states for which <m|L_k|n> is not zero.
                boost::dynamic_bitset<> bn(N, m);
                if ((k-1) < N){
                    jump_k[m*num_states+m] = Complex((bn.test(k-1)==true) ? -r*1:r*1, 0);
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