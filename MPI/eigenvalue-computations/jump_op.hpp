#ifndef JUMP_OP
#define JUMP_OP

#include <vector>
#include <iostream>

enum class Scope {Local,Collective};
enum class DecayType {Damping, Pumping, Dephasing}; // Sigma Minus, Sigma Plus, Sigma Z

//type alias complex data type
using Complex = std::complex<double>;
//using JumpOpsVec = std::vector<std::vector<Complex>>;

std::vector<Complex> constructJumpOperator(int N, int num_states, int k, double rate, DecayType type, Scope scope);

#endif