#ifndef JUMP_OP
#define JUMP_OP

#include <vector>
#include <iostream>

enum class Scope {Local,Collective};
enum class DecayType {SigmaMinus, SigmaPlus, Dephasing};

//type alias complex data type
using Complex = std::complex<double>;
//using JumpOpsVec = std::vector<std::vector<Complex>>;

std::vector<Complex> constructJumpOperator(int N, int num_states,int k, DecayType type, Scope scope);

#endif