#ifndef JUMP_OP
#define JUMP_OP

#include <vector>
#include <iostream>
#include <complex>
#include <boost/dynamic_bitset.hpp>

enum class Scope {Local,Collective};
enum class DecayType {Damping, Pumping, Dephasing}; // Sigma Minus, Sigma Plus, Sigma Z

//type alias complex data type
using Complex = std::complex<double>;
//using JumpOpsVec = std::vector<std::vector<Complex>>;

std::vector<Complex> constructJumpOperator(int N, int num_states, int k, double rate = 1.0, DecayType type = DecayType::Damping , Scope scope = Scope::Local);

#endif
