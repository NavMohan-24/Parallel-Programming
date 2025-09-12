#ifndef JUMP_OP
#define JUMP_OP

#include <vector>
#include <iostream>

//type alias complex data type
using Complex = std::complex<double>;

std::vector<Complex> constructJumpOperator(int num_states, std::string type, double r);

#endif