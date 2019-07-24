#pragma once
#include "DiracOperator.hpp"
#include <string>
class Wavefunction;
class HartreeFock;
class UserInput;

namespace Module {

void matrixElements(const std::string &module, const UserInput &input,
                    const Wavefunction &wf);

DiracOperator generateOperator(const std::string &operator_str,
                               const UserInput &input, const Wavefunction &wf);

} // namespace Module
