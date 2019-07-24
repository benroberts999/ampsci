#pragma once
#include "DiracOperator.hpp"
#include <string>
class Wavefunction;
class HartreeFock;
class UserInputBlock;

namespace Module {

void matrixElements(const UserInputBlock &input, const Wavefunction &wf);

DiracOperator generateOperator(const std::string &operator_str,
                               const UserInputBlock &input,
                               const Wavefunction &wf);

} // namespace Module
