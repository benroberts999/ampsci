#pragma once
#include "DiracOperator.hpp"
#include <string>
#include <vector>
class Wavefunction;
class HartreeFock;
class UserInputBlock;

namespace Module {

// void
std::vector<double> matrixElements(const UserInputBlock &input,
                                   const Wavefunction &wf);

ScalarOperator_old generateOperator(const std::string &operator_str,
                               const UserInputBlock &input,
                               const Wavefunction &wf);

} // namespace Module
