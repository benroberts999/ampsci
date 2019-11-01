#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
class HartreeFock;
class UserInputBlock;
class DiracOperator;

namespace Module {

void matrixElements(const UserInputBlock &input, const Wavefunction &wf);

std::unique_ptr<DiracOperator> generateOperator(const std::string &operator_str,
                                                const UserInputBlock &input,
                                                const Wavefunction &wf);

} // namespace Module
