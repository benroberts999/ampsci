#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
// namespace HF {
// class HartreeFock;
// }
class UserInputBlock;
namespace DiracOperator {
class TensorOperator;
}

namespace Module {

void matrixElements(const UserInputBlock &input, const Wavefunction &wf);

std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const std::string &operator_str, const UserInputBlock &input,
                 const Wavefunction &wf);

} // namespace Module
