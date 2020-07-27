#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
namespace IO {
class UserInputBlock;
}
namespace DiracOperator {
class TensorOperator;
}

namespace Module {

void HFAnomaly(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
