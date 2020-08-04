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

void testFeynman(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
