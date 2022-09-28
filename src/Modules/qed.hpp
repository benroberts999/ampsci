#pragma once
#include <string>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
}
namespace DiracOperator {
class TensorOperator;
}

namespace Module {

//! Calculates QED corrections to energies and matrix elements
void QED(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
