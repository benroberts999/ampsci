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
namespace MBPT {
class FeynmanSigma;
}

namespace Module {

//! Calculates effective exchange screening factors (Correlation Potential)
void screeningFactors(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
