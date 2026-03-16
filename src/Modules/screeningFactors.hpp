#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculates effective exchange screening factors (Correlation Potential)
void screeningFactors(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
