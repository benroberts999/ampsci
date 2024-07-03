#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

void IsotopeShift(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
