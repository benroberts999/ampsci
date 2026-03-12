#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

void GreenQED(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
