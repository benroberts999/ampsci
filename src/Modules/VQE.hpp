#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Just a test: example for playing with VQE
void VQE(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
