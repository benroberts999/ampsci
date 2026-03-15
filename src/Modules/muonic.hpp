#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Muonic parity violation tests
void muonPV(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
