#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculating muonic wavefunctions, energies, matrix elements
void muon(const IO::InputBlock &input, const Wavefunction &wf);

//! Muonic parity violation tests
void muonPV(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
