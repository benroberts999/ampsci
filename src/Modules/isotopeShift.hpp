#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculates field shift: F = d(E)/d(<r^2>)
void fieldShift(const IO::InputBlock &input, const Wavefunction &wf);

void fieldShift2(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
