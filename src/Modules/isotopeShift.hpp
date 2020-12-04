#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class UserInputBlock;
}

namespace Module {

//! Calculates field shift: F = d(E)/d(<r^2>)
void fieldShift(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
