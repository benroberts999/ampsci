#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculates lifetimes (E1, E2, M1)
void lifetimes(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
