#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Module to output continuum state wavefunctions to disk
void continuum(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
