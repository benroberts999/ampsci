#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Example module, designed as a "template" to help you add a new module
void continuum(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
