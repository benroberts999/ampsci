#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class UserInputBlock;
}

namespace Module {

//! Example module, designed as a "template" to help you add a new module
void exampleModule(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
