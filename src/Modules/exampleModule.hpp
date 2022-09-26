#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Example module, designed as a "template" to help you add a new module. Note:
//! if you add a new Module, you must also update module_lists.hpp, so it will
//! be compiled into the rest of the code.
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
