#pragma once
#include <string>
class Wavefunction;
namespace IO {
class InputBlock;
} // namespace IO

//! @brief Modules are run using calculated atomic wavefunctions
/*! @details
After the pogram has generated wavefunctions (Hartree-Fock etc), any number of
modules can then be run (see input documentation for how to run them).
Examples include: calculating matrix elements, lifetimes, PNC etc.

To add a new module:
  - a) Write the function (in the Module namespace). You can either add this to
one of the existing module cpp/hpp files, or make a new file.
  - b) Add your module function name to the module_list vector [in
Modules/modules_list.hpp]
  - c) If you added a new file, add the new .hpp file to the #include's list in
Modules/modules_list.hpp
  - ** See "Modules/exampleModule.{c/h}pp" for an example: you may start from
there using that one as a template
*/
namespace Module {

//! Loops through all given modules, runs them one at a time
void runModules(const IO::InputBlock &input, const Wavefunction &wf);

//! Figures out which module to run (Must be updated for each new module!)
void runModule(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
