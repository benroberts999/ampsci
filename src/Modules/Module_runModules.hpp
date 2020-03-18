#pragma once
#include <string>
class Wavefunction;
class UserInputBlock;
class UserInput;

//! @brief Modules are run using calculated atomic wavefunctions
/*! @details
After the pogram has generated wavefunctions (Hartree-Fock etc), any number of
modules can then be run (see input documentation for how to run them).
Examples include: calculating matrix elements, lifetimes, PNC etc.

To add a new module:
  - a) Write the function (in the Module namespace). You can either add this to
one of the existing module cpp/hpp files, or make a new file.
  - b) update the runModule() function [in Modules/Module_runModules.cpp] to
call your module
  - c) If you have written youyr module in a new file, you'll also need to
update Module.mk file. Add both the compile/Dependencies option, and add the
'.o' name to the MODULELIST (should be fairly obvious by copying existing
examples). If you added your function to an existing file, this isn't needed.
*/
namespace Module {

//! Loops through all given modules, runs them one at a time
void runModules(const UserInput &input, const Wavefunction &wf);

//! Figures out which module to run (Must be updated for each new module!)
void runModule(const UserInputBlock &input, const Wavefunction &wf);

//! Module: writes orbitals to text file (gnuplot format)
void writeOrbitals(const UserInputBlock &input, const Wavefunction &wf);

//! Module: Calculates 2nd order energy shift for valence states
void SecondOrder(const UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
