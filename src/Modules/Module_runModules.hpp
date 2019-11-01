#pragma once
#include <string>
class Wavefunction;
class UserInputBlock;
class UserInput;

namespace Module {

void runModules(const UserInput &input, const Wavefunction &wf);

void runModule(const UserInputBlock &input, const Wavefunction &wf);

void Module_tests(const UserInputBlock &input, const Wavefunction &wf);
void Module_Tests_orthonormality(const Wavefunction &wf,
                                 const bool print_all = true);
void Module_Tests_Hamiltonian(const Wavefunction &wf);
void Module_WriteOrbitals(const UserInputBlock &input, const Wavefunction &wf);
void Module_BohrWeisskopf(const UserInputBlock &input, const Wavefunction &wf);
void Module_testPNC(const UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
