#pragma once
#include <string>
class Wavefunction;
class UserInput;
class HartreeFock;

namespace Module {

void runModule(const std::string &module, const UserInput &input,
               const Wavefunction &wf, const HartreeFock &hf);

void Module_tests(const UserInput &input, const Wavefunction &wf,
                  const HartreeFock &hf);
void Module_Tests_orthonormality(const Wavefunction &wf);
void Module_Tests_Hamiltonian(const Wavefunction &wf, const HartreeFock &hf);

void Module_WriteOrbitals(const UserInput &input, const Wavefunction &wf,
                          const HartreeFock &hf);

} // namespace Module
