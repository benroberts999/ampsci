#pragma once
class Wavefunction;
namespace IO {
class UserInputBlock;
}

namespace Module {

//! Runs a range of tests (orthonorm, <H>, sum rules etc.) -- see input options
void Module_tests(const IO::UserInputBlock &input, const Wavefunction &wf);

namespace Tests {
void Module_Tests_orthonormality(const Wavefunction &wf,
                                 const bool print_all = true);
void Module_Tests_Hamiltonian(const Wavefunction &wf);
void Module_test_r0pinf(const Wavefunction &wf);
void basisTests(const Wavefunction &wf);
} // namespace Tests

} // namespace Module
