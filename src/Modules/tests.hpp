#pragma once
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! @brief A range of run-time tests (orthonorm, <H>, sum rules etc.) -- see
//! input options
//! @details Note: These are run-time tests to ensure input params were OK, not
//! unit tests! (see unitTests.cpp for unit tests)
void Module_tests(const IO::InputBlock &input, const Wavefunction &wf);

namespace Tests {
void Module_Tests_orthonormality(const Wavefunction &wf,
                                 const bool print_all = true);
void Module_Tests_Hamiltonian(const Wavefunction &wf);
void Module_test_r0pinf(const Wavefunction &wf);
void basisTests(const Wavefunction &wf);
} // namespace Tests

} // namespace Module
