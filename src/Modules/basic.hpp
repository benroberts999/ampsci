#pragma once
#include <string>
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Module: writes orbitals to text file (gnuplot format)
void writeOrbitals(const IO::InputBlock &input, const Wavefunction &wf);

//! @brief A range of run-time tests (orthonorm, <H>, sum rules etc.) -- see
//! input options
//! @details Note: These are run-time tests to ensure input params were OK, not
//! unit tests! (see unitTests.cpp for unit tests)
void tests(const IO::InputBlock &input, const Wavefunction &wf);

//! Module to output continuum state wavefunctions to disk, and calculate matrix elements between them and bound states
void continuum(const IO::InputBlock &input, const Wavefunction &wf);

// Sub-functions for "tests"
namespace Tests {
void orthonormality(const Wavefunction &wf, const bool print_all = true);
void Hamiltonian(const Wavefunction &wf);
void r0pinf(const Wavefunction &wf);
void basis(const Wavefunction &wf);
} // namespace Tests

} // namespace Module