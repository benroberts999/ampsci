#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <string>
#include <vector>

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Just a test: example for playing with VQE
void VQE(const IO::InputBlock &input, const Wavefunction &wf);

//! Takes a subset of input basis according to subset_string.
//! Only states with energy>eFermi are included (to exclude core).
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      double eFermi);

} // namespace Module
