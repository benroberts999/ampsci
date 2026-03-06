#pragma once
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Example module, designed as a "template" to help you add a new module. Note:
//! if you add a new Module, you must also update module_lists.hpp, so it will
//! be compiled into the rest of the code.

// struct SigmaData {
//   int kappa;
//   double en;
//   SpinorMatrix<double> Sigma;
//   int n{0};
//   double lambda{1.0};
// };

// struct rgrid_params {
//   double r0{1.0e-4};
//   double rmax{30.0};
//   std::size_t stride{4};
// };

void StructRad(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
