#pragma once
#include <vector>
// Forward declare classes:
class Wavefunction;
class DiracSpinor;
namespace IO {
class InputBlock;
}
namespace Coulomb {
class CoulombTable;
}
namespace Angular {
class SixJTable;
}

namespace Module {

//! Module for testing ladder diagram implementation
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

void check_L_symmetry(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &valence,
                      const Coulomb::CoulombTable &qk,
                      const Angular::SixJTable &sj,
                      const Coulomb::CoulombTable *const lk = nullptr);

} // namespace Module
