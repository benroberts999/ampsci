#pragma once
#include "Coulomb/include.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace UnitTest::LadderTest {

//! Small set of hydrogen-like orbitals, with Coulomb integrals: a pure
//! angular/algebra test bed for the ladder iteration (no Hartree-Fock needed)
struct LadderTestData {
  //! core: n <= 2 (1s, 2s, 2p-, 2p+), so hole pairs have several J
  std::vector<DiracSpinor> core;
  //! excited: n = 3, 4 (up to f states)
  std::vector<DiracSpinor> excited;
  //! "valence": 3s, 3p+, 3d+ (members of excited, as in the ladder driver)
  std::vector<DiracSpinor> valence;
  //! core + valence
  std::vector<DiracSpinor> i_orbs;
  Coulomb::YkTable yk;
  Coulomb::QkTable qk;
};

//! Forms the test data; Qk is filled up to k_cut (<= 0: no cut-off). zeff:
//! nuclear charge of the orbitals; a larger zeff gives a ladder iteration that
//! converges (Coulomb integrals scale as zeff, energy differences as zeff^2)
inline LadderTestData ladder_test_data(int k_cut = -1, double zeff = 1.0) {
  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const auto orbitals = DiracSpinor::HlikeBasis(3, 3, radial_grid, zeff);
  LadderTestData data{{}, {}, {}, {}, Coulomb::YkTable(orbitals), {}};
  for (const auto &Fn : orbitals) {
    (Fn.n() <= 2 ? data.core : data.excited).push_back(Fn);
    const bool valence_kappa =
      Fn.kappa() == -1 || Fn.kappa() == -2 || Fn.kappa() == -3;
    if (Fn.n() == 3 && valence_kappa) {
      data.valence.push_back(Fn);
    }
  }
  data.i_orbs = qip::merge(data.core, data.valence);
  data.qk.fill(orbitals, data.yk, k_cut, false);
  return data;
}

//! Largest k of the stored integrals, as the ladder driver: the basis limit,
//! cut at the Coulomb table's maximum
inline int stored_max_k(const LadderTestData &data) {
  return std::min(DiracSpinor::max_tj(qip::merge(data.core, data.excited)),
                  data.qk.max_k());
}

} // namespace UnitTest::LadderTest
