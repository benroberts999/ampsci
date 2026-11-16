#pragma once
#include "Coulomb/include.hpp"
#include "MBPT/LadderIteration.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <memory>
#include <random>
#include <vector>

namespace UnitTest::LadderTest {

//! Small set of hydrogen-like orbitals, with Coulomb integrals: a pure
//! angular/algebra test bed for the ladder iteration (no Hartree-Fock needed)
struct LadderTestData {
  //! core: n <= 2 (1s, 2s, 2p-, 2p+), so hole pairs have several J
  std::vector<DiracSpinor> core;
  //! excited: n >= 3, from 3 states per kappa starting at n = l+1 (3s; 3p,
  //! 4p; 3d to 5d; 4f to 6f)
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

//! Ladder table with every L^k_mnib (m,n excited; i in i_orbs; b core; k <=
//! stored_max_k() [add_Lk_mnib_entries()]) a random number in [-0.05, 0.05]:
//! a dressed integral with both parities of k, unlike the Coulomb integrals
inline Coulomb::LkTable random_ladder_table(const LadderTestData &data) {
  Coulomb::LkTable lk;
  MBPT::add_Lk_mnib_entries(&lk, data.excited, data.core, data.i_orbs,
                            stored_max_k(data));
  std::mt19937 generator(1);
  std::uniform_real_distribution<double> random_value(-0.05, 0.05);
  for (std::size_t k = 0; k < lk->size(); ++k) {
    for (auto &entry : lk->at(k)) {
      entry.second = random_value(generator);
    }
  }
  return lk;
}

//! Closing table of x_F with every L'^k_vnab (v valence; n excited; a,b core;
//! k <= stored_max_k() [add_Lk_vnab_F_entries()]) a random number in
//! [-0.05, 0.05]: both parities of k, and no symmetry
inline Coulomb::NkTable random_closing_table(const LadderTestData &data) {
  Coulomb::NkTable L_vnab;
  MBPT::add_Lk_vnab_F_entries(&L_vnab, data.excited, data.core, data.valence,
                              stored_max_k(data));
  std::mt19937 generator(2);
  std::uniform_real_distribution<double> random_value(-0.05, 0.05);
  for (std::size_t k = 0; k < L_vnab->size(); ++k) {
    for (auto &entry : L_vnab->at(k)) {
      entry.second = random_value(generator);
    }
  }
  return L_vnab;
}

} // namespace UnitTest::LadderTest
