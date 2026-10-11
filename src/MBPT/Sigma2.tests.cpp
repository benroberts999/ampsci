#include "Sigma2.hpp"
#include "CI/CI_Integrals.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include <array>
#include <memory>
#include <utility>
#include <vector>

//==============================================================================
TEST_CASE(
  "MBPT: Sigma2 agrees with finite magnetic-substate perturbation theory",
  "[Sigma2][Sigma2Phase][MBPT][CI][unit]") {
  const auto grid = std::make_shared<const Grid>(
    GridParameters{6000, 1.0e-8, 80.0, 1.0, GridType::loglinear});
  std::vector<DiracSpinor> basis;
  for (const auto [n, kappa] : std::array<std::pair<int, int>, 4>{
         {{1, -1}, {2, -1}, {2, 1}, {2, -2}}}) {
    auto f = DiracSpinor::exactHlike(n, kappa, grid, 3.0);
    // Equal valence energies remove denominator-convention ambiguities.
    f.en() = n == 1 ? -2.0 : -0.4;
    basis.push_back(f);
  }
  const std::vector<DiracSpinor> core{basis.front()};
  const std::vector<DiracSpinor> excited(basis.begin() + 1, basis.end());
  const Coulomb::YkTable yk(basis);
  Coulomb::QkTable qk;
  qk.fill(basis, yk);
  const Angular::SixJTable sixj(5);
  Coulomb::LkTable sigma;
  for (const auto &v : excited) {
    for (const auto &w : excited) {
      for (const auto &x : excited) {
        for (const auto &y : excited) {
          const auto [kmin, kmax] = MBPT::k_minmax_S(v, w, x, y);
          for (int k = kmin; k <= kmax; ++k) {
            sigma.add(k, v, w, x, y,
                      MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited, sixj,
                                    MBPT::Denominators::Fermi0));
          }
        }
      }
    }
  }

  // References from direct second-order perturbation theory in the 10
  // magnetic orbitals of this four-spinor basis (no reduced Sigma2 formulas):
  //   Heff^(2) = P V Q (E0 - Q H0 Q)^(-1) Q V P.
  // V is Coulomb normal ordered about 1s^2. Subtract the closed-core term
  // and the lifted one-valence Sigma1 before projecting onto two-electron
  // CSFs. Indices below refer to excited = {2s1/2, 2p1/2, 2p3/2}.
  // Both parities and J=0,1,2 are included, with mixed-j and unchanged controls.
  struct Reference {
    int twoJ;
    std::size_t v, w, x, y;
    double value;
  };
  const std::array<Reference, 20> reference{
    {{0, 0, 1, 0, 1, 0.0036264422250935329},
     {0, 0, 0, 0, 0, -0.00052655138788083432},
     {0, 0, 0, 1, 1, 0.0054072802580107937},
     {0, 0, 0, 2, 2, 0.0076433585756673431},
     {0, 1, 1, 0, 0, 0.0054072802580107940},
     {0, 1, 1, 1, 1, -0.0055889833124979820},
     {0, 1, 1, 2, 2, 0.00059899231113111980},
     {0, 2, 2, 0, 0, 0.0076433585756673413},
     {0, 2, 2, 1, 1, 0.00059899231113111969},
     {0, 2, 2, 2, 2, -0.0051605006333440059},
     {2, 0, 1, 0, 1, 0.00072470246367442831},
     {2, 0, 1, 0, 2, -0.0041030630325135254},
     {2, 0, 2, 0, 1, -0.0041030630325135245},
     {2, 0, 2, 0, 2, -0.0021766456723840608},
     {2, 1, 2, 1, 2, -0.0060100909018721109},
     {4, 0, 2, 0, 2, 0.0036250896827350531},
     {4, 1, 2, 1, 2, -0.0035237884604664305},
     {4, 1, 2, 2, 2, -0.0017580182415392327},
     {4, 2, 2, 1, 2, -0.0017580182415392331},
     {4, 2, 2, 2, 2, -0.0047646303898725133}}};
  for (const auto &r : reference) {
    CAPTURE(r.twoJ, r.v, r.w, r.x, r.y);
    const auto actual =
      CI::CSF2_Sigma2(sigma, excited[r.v].nk_index(), excited[r.w].nk_index(),
                      excited[r.x].nk_index(), excited[r.y].nk_index(), r.twoJ);
    CHECK(actual == Approx(r.value).epsilon(0.0).margin(1.0e-12));
  }
}
