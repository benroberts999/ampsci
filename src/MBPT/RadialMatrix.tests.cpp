#include "RadialMatrix.hpp"
#include "DiracODE/include.hpp"
#include "Maths/Grid.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <cassert>
#include <complex>
#include <numeric>

//==============================================================================
TEST_CASE("MBPT: RadialMatrix", "[MBPT][RadialMatrix][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: RadialMatrix, [MBPT]\n";

  // Set up radial grid:
  const auto r0{1.0e-7};
  const auto rmax{100.0}; // NB: rmax depends on Zeff
  const auto num_grid_points{1000ul};
  const auto b{10.0};
  const auto grid = std::make_shared<const Grid>(r0, rmax, num_grid_points,
                                                 GridType::loglinear, b);

  // First, just create a list of Hydgrogenlike (s) orbitals:
  // States to solve for:
  const int n_max = 10;
  const double Zeff = 5.0;
  // Sperical potential w/ R_nuc = 0.0 is a pointlike potential
  const auto v_nuc = Nuclear::sphericalNuclearPotential(Zeff, 0.0, grid->r());
  // Solve Dirac ODE for each state, store in 'orbitals' vector:
  std::vector<DiracSpinor> orbitals;
  for (int n = 1; n <= n_max; ++n) {
    auto &Fnk = orbitals.emplace_back(n, -1, grid);
    const auto en_guess = -(Zeff * Zeff) / (2.0 * n * n);
    DiracODE::boundState(Fnk, en_guess, v_nuc);
  }
  // What the orbitals are is not important; but we will use the fact that they
  // are an orthonormal set (tested in DiracODE)

  //============================================================================
  // Test basics using "small" matrix (stride=5, no G)
  {
    const std::size_t i0 = 100;
    const std::size_t stride = 5;
    const std::size_t size = 150;
    // => max point is 850
    MBPT::RadialMatrix<double> m(i0, stride, size, grid);

    // Construct Green function, G(e) = |a><a|/(e-e_a)
    double e = -0.1;
    for (const auto &Fs : orbitals) {
      for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
          const auto i_f = m.index_to_fullgrid(i);
          const auto j_f = m.index_to_fullgrid(j);
          m(i, j) += Fs.f(i_f) * Fs.f(j_f) / (e - Fs.en());
        }
      }
    }

    // Test "standard" operations: +, -, scalar multiplication
    auto m2 = 2.0 * m;
    const auto m3 = m + m;
    const auto m4 = 3.0 * m - m;
    REQUIRE(equal(m2, m3));
    REQUIRE(!equal(m2, 1.1 * m3));
    REQUIRE(equal(m2, m4));

    // Test inversion:
    // m not invertable, mM is:
    const auto mM = m + 100.0;
    auto m5 = mM;
    m5.invert_in_place();
    const auto m6 = mM.inverse();
    REQUIRE(equal(m5, m6));
    REQUIRE(!equal(m, m5));

    auto m7 = m6 * mM;
    auto m8 = mM * m6;
    const auto det8 = m8.Rmatrix().determinant();
    const auto det7 = m7.Rmatrix().determinant();
    REQUIRE(det8 == Approx(1.0));
    REQUIRE(det8 == Approx(det7));

    m7 -= 1.0;
    m8 -= 1.0;
    const auto max_element2 = MBPT::max_element(m7);
    const auto max_element3 = MBPT::max_element(m8);
    REQUIRE(max_element2 < 1.0e-12);
    REQUIRE(max_element3 < 1.0e-12);
  }
}