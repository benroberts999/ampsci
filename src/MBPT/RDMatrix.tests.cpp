#include "DiracODE/DiracODE.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <cassert>
#include <complex>
#include <numeric>

//==============================================================================
TEST_CASE("MBPT: RDMatrix", "[MBPT][RDMatrix][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: RDMatrix, [MBPT]\n";

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
    const bool include_g = false;
    const std::size_t i0 = 100;
    const std::size_t stride = 5;
    const std::size_t size = 150;
    // => max point is 850
    MBPT::RDMatrix<double> m(i0, stride, size, include_g, grid);

    // Construct Green function, G(e) = |a><a|/(e-e_a)
    double e = -0.1;
    for (const auto &Fs : orbitals) {
      m.add(Fs, Fs, 1.0 / (e - Fs.en()));
    }

    // Test "standard" operations: +, -, scalar multiplication
    auto m2 = 2.0 * m;
    const auto m3 = m + m;
    const auto m4 = 3.0 * m - m;
    assert(equal(m2, m3));
    assert(!equal(m2, 1.1 * m3));
    assert(equal(m2, m4));

    // Calculate <Fs|G|Fs> - should be equal to denominator for Fs part of M
    // <Fs|G|Fs> = 1.0 / (e - e_s)
    double worst = 0.0;
    for (const auto &Fs : orbitals) {
      const auto expected = 1.0 / (e - Fs.en());
      const auto res = Fs * (m.drj() * Fs);
      const auto eps = std::abs((res - expected) / expected);
      if (eps > worst)
        worst = eps;
    }
    // pass &= qip::check_value(&obuff, "<F|G|F> (small)", worst, 0.0, 1.0e-2);
    REQUIRE(std::abs(worst) < 1.0e-2);

    // Calculate <Fs|G*G|Fs> - should be equal to denominator^2 for Fs part of M
    // <Fs|G*S|Fs> = 1.0 / (e - e_s)^2
    // This tests matrix multiplication
    const auto mm = m.drj() * m;
    double worst_2 = 0.0;
    for (const auto &Fs : orbitals) {
      const auto expected = 1.0 / (e - Fs.en()) / (e - Fs.en());
      const auto res = Fs * (mm.drj() * Fs);
      const auto eps = std::abs((res - expected) / expected);
      if (eps > worst_2)
        worst_2 = eps;
    }
    // pass &= qip::check_value(&obuff, "<F|G*G|F> (small)", worst_2,
    // 0.0, 1.0e-1);
    REQUIRE(std::abs(worst_2) < 1.0e-1);

    // Test inverse (real, no g): two ways:
    {
      // Test inverting matrix, using Neumann expansion:
      // a = (1+lam*m)^-1
      // a =~ 1 + lam*m - lam^2*m*m + ...
      // Note: can't directly check a.inv()*a, since * involves integration!
      // (We check direct inversion below, using a "direct" matrix
      // multiplication)
      const double lambda = 1.0e-5;
      const auto a = (lambda * m + 1.0).inverse();
      // Since Neumann expansion is approximate, and only valid for small lambda
      // we check two properties
      // 1) Check (1+lam*m)^-1 is close to 1 + lam*m - lam^2*m*m
      // 2) Check that each term (up to lam^2) _improves_ neumann expansion

      // First three Neuman approximations
      const auto b0 = 0.0 * m + 1.0; // just identity
      const auto b1 = b0 + (-lambda) * m;
      const auto b2 = b1 + lambda * lambda * m * m;
      // Differences from expected: [nb: eps is bad, since values small]
      const auto del0 = max_delta(a, b0);
      const auto del1 = max_delta(a, b1);
      const auto del2 = max_delta(a, b2);

      // pass &= qip::check_value(&obuff, "G^-1 (small)", del2, 0.0, 1.0e-8);
      // pass &= qip::check(&obuff, "G^-1 Neumann converge",
      //                    (del2 < del1 && del1 < del0));
      REQUIRE(std::abs(del2) < 1.0e-8);
      REQUIRE((del2 < del1 && del1 < del0));

      // Now, use 'raw_mat_mul' (doesn't contain integration measure) to
      // directly test inverse. Note that we don't use 'raw_mat_mul' in any
      // physics..so only a test of method
      // +1.0 to ensure inverse exists..
      // const auto prod = raw_mat_mul((m + 1.0).inverse(), m + 1.0);
      const auto prod = (m + 1.0).inverse() * (m + 1.0);
      const auto ident = 0.0 * m + 1.0;
      const auto del_d = max_delta(prod, ident);
      // pass &=
      //     qip::check_value(&obuff, "G^-1 direct (small)", del_d,
      //     0.0, 1.0e-13);
      REQUIRE(std::abs(del_d) < 1.0e-13);
    }
  }

  //============================================================================
  // Test basics using "small" matrix (but, include G)
  {
    const bool include_g = true;
    const std::size_t i0 = 100;
    const std::size_t stride = 5;
    const std::size_t size = 150;
    MBPT::RDMatrix<double> m(i0, stride, size, include_g, grid);

    // Construct Green function, G(e) = |a><a|/(e-e_a)
    const double e = -0.1;
    for (const auto &Fs : orbitals) {
      m.add(Fs, Fs, 1.0 / (e - Fs.en()));
    }

    // Calculate <Fs|G*G|Fs> - should be equal to denominator^2 for Fs part of M
    // This tests matrix multiplication, including g
    double worst = 0.0;
    for (const auto &Fs : orbitals) {
      const auto expected = 1.0 / (e - Fs.en());
      const auto res = Fs * (m.drj() * Fs);
      const auto eps = std::abs((res - expected) / expected);
      if (eps > worst)
        worst = eps;
    }
    // pass &= qip::check_value(&obuff, "<F|G|F> (small+g)", worst,
    // 0.0, 1.0e-4);
    REQUIRE(std::abs(worst) < 1.0e-4);

    // Matrix inversion (real, including g):
    {
      // Test inverting matrix, using Neumann expansion:
      const double lambda = 1.0e-5;
      const auto a = (lambda * m + 1.0).inverse();
      // First three Neuman approximations
      const auto b0 = 0.0 * m + 1.0; // just identity
      const auto b1 = b0 + (-lambda) * m;
      const auto b2 = b1 + lambda * lambda * m * m;
      const auto del0 = max_delta(a, b0);
      const auto del1 = max_delta(a, b1);
      const auto del2 = max_delta(a, b2);

      // pass &= qip::check_value(&obuff, "G^-1 (small+g)", del2, 0.0, 1.0e-8);
      // pass &= qip::check(&obuff, "G^-1 Neumann converge",
      //                    (del2 < del1 && del1 < del0));
      REQUIRE(std::abs(del2) < 1.0e-8);
      REQUIRE((del2 < del1 && del1 < del0));

      // Now, use 'raw_mat_mul' (doesn't contain integration measure) to
      // directly test inverse.
      const auto prod = ((m + 1.0).inverse()) * (m + 1.0);
      const auto ident = 0.0 * m + 1.0;
      const auto del_d = max_delta(prod, ident);
      // pass &= qip::check_value(&obuff, "G^-1 direct (small+g)", del_d, 0.0,
      //  1.0e-13);
      REQUIRE(std::abs(del_d) < 1.0e-13);
    }
  }

  //============================================================================
  // Test basics using "small" matrix (but, include G), and complex!
  {
    using namespace std::complex_literals;
    const bool include_g = true;
    const std::size_t i0 = 100;
    const std::size_t stride = 5;
    const std::size_t size = 150;
    MBPT::RDMatrix<std::complex<double>> m(i0, stride, size, include_g, grid);

    // Construct Green function, G(e) = |a><a|/(e-e_a)
    const double e = -0.1;
    for (const auto &Fs : orbitals) {
      m.add(Fs, Fs, (2.0 + 3.0i) / (e - Fs.en()));
    }

    const auto mr = m.real();
    const auto mi = m.imag();

    // Calculate <Fs|G|Fs> - should be equal to denominator for Fs part of M
    double worst = 0.0;
    for (const auto &Fs : orbitals) {
      const auto expected_re = 2.0 / (e - Fs.en());
      const auto expected_im = 3.0 / (e - Fs.en());
      const auto res_re = Fs * (mr.drj() * Fs);
      const auto res_im = Fs * (mi.drj() * Fs);
      const auto eps_re = std::abs((res_re - expected_re) / expected_re);
      const auto eps_im = std::abs((res_im - expected_im) / expected_im);
      const auto eps = std::max(eps_re, eps_im);
      if (eps > worst)
        worst = eps;
    }
    // pass &= qip::check_value(&obuff, "<F|iG|F> (small+g)", worst,
    // 0.0, 1.0e-4);
    REQUIRE(std::abs(worst) < 1.0e-4);

    // Test multiplication with imag. matrix:
    // const auto mm = m.drj() * m;
    const auto mm = m * m.dri();
    const auto mm_re = mm.real();
    double worst2 = 0.0;
    for (const auto &Fs : orbitals) {
      const auto v = (2.0 + 3.0i) / (e - Fs.en());
      const auto expected = std::real(v * v);
      const auto res = Fs * (mm_re.drj() * Fs);
      const auto eps = std::abs((res - expected) / expected);
      if (eps > worst2)
        worst2 = eps;
    }
    // pass &=
    //     qip::check_value(&obuff, "<F|iG*iG|F> (small+g)", worst2,
    //     0.0, 1.0e-4);
    REQUIRE(std::abs(worst2) < 1.0e-4);

    // Matrix inversion (complex, including g):
    {
      // Test inverting matrix, using Neumann expansion:
      const double lambda = 1.0e-5;
      const auto a = (lambda * m + 1.0).inverse();
      // First three Neuman approximations
      const auto b0 = 0.0 * m + 1.0; // just identity
      const auto b1 = b0 + (-lambda) * m;
      const auto b2 = b1 + lambda * lambda * m * m;
      const auto del0 = max_delta(a, b0);
      const auto del1 = max_delta(a, b1);
      const auto del2 = max_delta(a, b2);

      // pass &= qip::check_value(&obuff, "iG^-1 (small+g)", del2, 0.0, 1.0e-6);
      // pass &= qip::check(&obuff, "iG^-1 Neumann converge",
      //                    (del2 < del1 && del1 < del0));
      REQUIRE(std::abs(del2) < 1.0e-6);
      REQUIRE((del2 < del1 && del1 < del0));

      // Now, use 'raw_mat_mul' (doesn't contain integration measure) to
      // directly test inverse.
      const auto prod = ((m + 1.0).inverse()) * (m + 1.0);
      const auto ident = 0.0 * m + 1.0;
      const auto del_d = max_delta(prod, ident);
      // pass &= qip::check_value(&obuff, "iG^-1 direct (small+g)", del_d, 0.0,
      //                          1.0e-13);
      REQUIRE(std::abs(del_d) < 1.0e-13);
    }
  }

  //============================================================================
  // Test basics using "full" matrix (no stride, include G)
  {
    const bool include_g = true;
    const std::size_t i0 = 0;
    const std::size_t stride = 1;
    const std::size_t size = grid->num_points();
    MBPT::RDMatrix<double> m(i0, stride, size, include_g, grid);

    // Construct Green function, G(e) = |a><a|/(e-e_a)
    const double e = -0.1;
    for (const auto &Fs : orbitals) {
      m.add(Fs, Fs, 1.0 / (e - Fs.en()));
    }

    // Test "standard" operations: +, -, scalar multiplication, for matrix w/ G
    auto m2 = 2.0 * m;
    const auto m3 = m + m;
    const auto m4 = 3.0 * m - m;
    assert(equal(m2, m3));
    assert(!equal(m2, 1.1 * m3));
    assert(equal(m2, m4));

    // Calculate <Fs|G|Fs> - should be equal to denominator for Fs part of M
    double worst = 0.0;
    for (const auto &Fs : orbitals) {
      const auto expected = 1.0 / (e - Fs.en());
      const auto res = Fs * (m.drj() * Fs);
      const auto eps = std::abs((res - expected) / expected);
      if (eps > worst)
        worst = eps;
    }
    // pass &= qip::check_value(&obuff, "<F|G|F> (full)", worst, 0.0, 1.0e-14);
    REQUIRE(std::abs(worst) < 1.0e-14);

    // Test multiplication, including G (full matrix)
    const auto mm = m.drj() * m;
    double worst_2 = 0.0;
    for (const auto &Fs : orbitals) {
      const auto expected = 1.0 / (e - Fs.en()) / (e - Fs.en());
      const auto res = Fs * (mm.drj() * Fs);
      const auto eps = std::abs((res - expected) / expected);
      if (eps > worst_2)
        worst_2 = eps;
    }
    // pass &= qip::check_value(&obuff, "<F|G*G|F> (full)", worst_2,
    // 0.0, 1.0e-12);
    REQUIRE(std::abs(worst_2) < 1.0e-12);
  }
}
