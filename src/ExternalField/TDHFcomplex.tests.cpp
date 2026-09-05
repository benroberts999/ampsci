#include "TDHFcomplex.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <cmath>
#include <complex>
#include <string>
#include <vector>

//==============================================================================
//! Interface checks: channel bookkeeping, amplitudes, warm restart, clear.
//! Ne: the 2p threshold is 0.85 au and 2s is 1.93 au, so omega = 0.5 has no
//! open channel and omega = 1.5 ionises 2p only.
TEST_CASE("TDHFcntm: basic unit tests",
          "[ExternalField][TDHF][TDHFcntm][unit]") {

  // The grid must resolve the continuum at en_+ = en_2p + omega ~ 0.65 au
  // out to rmax (wavelength ~5.5 au)
  Wavefunction wf({2000, 1.0e-6, 30.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-10, false);

  const auto E1 = DiracOperator::E1(wf.grid());
  auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());

  const auto *F1s = wf.getState("1s");
  const auto *F2p = wf.getState("2p+");
  REQUIRE(F1s != nullptr);
  REQUIRE(F2p != nullptr);

  // Before any solve: no open channels, no amplitudes, no correction
  REQUIRE(rpa.channel_list().empty());
  REQUIRE(rpa.A_phys().empty());
  REQUIRE(rpa.dV_complex(*F1s, *F2p) == std::complex<double>{0.0, 0.0});

  // Below every threshold: every channel is closed, and the imaginary parts
  // are exactly zero (nothing ever feeds them)
  rpa.solve_core(0.5, 100, false);
  REQUIRE(rpa.channel_list().empty());
  REQUIRE(rpa.A_phys().empty());
  REQUIRE(rpa.A_phys(*F2p, -1) == std::complex<double>{0.0, 0.0});
  for (const auto &Fa : wf.core()) {
    for (const auto &Fb : wf.core()) {
      if (E1.isZero(Fa, Fb))
        continue;
      const auto dv = rpa.dV_complex(Fa, Fb);
      REQUIRE(dv.real() != 0.0);
      REQUIRE(dv.imag() == 0.0);
      REQUIRE(rpa.dV(Fa, Fb) == dv.real());
    }
  }

  // Above the 2p threshold only. E1 channels of 2p- (j=1/2) are s1/2 and
  // d3/2 (kappa -1, 2); of 2p+ (j=3/2) are s1/2, d3/2, d5/2 (kappa -1, 2,
  // -3): five open channels
  const double omega = 1.5;
  rpa.solve_core(omega, 100, false);
  REQUIRE(rpa.last_eps() < rpa.eps_target());
  const auto channels = rpa.channel_list();
  const auto A = rpa.A_phys();
  REQUIRE(channels.size() == 5);
  REQUIRE(A.size() == 5);
  for (std::size_t i = 0; i < channels.size(); ++i) {
    const auto &Fa = wf.core().at(channels[i].i_core);
    REQUIRE(Fa.l() == 1);
    REQUIRE(channels[i].en == Approx(Fa.en() + omega));
    const auto kappa = channels[i].kappa;
    REQUIRE((kappa == -1 || kappa == 2 || kappa == -3));
    REQUIRE(std::isfinite(std::abs(A[i])));
    REQUIRE(std::abs(A[i]) > 0.0);
    // The (orbital, kappa) accessor returns the same amplitude
    REQUIRE(rpa.A_phys(Fa, kappa) == A[i]);
  }
  // Closed orbitals have no amplitude (1s: E1 channels p1/2, p3/2)
  REQUIRE(rpa.A_phys(*F1s, 1) == std::complex<double>{0.0, 0.0});
  REQUIRE(rpa.A_phys(*F1s, -2) == std::complex<double>{0.0, 0.0});

  // Above threshold the corrections are complex: some bound-bound dV has an
  // imaginary part
  bool any_imaginary = false;
  for (const auto &Fa : wf.core()) {
    for (const auto &Fb : wf.core()) {
      if (E1.isZero(Fa, Fb))
        continue;
      if (rpa.dV_complex(Fa, Fb).imag() != 0.0) {
        any_imaginary = true;
      }
    }
  }
  REQUIRE(any_imaginary);

  // Re-solve at the same omega: warm start from the converged solution, so
  // it converges at once, to the same amplitudes (to the convergence target)
  const auto its_first = rpa.last_its();
  REQUIRE(its_first >= 1);
  rpa.solve_core(omega, 100, false);
  REQUIRE(rpa.last_its() < its_first);
  const auto A_again = rpa.A_phys();
  REQUIRE(A_again.size() == A.size());
  for (std::size_t i = 0; i < A.size(); ++i) {
    REQUIRE(std::abs(A_again[i] - A[i]) < 1.0e-4 * std::abs(A[i]));
  }

  // clear() resets everything
  rpa.clear();
  REQUIRE(rpa.channel_list().empty());
  REQUIRE(rpa.A_phys().empty());
  REQUIRE(rpa.dV_complex(*F1s, *F2p) == std::complex<double>{0.0, 0.0});
}

//==============================================================================
//! Below every ionisation threshold no channel is open, so TDHFcntm must
//! reproduce bound TDHF: the same fixed point (dV matrix elements), with the
//! imaginary parts exactly zero. E1 (odd parity) static and at finite
//! frequency, and E2 (even parity: diagonal channels, de projection).
TEST_CASE("TDHFcntm: matches TDHF below threshold",
          "[ExternalField][TDHF][TDHFcntm][unit]") {

  Wavefunction wf({2000, 1.0e-6, 30.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-10, false);

  const auto E1 = DiracOperator::E1(wf.grid());
  const auto E2 = DiracOperator::Ek(wf.grid(), 2);

  struct Case {
    std::string name;
    const DiracOperator::TensorOperator *h;
    double omega;
  };
  const std::vector<Case> cases{
    {"E1", &E1, 0.0}, {"E1", &E1, 0.5}, {"E2", &E2, 0.0}};

  fmt::print("\nTDHFcomplex vs TDHF below threshold (Ne):\n");
  fmt::print("{:>4s} {:>5s} {:>4s} {:>4s} {:>14s} {:>14s} {:>9s} {:>9s}\n", "h",
             "omega", "a", "b", "TDHF", "complex", "Im", "rel.diff");

  for (const auto &[name, h, omega] : cases) {
    auto rpa_bound = ExternalField::TDHF(h, wf.vHF());
    auto rpa_c = ExternalField::TDHFcntm(h, wf.vHF());
    rpa_bound.eps_target() = 1.0e-14;
    rpa_c.eps_target() = 1.0e-14;
    rpa_bound.solve_core(omega, 200, false);
    rpa_c.solve_core(omega, 200, false);
    REQUIRE(rpa_c.channel_list().empty());

    // Scale for the comparison: the largest dV in this case
    double dv_max = 0.0;
    for (const auto &Fa : wf.core()) {
      for (const auto &Fb : wf.core()) {
        if (h->isZero(Fa, Fb))
          continue;
        dv_max = std::max(dv_max, std::abs(rpa_bound.dV(Fa, Fb)));
      }
    }
    REQUIRE(dv_max > 0.0);

    for (const auto &Fa : wf.core()) {
      for (const auto &Fb : wf.core()) {
        if (h->isZero(Fa, Fb))
          continue;
        const auto dv_b = rpa_bound.dV(Fa, Fb);
        const auto dv_c = rpa_c.dV_complex(Fa, Fb);
        const auto rel = std::abs(dv_c.real() - dv_b) / dv_max;
        fmt::print("{:>4s} {:5.2f} {:>4s} {:>4s} {:14.7e} {:14.7e} {:9.1e} "
                   "{:9.1e}\n",
                   name, omega, Fa.shortSymbol(), Fb.shortSymbol(), dv_b,
                   dv_c.real(), dv_c.imag(), rel);
        REQUIRE(dv_c.imag() == 0.0);
        REQUIRE(rel < 1.0e-6);
        REQUIRE(rpa_c.dV(Fa, Fb) == dv_c.real());
      }
    }
  }
}

//==============================================================================
//! Gauge invariance above threshold: the E1 ionisation amplitude in length
//! and velocity form (frequency dependent: t_+ at +omega, t_- at -omega)
//! must agree channel by channel as a COMPLEX number (same channel, same
//! phase reference), modulus and phase. Ne at omega = 4 au: 2s and 2p open.
TEST_CASE("TDHFcntm: E1 gauge invariance above threshold",
          "[ExternalField][TDHF][TDHFcntm][integration]") {

  Wavefunction wf({5000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-10, false);

  const double omega = 4.0;
  const auto E1 = DiracOperator::E1(wf.grid());
  const auto E1v = DiracOperator::E1v(wf.alpha(), omega);
  const auto E1v_minus = DiracOperator::E1v(wf.alpha(), -omega);

  auto rpa_L = ExternalField::TDHFcntm(&E1, wf.vHF());
  rpa_L.eps_target() = 1.0e-7;
  rpa_L.solve_core(omega, 60, false);
  const auto channels = rpa_L.channel_list();
  const auto A_L = rpa_L.A_phys();

  auto rpa_V = ExternalField::TDHFcntm(&E1v, wf.vHF(), &E1v_minus);
  rpa_V.eps_target() = 1.0e-7;
  rpa_V.solve_core(omega, 60, false);
  const auto A_V = rpa_V.A_phys();

  REQUIRE(channels.size() == 7);
  REQUIRE(A_L.size() == channels.size());
  REQUIRE(A_V.size() == channels.size());

  fmt::print("\nE1 gauge test (outgoing-wave TDHF), Ne omega = {:.1f} au, "
             "{} open channels:\n"
             "{:>6s} {:>6s} {:>11s} {:>11s} {:>9s} {:>9s} {:>9s}\n",
             omega, channels.size(), "shell", "kappa", "|A_L|", "|A_V|",
             "V/L - 1", "dphase", "|dA|/A");
  double sum_L = 0.0;
  double sum_V = 0.0;
  for (std::size_t i = 0; i < channels.size(); ++i) {
    const auto &Fa = wf.core().at(channels[i].i_core);
    const auto aL = std::abs(A_L[i]);
    const auto aV = std::abs(A_V[i]);
    const auto dphase = std::arg(A_V[i] / A_L[i]);
    const auto dA = std::abs(A_V[i] - A_L[i]) / aL;
    fmt::print("{:>6s} {:>6d} {:11.6f} {:11.6f} {:9.1e} {:9.1e} {:9.1e}\n",
               Fa.shortSymbol(), channels[i].kappa, aL, aV, aV / aL - 1.0,
               dphase, dA);
    REQUIRE(aV == Approx(aL).epsilon(1.0e-2));
    REQUIRE(dA < 1.0e-2);
    sum_L += aL * aL;
    sum_V += aV * aV;
  }
  // Total (cross-section level): tighter
  REQUIRE(sum_V == Approx(sum_L).epsilon(2.0e-3));
}
