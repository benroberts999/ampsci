#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFcntm.hpp"
#include "ExternalField/TDHFcomplex.hpp"
#include "HF/HartreeFock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <complex>

//==============================================================================
//! Below every ionisation threshold no channel is open, so the outgoing-wave
//! class has nothing to add: the imaginary sets must stay EXACTLY zero (no
//! source ever feeds them) and the real sets iterate to the bound TDHF fixed
//! point (dV matrix elements agree with bound TDHF).
TEST_CASE("cntmRPAcomplex: below threshold is bound TDHF",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  Wavefunction wf({4000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");

  const auto E1 = DiracOperator::E1(wf.grid());
  // 2p threshold is ~0.85 au: at omega = 0.5 all channels are closed
  const double omega = 0.5;

  auto rpa_bound = ExternalField::TDHF(&E1, wf.vHF());
  rpa_bound.solve_core(omega, 100, false);
  auto rpa_c = ExternalField::TDHFcomplex(&E1, wf.vHF());
  rpa_c.solve_core(omega, 100, false);

  REQUIRE(rpa_c.imaginary_norm2() == 0.0);
  REQUIRE(rpa_c.channel_list().empty());

  fmt::print("\nBound TDHF vs TDHFcomplex at omega = {:.2f} au (all closed):\n",
             omega);
  fmt::print("{:>9s} {:>13s} {:>13s} {:>13s}\n", "<a|dV|b>", "TDHF", "complex",
             "Im");
  for (const auto &Fa : wf.core()) {
    for (const auto &Fb : wf.core()) {
      if (E1.isZero(Fa, Fb))
        continue;
      const auto dv_b = rpa_bound.dV(Fa, Fb);
      const auto dv_c = rpa_c.dV_complex(Fa, Fb);
      fmt::print("{:>4s} {:>4s} {:13.6e} {:13.6e} {:13.1e}\n", Fa.shortSymbol(),
                 Fb.shortSymbol(), dv_b, dv_c.real(), dv_c.imag());
      REQUIRE(dv_c.imag() == 0.0);
      // Anderson outer driver: same fixed point, to its residual target
      REQUIRE(dv_c.real() == Approx(dv_b).epsilon(1.0e-3));
      REQUIRE(rpa_c.dV(Fa, Fb) == dv_c.real());
    }
  }
}

//==============================================================================
//! THE validation of the outgoing-wave method: its converged complex
//! amplitudes must equal the Johnson-route physical amplitudes of TDHFcntm
//! (standing-wave driven solve + on-shell K matrix from N_open seeded
//! solves), channel by channel, as COMPLEX numbers: modulus AND phase. The
//! two routes share nothing but the channel solves (no K matrix is built
//! here), so agreement confirms the outgoing Green's-function
//! decomposition, the conj(Y) bookkeeping of the complex map, the
//! exchange-dressed incident wave, and the A = conj(K_+) convention. Also
//! checks the complex K_+ = pi*D identity of the converged channels.
TEST_CASE("cntmRPAcomplex: outgoing amplitudes vs Johnson K matrix",
          "[ExternalField][TDHF][cntmrpa][integration]") {

  Wavefunction wf({4000, 1.0e-6, 30.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  const auto E1 = DiracOperator::E1(wf.grid());

  // Smooth region: only 2p open (5 channels)
  const double omega = 1.5;

  auto rpa_J = ExternalField::TDHFcntm(&E1, wf.vHF());
  rpa_J.eps_target() = 1.0e-7;
  rpa_J.solve_core(omega, 60, false);
  REQUIRE(rpa_J.kmatrix().has_value());
  const auto &kmat = *rpa_J.kmatrix();
  const auto A_J = rpa_J.A_phys();

  auto rpa_c = ExternalField::TDHFcomplex(&E1, wf.vHF());
  rpa_c.eps_target() = 1.0e-7;
  rpa_c.solve_core(omega, 60, false);
  const auto channels = rpa_c.channel_list();
  const auto A_c = rpa_c.A_phys();

  REQUIRE(channels.size() == 5);
  REQUIRE(A_J.size() == channels.size());
  REQUIRE(A_c.size() == channels.size());

  fmt::print("\nNe omega = {:.2f} au: Johnson (K matrix, asym {:.1e}) vs "
             "outgoing-wave amplitudes A/pi\n"
             "TDHFcomplex: its {} eps {:.1e} KpiD {:.1e} [{}]\n",
             omega, kmat.asymmetry, rpa_c.last_its(), rpa_c.last_eps(),
             rpa_c.KpiD_dev(), rpa_c.KpiD_worst_channel());
  fmt::print("{:>6s} {:>6s} {:>12s} {:>12s} {:>12s} {:>12s} {:>9s} {:>9s}\n",
             "shell", "kappa", "Re A_J", "Im A_J", "Re A_c", "Im A_c", "|dA|/A",
             "|dA*|/A");
  double A_max = 0.0;
  for (const auto &a : A_J) {
    A_max = std::max(A_max, std::abs(a));
  }
  double worst = 0.0;
  for (std::size_t i = 0; i < channels.size(); ++i) {
    REQUIRE(channels[i].i_core == kmat.channels[i].i_core);
    REQUIRE(channels[i].kappa == kmat.channels[i].kappa);
    const auto &Fa = wf.core()[channels[i].i_core];
    const auto dA = std::abs(A_c[i] - A_J[i]) / A_max;
    // the conjugate: distinguishes the phase convention (must NOT agree)
    const auto dAc = std::abs(A_c[i] - std::conj(A_J[i])) / A_max;
    fmt::print("{:>6s} {:>6d} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:9.1e} "
               "{:9.1e}\n",
               Fa.shortSymbol(), channels[i].kappa, A_J[i].real(),
               A_J[i].imag(), A_c[i].real(), A_c[i].imag(), dA, dAc);
    worst = std::max(worst, dA);
    REQUIRE(std::isfinite(std::abs(A_c[i])));
    // The per-channel entry of the (Fa, kappa) accessor is the same number
    REQUIRE(rpa_c.A_phys(Fa, channels[i].kappa) == A_c[i]);
    REQUIRE(rpa_c.D_phys(Fa, channels[i].kappa) == std::abs(A_c[i]));
  }
  // Modulus and phase agree to the numerical floor of the two routes (the
  // K-matrix route carries the Kbar asymmetry ~1e-3; the SCF targets 1e-7)
  REQUIRE(worst < 5.0e-3);
  // and the phase is NOT the conjugate: the dominant channels have a
  // non-negligible imaginary part, so the conjugate must differ clearly
  {
    std::size_t i_dom = 0;
    for (std::size_t i = 0; i < channels.size(); ++i) {
      if (std::abs(A_J[i]) > std::abs(A_J[i_dom]))
        i_dom = i;
    }
    REQUIRE(std::abs(A_J[i_dom].imag()) > 0.02 * A_max);
    REQUIRE(std::abs(A_c[i_dom] - std::conj(A_J[i_dom])) >
            10.0 * std::abs(A_c[i_dom] - A_J[i_dom]));
  }

  // Complex K_+ = pi*D identity of the converged solve
  REQUIRE(rpa_c.KpiD_dev() < 0.02);
  fmt::print("{:>6s} {:>6s} {:>12s} {:>12s} {:>12s} {:>12s}\n", "shell",
             "kappa", "Re K", "Im K", "Re piD", "Im piD");
  for (const auto &Fa : wf.core()) {
    for (const auto &oc : rpa_c.outgoing_channels(Fa)) {
      fmt::print("{:>6s} {:>6d} {:12.5e} {:12.5e} {:12.5e} {:12.5e}\n",
                 Fa.shortSymbol(), oc.kappa, oc.K.real(), oc.K.imag(),
                 M_PI * oc.D.real(), M_PI * oc.D.imag());
    }
  }

  // The TDHF-like route: the dressed matrix element with the V^{N-1} HF
  // continuum bra (solved exactly as for the bare/standing columns:
  // hole_particle and force_orthog) is the physical amplitude too,
  // |<Fe|t|a> + dV_complex(Fe, Fa)| = |A|, with the same Fe-quality limit
  // as those columns (its inner region comes from the local-exchange
  // trick; largest for the s channels). Modulus compared; the phase sits
  // in the HF reference (K_+ exp(i delta_ex)), not the KS one.
  fmt::print("\nFe route: |<Fe|t|a> + dV_complex(Fe, Fa)| vs |A_phys|\n"
             "{:>6s} {:>6s} {:>12s} {:>12s} {:>9s}\n",
             "shell", "kappa", "|D|", "|A|", "D/A - 1");
  double dom_A = 0.0, dom_ratio = 0.0;
  for (const auto &Fa : wf.core()) {
    const auto ec = omega + Fa.en();
    if (ec <= 0.0)
      continue;
    ContinuumOrbitals bra(wf.vHF());
    bra.solveContinuumHF(ec, std::max(Fa.l() - 1, 0), Fa.l() + 1, &Fa, false,
                         true, true);
    for (const auto &Fe : bra.orbitals) {
      if (E1.isZero(Fe, Fa) || Fe.norm2() == 0.0)
        continue;
      const auto D = E1.reducedME(Fe, Fa) + rpa_c.dV_complex(Fe, Fa);
      const auto A = rpa_c.A_phys(Fa, Fe.kappa());
      const auto ratio = std::abs(D) / std::abs(A);
      fmt::print("{:>6s} {:>6d} {:12.5e} {:12.5e} {:9.1e}\n", Fa.shortSymbol(),
                 Fe.kappa(), std::abs(D), std::abs(A), ratio - 1.0);
      if (std::abs(A) > dom_A) {
        dom_A = std::abs(A);
        dom_ratio = ratio;
      }
      // Fe-quality limited (few % for p, ~10% for s channels)
      REQUIRE(ratio == Approx(1.0).epsilon(0.15));
    }
  }
  // Dominant (d) channel: Fe is accurate there
  REQUIRE(dom_ratio == Approx(1.0).epsilon(0.03));
}

//==============================================================================
//! Gauge invariance with the outgoing-wave class: E1 length vs velocity
//! (frequency-dependent, t_+ at +omega and t_- at -omega). The physical
//! amplitude is gauge invariant as a COMPLEX number (same channel, same
//! phase reference), so both modulus and phase are compared. Also exercises
//! set_operator re-use of the outgoing caches within the (rank, parity)
//! class (velocity solved on a copy of the length instance).
TEST_CASE("cntmRPAcomplex: gauge invariance (E1 vs E1v)",
          "[ExternalField][TDHF][cntmrpa][integration]") {

  Wavefunction wf({5000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  const auto E1 = DiracOperator::E1(wf.grid());

  // 2s + 2p open
  const double omega = 4.0;
  auto E1v = DiracOperator::E1v(wf.alpha(), omega);
  auto E1v_minus = DiracOperator::E1v(wf.alpha(), -omega);

  auto rpa_L = ExternalField::TDHFcomplex(&E1, wf.vHF());
  rpa_L.eps_target() = 1.0e-7;
  rpa_L.solve_core(omega, 60, false);
  const auto channels = rpa_L.channel_list();
  const auto A_L = rpa_L.A_phys();

  // Velocity on a copy of the (prepared, solved) length instance: keeps the
  // channel pairs and outgoing caches, restarts the corrections
  auto rpa_V = rpa_L;
  rpa_V.set_operator(&E1v, &E1v_minus);
  REQUIRE(rpa_V.channel_list().size() == channels.size());
  rpa_V.solve_core(omega, 60, false);
  const auto A_V = rpa_V.A_phys();

  REQUIRE(!channels.empty());
  REQUIRE(A_L.size() == channels.size());
  REQUIRE(A_V.size() == channels.size());

  fmt::print("\nGauge test (outgoing-wave) at omega = {:.1f} au ({} open "
             "channels):\n{:>6s} {:>6s} {:>11s} {:>11s} {:>9s} {:>9s} {:>9s}\n",
             omega, channels.size(), "shell", "kappa", "|A_L|/pi", "|A_V|/pi",
             "V/L - 1", "dphase", "|dA|/A");
  double sum_L = 0.0, sum_V = 0.0;
  for (std::size_t i = 0; i < channels.size(); ++i) {
    const auto &Fa = wf.core()[channels[i].i_core];
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

//==============================================================================
//! The batched dV builder of TDHFcomplex (a copy of TDHFcntm's, applied to
//! explicit correction sets) must reproduce the per-task TDHF::dV_rhs
//! exactly. Checked on a mid-convergence state above threshold (nonzero X
//! and Y, open and closed channels), for both task types (conj).
TEST_CASE("cntmRPAcomplex: batched dV_rhs (dV_rhs_all vs dV_rhs)",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  Wavefunction wf({4000, 1.0e-6, 30.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  const auto E1 = DiracOperator::E1(wf.grid());

  auto rpa = ExternalField::TDHFcomplex(&E1, wf.vHF());
  rpa.solve_core(1.5, 3, false);
  REQUIRE(rpa.imaginary_norm2() > 0.0);

  // The real sets (nonzero above threshold), via the public accessors
  std::vector<std::vector<DiracSpinor>> X, Y;
  for (const auto &Fc : wf.core()) {
    X.push_back(rpa.get_dPsis(Fc, ExternalField::dPsiType::X));
    Y.push_back(rpa.get_dPsis(Fc, ExternalField::dPsiType::Y));
  }
  const auto batched = rpa.dV_rhs_all(X, Y, true);

  double worst = 0.0;
  for (std::size_t ib = 0; ib < wf.core().size(); ++ib) {
    const auto &Fa = wf.core()[ib];
    for (std::size_t be = 0; be < X[ib].size(); ++be) {
      const auto kappa = X[ib][be].kappa();
      for (const bool conj : {false, true}) {
        const auto per_task = rpa.dV_rhs(kappa, Fa, conj);
        const auto &b = conj ? batched.Y[ib][be] : batched.X[ib][be];
        const auto diff = (b - per_task).norm2();
        const auto scale = per_task.norm2();
        const auto rel =
          scale == 0.0 ? std::sqrt(diff) : std::sqrt(diff / scale);
        worst = std::max(worst, rel);
        REQUIRE(rel < 1.0e-10);
      }
    }
  }
  fmt::print("\nTDHFcomplex batched dV_rhs vs per-task: worst rel. diff "
             "{:.1e}\n",
             worst);
}
