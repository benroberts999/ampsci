#include "DiracOperator/include.hpp"
#include "ExternalField/TDHFcntm.hpp"
#include "HF/HartreeFock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>

//==============================================================================
//! Below all ionisation thresholds every channel is bound, so TDHFcntm
//! (Anderson bound solves) and bound TDHF (damped solves) iterate to the SAME
//! TDHF fixed point: the dV matrix elements must agree. Verifies both that
//! the continuum class reproduces the bound physics, and that the bound TDHF
//! class is untouched by the refactor.
TEST_CASE("cntmRPA: below-threshold matches bound TDHF",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  Wavefunction wf({4000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");

  const auto E1 = DiracOperator::E1(wf.grid());
  // 2p threshold is ~0.85 au: at omega = 0.5 all channels are closed
  const double omega = 0.5;

  auto rpa_bound = ExternalField::TDHF(&E1, wf.vHF());
  rpa_bound.solve_core(omega, 100, false);
  // damped driver: same outer iteration as bound TDHF (tight agreement)
  auto rpa_cntm = ExternalField::TDHFcntm(&E1, wf.vHF());
  rpa_cntm.set_anderson(false);
  rpa_cntm.solve_core(omega, 100, false);
  // Anderson (default) driver: same fixed point, different iteration path;
  // its residual target maps to a slightly looser effective tolerance
  auto rpa_and = ExternalField::TDHFcntm(&E1, wf.vHF());
  rpa_and.solve_core(omega, 100, false);

  fmt::print("\nBound TDHF vs TDHFcntm at omega = {:.2f} au (all closed):\n",
             omega);
  fmt::print("{:>9s} {:>13s} {:>13s} {:>13s}\n", "<a|dV|b>", "TDHF", "damped",
             "Anderson");
  for (const auto &Fa : wf.core()) {
    for (const auto &Fb : wf.core()) {
      if (E1.isZero(Fa, Fb))
        continue;
      const auto dv_b = rpa_bound.dV(Fa, Fb);
      const auto dv_c = rpa_cntm.dV(Fa, Fb);
      const auto dv_a = rpa_and.dV(Fa, Fb);
      fmt::print("{:>4s} {:>4s} {:13.6e} {:13.6e} {:13.6e}\n", Fa.shortSymbol(),
                 Fb.shortSymbol(), dv_b, dv_c, dv_a);
      // Damped driver: same outer scheme as bound TDHF (only the inner
      // mixed-states solver differs: damped vs Anderson) -- tight:
      REQUIRE(dv_c == Approx(dv_b).epsilon(1.0e-5));
      // Anderson outer driver: same fixed point, to its residual target:
      REQUIRE(dv_a == Approx(dv_b).epsilon(1.0e-3));
    }
  }
}

//==============================================================================
//! Internal consistency of the continuum solve: the standing-wave (K-matrix)
//! amplitude of the converged correction must equal pi * D, with
//! D = <F_reg|(t + dV')phi_a> the dressed amplitude built from the SAME
//! (local-KS) F_reg (Methods eq. "c-K-overlap"). Stringent check that the
//! solver applies the boundary condition consistently AND that the source is
//! short-ranged (a long-range source residual breaks the surface-term
//! derivation of the relation). The tolerance is physics-limited: the
//! iterated source carries the full nonlocal exchange while F_reg is a
//! local-KS eigenstate, plus the O(1/(p r_box)) irregular-seed slip.
TEST_CASE("cntmRPA: K = pi*D consistency",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  Wavefunction wf({10000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");

  const auto E1 = DiracOperator::E1(wf.grid());

  for (const auto omega : {1.7, 3.0}) {
    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-6;
    rpa.set_unitarise(false); // standing-wave quantities only
    rpa.solve_core(omega, 60, false);

    fmt::print("\nK = pi*D consistency, omega = {:.2f} au:\n", omega);
    fmt::print("{:>4s} {:>6s} {:>12s} {:>12s} {:>8s}\n", "a", "kappa", "K",
               "pi*D", "K/piD-1");
    int n_channels = 0;
    for (const auto &Fb : wf.core()) {
      for (const auto &ch : rpa.open_channels(Fb)) {
        if (ch.K == 0.0)
          continue; // grid too sparse to resolve this channel
        const auto piD = M_PI * ch.D;
        fmt::print("{:>4s} {:>6d} {:12.5e} {:12.5e} {:8.1e}\n",
                   Fb.shortSymbol(), ch.kappa, ch.K, piD, ch.K / piD - 1.0);
        REQUIRE(ch.K == Approx(piD).epsilon(5.0e-2));
        ++n_channels;
      }
    }
    REQUIRE(n_channels > 0);
  }
}

//==============================================================================
//! High-frequency limit: far above the ionisation thresholds the core has no
//! time to respond (the response corrections fall off as ~1/omega), so the
//! PHYSICAL (unitarised) RPA cross section must approach the bare (V^{N-1}
//! tree) one. nb: this limit belongs to the unitarised amplitudes |A|: the
//! standing-wave (principal-value) amplitudes differ from |A| by the
//! on-shell phase rotation, and their per-channel ratios to bare do NOT
//! tend to 1 (the occupied components of the corrections carry a
//! bare-scaled contribution through the orthogonalised-bra matrix
//! elements). The approach is slow (~1/omega), so the test checks the
//! TREND -- |sigma_U/sigma_bare - 1| decreasing above the last threshold --
//! plus a loose bound at the largest omega.
TEST_CASE("cntmRPA: high-omega limit", "[ExternalField][TDHF][cntmrpa][unit]") {

  Wavefunction wf({5000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  const auto E1 = DiracOperator::E1(wf.grid());

  // Ne thresholds: 2p ~0.85, 2s ~1.93, 1s ~32.8 au. All shells are open at
  // these omega (the 15 -> 60 au step crosses the 1s threshold, so monotone
  // decay is only required for the 60 -> 240 au step):
  const auto omegas = std::vector{15.0, 60.0, 240.0};

  // |sigma_U/sigma_bare - 1| per omega
  std::vector<double> dev;

  for (const auto omega : omegas) {
    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-6;
    rpa.solve_core(omega, 60, false); // unitarised (K matrix built)

    // Bare cross-section from the (orthogonalised) V^{N-1} HF bra states
    double sum_D0_2 = 0.0;
    for (const auto &Fb : wf.core()) {
      const auto ec = Fb.en() + omega;
      if (ec <= 0.0)
        continue;
      ContinuumOrbitals cntm(wf.vHF());
      cntm.solveContinuumHF(ec, std::max(Fb.l() - 1, 0), Fb.l() + 1, &Fb, false,
                            true, true);
      for (const auto &Fe : cntm.orbitals) {
        if (E1.isZero(Fe, Fb) || Fe.norm2() == 0.0)
          continue;
        const auto D0 = E1.reducedME(Fe, Fb);
        sum_D0_2 += D0 * D0;
      }
    }
    double sum_A_2 = 0.0;
    for (const auto &Fb : wf.core()) {
      for (const auto &oc : rpa.open_channels(Fb)) {
        const auto Dp = rpa.D_phys(Fb, oc.kappa);
        sum_A_2 += Dp * Dp;
      }
    }
    REQUIRE(sum_D0_2 > 0.0);
    REQUIRE(sum_A_2 > 0.0);

    const auto ratio = sum_A_2 / sum_D0_2;
    fmt::print("High-omega limit, omega = {:>3.0f} au: sigma_U/sigma_bare = "
               "{:.4f}\n",
               omega, ratio);
    dev.push_back(std::abs(ratio - 1.0));
  }

  // Monotone decay above the last (1s) threshold, and small at the top:
  REQUIRE(dev[2] < dev[1]);
  REQUIRE(dev[2] < 0.1);
}

//==============================================================================
//! The hole-particle acceptance test. (a) The compensated source
//! S = [(t + dV')phi_a]_beta for an ionised orbital must be SHORT-RANGED:
//! the +y^0_aa*chi hole compensation must cancel, pointwise, the 1/r tail
//! hidden in the b=a term of dV*phi_a (the self-interaction moved into the
//! static V^{N-1} Hamiltonian). A residual tail makes the response
//! box-dependent (log(r_box)-divergent matrix elements) and the RPA
//! correction spuriously large. (b) The converged dressed amplitudes must
//! then be box-independent.
TEST_CASE("cntmRPA: source short-rangedness, box-independence",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  const double omega = 1.7;

  // --- (a) tail fraction of the source, with and without compensation
  {
    Wavefunction wf({10000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                    {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Ne]");
    const auto E1 = DiracOperator::E1(wf.grid());

    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-6;
    rpa.set_unitarise(false); // the source and standing D only
    rpa.solve_core(omega, 60, false);

    const auto &gr = wf.grid();
    // integral of |f| over [r0, rmax] (relative measures only: du omitted)
    const auto int_abs_f = [&gr](const DiracSpinor &F, double r0) {
      double sum = 0.0;
      for (std::size_t i = 0; i < F.max_pt(); ++i) {
        if (gr.r(i) >= r0) {
          sum += std::abs(F.f(i)) * gr.drdu(i);
        }
      }
      return sum;
    };
    const auto r_tail = 0.7 * gr.r().back();

    fmt::print("\nSource tail fraction (r > {:.0f} au), omega = {:.2f} au:\n",
               r_tail, omega);
    fmt::print("{:>4s} {:>6s} {:>13s} {:>13s}\n", "a", "kappa", "compensated",
               "bare dV");
    for (const auto &Fb : wf.core()) {
      if (Fb.en() + omega <= 0.0)
        continue; // closed at this omega
      for (const auto &chi : rpa.get_dPsis(Fb, ExternalField::dPsiType::X)) {
        const auto kappa = chi.kappa();
        const auto hFb = E1.reduced_rhs(kappa, Fb);
        const auto S_comp = hFb + rpa.dVprime_rhs(kappa, Fb, false, chi);
        const auto S_bare = hFb + rpa.dV_rhs(kappa, Fb, false);
        const auto tail_comp =
          int_abs_f(S_comp, r_tail) / int_abs_f(S_comp, 0.0);
        const auto tail_bare =
          int_abs_f(S_bare, r_tail) / int_abs_f(S_bare, 0.0);
        fmt::print("{:>4s} {:>6d} {:13.2e} {:13.2e}\n", Fb.shortSymbol(), kappa,
                   tail_comp, tail_bare);
        // The compensated source is short-ranged (tail negligible), and the
        // compensation is doing real work (>= 10x smaller than bare):
        REQUIRE(tail_comp < 1.0e-3);
        REQUIRE(10.0 * tail_comp < tail_bare);
      }
    }
  }

  // --- (b) box-independence of the dressed production amplitude
  {
    // Two boxes, comparable point density at large r
    const auto grids = {std::pair{10000ul, 40.0}, std::pair{18750ul, 75.0}};
    // D indexed by [orbital symbol + channel kappa]
    std::vector<std::vector<double>> Ds;
    std::vector<std::vector<std::string>> labels;
    for (const auto &[npts, rmax] : grids) {
      Wavefunction wf({npts, 1.0e-6, rmax, 1.0, "loglinear", -1.0},
                      {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
      wf.solve_core("HartreeFock", "[Ne]");
      const auto E1 = DiracOperator::E1(wf.grid());
      auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
      rpa.eps_target() = 1.0e-6;
      rpa.set_unitarise(false);
      rpa.solve_core(omega, 60, false);

      std::vector<double> D_box;
      std::vector<std::string> lab_box;
      for (const auto &Fb : wf.core()) {
        const auto ec = Fb.en() + omega;
        if (ec <= 0.0)
          continue;
        ContinuumOrbitals cntm(wf.vHF());
        cntm.solveContinuumHF(ec, std::max(Fb.l() - 1, 0), Fb.l() + 1, &Fb,
                              false, true, true);
        for (const auto &Fe : cntm.orbitals) {
          if (E1.isZero(Fe, Fb) || Fe.norm2() == 0.0)
            continue;
          D_box.push_back(E1.reducedME(Fe, Fb) + rpa.dV_cntm(Fe, Fb));
          lab_box.push_back(Fb.shortSymbol() + "->" + Fe.shortSymbol());
        }
      }
      Ds.push_back(D_box);
      labels.push_back(lab_box);
    }

    fmt::print("\nBox-independence of dressed D, omega = {:.2f} au:\n", omega);
    fmt::print("{:>11s} {:>13s} {:>13s} {:>8s}\n", "channel", "40 au", "75 au",
               "eps");
    REQUIRE(Ds[0].size() == Ds[1].size());
    for (std::size_t i = 0; i < Ds[0].size(); ++i) {
      REQUIRE(labels[0][i] == labels[1][i]);
      const auto eps = std::abs(Ds[1][i] / Ds[0][i] - 1.0);
      fmt::print("{:>11s} {:13.6e} {:13.6e} {:8.1e}\n", labels[0][i], Ds[0][i],
                 Ds[1][i], eps);
      // Energy-normalised continuum bra and short-ranged source: the
      // physical amplitude must not depend on the box size:
      REQUIRE(Ds[1][i] == Approx(Ds[0][i]).epsilon(1.0e-2));
    }
  }
}

//==============================================================================
//! Validation against Johnson & Cheng, PRA 20, 978 (1979), Fig 1: Ne total
//! E1 (length) photoionisation cross-section, RRPA vs HF. Expected values
//! are read off the published figure (so carry ~10-15% reading uncertainty).
//! The comparable quantity is the UNITARISED (physical) cross-section --
//! Johnson's RRPA amplitudes are the outgoing-wave (eigenchannel) ones.
//! Chosen omega sit in the smooth region above the 2s threshold, away from
//! the 2s->np autoionising resonances (cf. Johnson & Cheng: "we avoid these
//! resonances altogether").
//! nb: tolerance 25%: our (gauge-invariant, L=V) result sits systematically
//! ~15-20% above the figure-read values in this 60-120 eV window (while
//! matching experiment near the 2p maximum and at high omega) -- open
//! question, tracked in the status notes.
TEST_CASE("cntmRPA: Johnson 1979 Ne cross-section",
          "[ExternalField][TDHF][cntmrpa][integration]") {

  Wavefunction wf({10000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  const auto E1 = DiracOperator::E1(wf.grid());

  // {omega (au), RRPA sigma (Mb) read from Fig 1 of Johnson & Cheng 1979}
  const auto expected = std::vector{std::pair{2.31, 6.0}, std::pair{2.91, 4.7},
                                    std::pair{3.67, 3.5}};

  fmt::print("\nNe total E1 cross-section vs Johnson & Cheng (1979):\n");
  fmt::print("{:>7s} {:>10s} {:>10s} {:>12s} {:>9s}\n", "w (au)", "bare(Mb)",
             "RPA-U(Mb)", "Johnson(Mb)", "eps");

  for (const auto &[omega, sigma_expct] : expected) {
    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-5;
    rpa.solve_core(omega, 90, false); // unitarised (K matrix built)
    const auto &kmat = *rpa.kmatrix();

    // sigma = (4 pi^2 alpha / 3) * omega * sum_channels |D|^2, in Mb
    const auto Ksigma = 4.0 * M_PI * M_PI * PhysConst::alpha * omega / 3.0 *
                        PhysConst::aB_cm * PhysConst::aB_cm * 1.0e18;
    double sigma_0 = 0.0;
    for (const auto &Fb : wf.core()) {
      const auto ec = Fb.en() + omega;
      if (ec <= 0.0)
        continue;
      ContinuumOrbitals cntm(wf.vHF());
      cntm.solveContinuumHF(ec, std::max(Fb.l() - 1, 0), Fb.l() + 1, &Fb, false,
                            true, true);
      for (const auto &Fe : cntm.orbitals) {
        if (E1.isZero(Fe, Fb) || Fe.norm2() == 0.0)
          continue;
        const auto D0 = E1.reducedME(Fe, Fb);
        sigma_0 += Ksigma * D0 * D0;
      }
    }
    double sigma_U = 0.0;
    for (const auto &Fb : wf.core()) {
      for (const auto &oc : rpa.open_channels(Fb)) {
        const auto Dp = rpa.D_phys(Fb, oc.kappa);
        sigma_U += Ksigma * Dp * Dp;
      }
    }

    fmt::print("{:7.2f} {:10.2f} {:10.2f} {:12.1f} {:9.1e}\n", omega, sigma_0,
               sigma_U, sigma_expct, rpa.last_eps());

    // TDHF converged (away from the resonances):
    REQUIRE(rpa.last_eps() < 1.0e-2);
    // Kbar consistency (reciprocity, not imposed):
    REQUIRE(kmat.asymmetry < 0.05);
    // Agreement with the published RRPA curve (see @note above):
    REQUIRE(sigma_U == Approx(sigma_expct).epsilon(0.25));
  }
}

//==============================================================================
//! Continuum TDHF demonstration: E1 photoionisation matrix elements for Ne,
//! across the available approximation levels and a scan of frequencies.
/*!
  For each open E1 channel (en_c = en_a + omega > 0), the reduced matrix element
  D = <ec,a' || E1 (+dV') || a> is printed for six methods:

    - noself : old continuum-HF bra, NO self-interaction subtracted (plain V^N)
    - V^N-1  : old continuum-HF bra, one-electron self-interaction subtracted
               (the V^{N-1} residual-ion tree level)
    - Green0 : new continuum TDHF, standing-wave forward solve, dV switched off
               (= V^{N-1} tree; checks the new machinery reproduces it)
    - 1-iter : new continuum TDHF, one iteration of dV (first-order core pol.)
    - closed : new continuum TDHF, converged, open channels suppressed
               (only the bound part of the core responds)
    - full   : new continuum TDHF, converged, full dV

  Scanned over a range of omega from just above the 2p threshold up to very high
  energy (where 2s and 1s also open). Channels the radial grid is too sparse to
  resolve at very high ec are skipped. The full RPA converges at all these
  frequencies except near a collective resonance (e.g. omega ~ 1.5 au), where
  the OUTER dV self-consistency converges only slowly (it stays bounded, but
  the eps floor is higher) -- a genuine RPA-level effect, not a solver failure.

  Not intended to be high-accuracy; just demonstrates the methods run and the
  RPA moves the matrix elements.
*/
TEST_CASE("cntmRPA: Ne E1 photoionisation table",
          "[ExternalField][TDHF][cntmrpa][integration]") {

  Wavefunction wf({10000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  wf.printCore();

  const auto E1 = DiracOperator::E1(wf.grid());

  for (const auto omega : {1.5, 1.7, 3.0, 10.0, 30.0}) {

    // New continuum TDHF, at the various levels (standing-wave matrix
    // elements only: no K matrix):
    auto rpa_g0 = ExternalField::TDHFcntm(&E1, wf.vHF()); // fwd solve, no dV
    rpa_g0.set_unitarise(false);
    rpa_g0.solve_core(omega, 0);                           // 0 iterations
    auto rpa_1it = ExternalField::TDHFcntm(&E1, wf.vHF()); // one dV iteration
    rpa_1it.set_unitarise(false);
    rpa_1it.solve_core(omega, 1);
    auto rpa_closed = ExternalField::TDHFcntm(&E1, wf.vHF()); // open chnls off
    rpa_closed.set_unitarise(false);
    rpa_closed.set_suppress_open(true);
    rpa_closed.solve_core(omega, 60);
    auto rpa_full = ExternalField::TDHFcntm(&E1, wf.vHF()); // full dV
    rpa_full.set_unitarise(false);
    // F_reg is re-normalised from its large-r envelope (~1e-6 noise) each
    // rebuild; that sets the achievable TDHF floor (well below the MEs).
    rpa_full.eps_target() = 1.0e-5;
    rpa_full.solve_core(omega, 60);

    fmt::print(
      "\nNe E1 photoionisation, omega = {:.2f} au  (TDHF eps = {:.1e})\n",
      omega, rpa_full.last_eps());
    fmt::print("D = <ec,a' || E1 (+dV') || a>\n");
    fmt::print(
      "{:>9s} {:>7s} {:>11s} {:>11s} {:>11s} {:>11s} {:>11s} {:>11s}\n",
      "a -> a'", "ec", "noself", "V^N-1", "Green0", "1-iter", "closed", "full");

    int n_channels = 0;
    for (const auto &Fb : wf.core()) {
      const auto ec = Fb.en() + omega;
      if (ec <= 0.0)
        continue; // closed shell at this omega

      // Continuum bra: WITHOUT self subtraction (plain V^N), and WITH (V^{N-1})
      ContinuumOrbitals cntm_bare(wf.vHF());
      cntm_bare.solveContinuumHF(ec, std::max(Fb.l() - 1, 0), Fb.l() + 1, &Fb,
                                 false, false, true);
      ContinuumOrbitals cntm_vn1(wf.vHF());
      cntm_vn1.solveContinuumHF(ec, std::max(Fb.l() - 1, 0), Fb.l() + 1, &Fb,
                                false, true, true);

      for (std::size_t i = 0; i < cntm_vn1.orbitals.size(); ++i) {
        const auto &Fe_bare = cntm_bare.orbitals.at(i);
        const auto &Fe = cntm_vn1.orbitals.at(i); // V^{N-1} state (RPA bra)
        if (E1.isZero(Fe, Fb))
          continue; // selection-rule-forbidden channel
        // Skip channels the grid is too sparse to resolve at this ec (the
        // continuum solver returns a zero orbital in that case):
        if (Fe.norm2() == 0.0 || Fe_bare.norm2() == 0.0)
          continue;

        const auto D_noself = E1.reducedME(Fe_bare, Fb);
        const auto D_vn1 = E1.reducedME(Fe, Fb); // tree, V^{N-1}
        const auto D_g0 = D_vn1 + rpa_g0.dV(Fe, Fb);
        const auto D_1it = D_vn1 + rpa_1it.dV(Fe, Fb);
        const auto D_closed = D_vn1 + rpa_closed.dV(Fe, Fb);
        const auto D_full = D_vn1 + rpa_full.dV(Fe, Fb);

        fmt::print("{:>4s} -> {:<3s} {:7.3f} {:11.4e} {:11.4e} {:11.4e} "
                   "{:11.4e} {:11.4e} {:11.4e}\n",
                   Fb.shortSymbol(), Fe.shortSymbol(), ec, D_noself, D_vn1,
                   D_g0, D_1it, D_closed, D_full);
        ++n_channels;

        // Tree level (V^{N-1}) is finite and nonzero:
        REQUIRE(std::isfinite(D_vn1));
        REQUIRE(std::abs(D_vn1) > 0.0);
        // The new continuum machinery (forward solve, dV off) reproduces the
        // V^{N-1} continuum-HF tree exactly:
        REQUIRE(std::abs(D_g0 - D_vn1) < 1.0e-8 * std::abs(D_vn1));
        // Subtracting the one-electron self-interaction (V^{N-1}) changes the
        // tree matrix element (the box-dependence fix matters):
        REQUIRE(std::abs(D_noself - D_vn1) > 1.0e-3 * std::abs(D_vn1));
        // The full RPA matrix element is finite (bounded) and shifted from tree:
        REQUIRE(std::isfinite(D_full));
        REQUIRE(std::abs(D_full - D_vn1) > 1.0e-6 * std::abs(D_vn1));
      }
    }

    // Some open E1 channels were found and processed at this omega:
    REQUIRE(n_channels > 0);
    // The TDHF stays bounded (no divergence); it converges tightly except near a
    // collective resonance, where it stalls at a higher eps (still finite).
    REQUIRE(std::isfinite(rpa_full.last_eps()));
  }
}

//==============================================================================
//! Unitarisation (on-shell rescattering, Johnson 1979 appendix): the
//! standing-wave RRPA amplitudes D have real poles where an eigenphase passes
//! through pi/2 (e.g. just above the 1s threshold); the physical amplitudes
//! A = (1 - i*Kbar)^{-1} pi*D must stay finite there. Checks: (a) Kbar is
//! symmetric (reciprocity -- not imposed anywhere by the construction);
//! (b) for the dominant channel the multichannel result reduces to the
//! single-channel form |A| ~ pi*D*cos(atan Kbar_ii); (c) through the Ne 1s
//! near-edge pole the standing-wave cross-section blows up while the
//! unitarised one stays at the physical (Henke-scale) value.
TEST_CASE("cntmRPA: unitarisation (rescattering)",
          "[ExternalField][TDHF][cntmrpa][integration]") {

  Wavefunction wf({4000, 1.0e-6, 30.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  const auto E1 = DiracOperator::E1(wf.grid());

  { // (a) + (b): consistency in the smooth region (only 2p open)
    const double omega = 1.5;
    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-7;
    rpa.solve_core(omega, 60, false);
    const auto &kmat = *rpa.kmatrix();

    fmt::print("\nRescattering at omega = {:.2f} au ({} open channels), "
               "Kbar asymmetry = {:.1e}\n",
               omega, kmat.channels.size(), kmat.asymmetry);
    fmt::print("{:>7s} {:>12s} {:>12s} {:>12s} {:>9s}\n", "kappa", "D",
               "D_phys", "D*cos(th)", "Kbar_ii");
    REQUIRE(kmat.channels.size() == 5);
    REQUIRE(kmat.asymmetry < 0.05);

    // Standing D = K/pi and physical |A|/pi per channel; open_channels()
    // iterates in the same (core, channel) order as kmat.channels
    std::vector<double> D, D_phys;
    for (const auto &Fb : wf.core()) {
      for (const auto &oc : rpa.open_channels(Fb)) {
        D.push_back(oc.K / M_PI);
        D_phys.push_back(rpa.D_phys(Fb, oc.kappa));
      }
    }
    REQUIRE(D.size() == kmat.channels.size());

    std::size_t i_dom = 0;
    for (std::size_t i = 0; i < kmat.channels.size(); ++i) {
      REQUIRE(std::isfinite(D_phys[i]));
      if (std::abs(D[i]) > std::abs(D[i_dom]))
        i_dom = i;
    }
    for (std::size_t i = 0; i < kmat.channels.size(); ++i) {
      const auto d_cos = std::abs(D[i]) * std::cos(std::atan(kmat.Kbar(i, i)));
      fmt::print("{:>7d} {:12.5e} {:12.5e} {:12.5e} {:9.5f}\n",
                 kmat.channels[i].kappa, D[i], D_phys[i], d_cos,
                 kmat.Kbar(i, i));
    }
    // Dominant channel: single-channel unitarisation dominates (off-diagonal
    // rescattering gives only a small correction here):
    const auto d_cos_dom =
      std::abs(D[i_dom]) * std::cos(std::atan(kmat.Kbar(i_dom, i_dom)));
    REQUIRE(D_phys[i_dom] == Approx(d_cos_dom).epsilon(0.10));
  }

  { // (c): the Ne 1s near-edge region. Just above the 1s threshold (893 eV)
    // an RRPA eigenphase sweeps through pi/2: the standing-wave amplitudes
    // pass through a real pole (their value there is arbitrary-large and
    // grid-sensitive), while the physical (unitarised) amplitude must stay
    // bounded at the smooth (Henke ~0.35 Mb) scale. Assert (i) standing and
    // physical amplitudes differ strongly in the sweep region -- the
    // standing column is NOT physical here -- and (ii) the physical one
    // remains at the Mb scale.
    const double omega = 894.2 / PhysConst::Hartree_eV;
    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-7;
    rpa.solve_core(omega, 60, false);
    const auto &kmat = *rpa.kmatrix();
    REQUIRE(kmat.asymmetry < 0.05);

    // 1s-channel partial cross-sections (cm^2), standing vs unitarised:
    const auto Ksigma = 4.0 * M_PI * M_PI * PhysConst::alpha *
                        PhysConst::aB_cm * PhysConst::aB_cm * omega / 3.0;
    double sig_std = 0.0, sig_uni = 0.0;
    const auto &F1s = wf.core().front();
    for (const auto &oc : rpa.open_channels(F1s)) {
      const auto D = oc.K / M_PI;
      const auto Dp = rpa.D_phys(F1s, oc.kappa);
      sig_std += Ksigma * D * D;
      sig_uni += Ksigma * Dp * Dp;
    }
    fmt::print("\nNe 1s channel at omega = 894.2 eV (near-edge sweep):\n"
               "sigma_standing = {:.3e} cm^2, sigma_unitarised = {:.3e} cm^2 "
               "(Henke ~3.5e-19)\n",
               sig_std, sig_uni);
    // Standing-wave vs physical: clearly different through the sweep (the
    // deviation grows without bound AT the pole itself, but the pole
    // position is grid-sensitive; assert a robust lower bound at this
    // sampled omega)
    REQUIRE(std::abs(sig_std / sig_uni - 1.0) > 0.15);
    // Unitarised: finite, bounded at the physical (Mb) scale
    REQUIRE(sig_uni > 1.0e-19);
    REQUIRE(sig_uni < 3.0e-18);
  }
}

//==============================================================================
//! Even-parity operators: the diagonal channel (kappa_e = kappa_a) carries
//! the norm-conservation (Lagrange multiplier) term de*phi_a in its source,
//! de = <a|(t + dV)phi_a>. The key application is scattering
//! (electron-impact ionisation), where the even K=0 temporal multipole
//! t^0 = j0(qr) usually dominates -- and for rank 0 EVERY channel is
//! diagonal. Uses the temporal vector multipole Phik (q = alpha*omega_op).
//! Checks:
//! (a) below all thresholds TDHFcntm matches bound TDHF for t^0 (the bound
//!     dispatch handles even parity via the conditioning projection);
//! (b) above threshold: TDHF converges, K = pi*D holds in the (all-diagonal)
//!     open channels, and the seeded (rescattering) solves -- which
//!     orthogonalise the diagonal seed -- give a symmetric Kbar and finite
//!     physical amplitudes;
//! (c) identity limit: j0(qr) -> 1 as q -> 0, and the identity operator
//!     must give ZERO response (norm conservation). The response is linear
//!     in the effective transition operator j0 - <j0> = O(q^2), so the K
//!     amplitudes must vanish as q^2. Without the de term the diagonal
//!     source is O(1) at small q and the response is spuriously large:
//!     this is the acceptance test for the diagonal-channel treatment.
TEST_CASE("cntmRPA: even-parity operator (temporal t0)",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  Wavefunction wf({8000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");

  // Momentum transfers (au^-1): identity limit, quadratic regime, moderate.
  // Phik takes the operator frequency: q = alpha*omega_op
  const auto t0_at_q = [&wf](double q) {
    return DiracOperator::Phik(wf.grid(), 0, q / PhysConst::alpha);
  };
  const double q_tiny = 1.0e-4;
  const double q_1 = 0.02;
  const double q_2 = 0.04;
  const double q_mod = 1.0;

  // --- (a) below-threshold: TDHFcntm matches bound TDHF (rank 0, even)
  {
    const double omega = 0.5; // 2p threshold ~0.85 au: all channels closed
    const auto t0 = t0_at_q(q_mod);
    auto rpa_bound = ExternalField::TDHF(&t0, wf.vHF());
    rpa_bound.solve_core(omega, 100, false);
    auto rpa_cntm = ExternalField::TDHFcntm(&t0, wf.vHF());
    rpa_cntm.solve_core(omega, 100, false);

    fmt::print("\nt0(qr), q = {:.3f}: bound TDHF vs TDHFcntm at omega = "
               "{:.2f} au (all closed):\n",
               q_mod, omega);
    fmt::print("{:>9s} {:>13s} {:>13s}\n", "<a|dV|b>", "TDHF", "TDHFcntm");
    for (const auto &Fa : wf.core()) {
      for (const auto &Fb : wf.core()) {
        if (t0.isZero(Fa, Fb))
          continue;
        const auto dv_b = rpa_bound.dV(Fa, Fb);
        const auto dv_c = rpa_cntm.dV(Fa, Fb);
        fmt::print("{:>4s} {:>4s} {:13.6e} {:13.6e}\n", Fa.shortSymbol(),
                   Fb.shortSymbol(), dv_b, dv_c);
        REQUIRE(dv_c == Approx(dv_b).epsilon(1.0e-2).margin(1.0e-6));
      }
    }
  }

  const double omega = 1.7; // 2p open (ec ~ 0.85 au); 2s, 1s closed

  // --- (b) above threshold: K = pi*D, and the seeded (diagonal) solves
  {
    const auto t0 = t0_at_q(q_mod);
    auto rpa = ExternalField::TDHFcntm(&t0, wf.vHF());
    rpa.eps_target() = 1.0e-6;
    rpa.solve_core(omega, 60, false);
    REQUIRE(rpa.last_eps() < 1.0e-5);

    fmt::print("\nt0(qr), q = {:.3f}: K = pi*D at omega = {:.2f} au "
               "(diagonal channels):\n",
               q_mod, omega);
    fmt::print("{:>4s} {:>6s} {:>12s} {:>12s} {:>8s}\n", "a", "kappa", "K",
               "pi*D", "K/piD-1");
    int n_channels = 0;
    for (const auto &Fb : wf.core()) {
      for (const auto &ch : rpa.open_channels(Fb)) {
        if (ch.K == 0.0)
          continue;
        const auto piD = M_PI * ch.D;
        fmt::print("{:>4s} {:>6d} {:12.5e} {:12.5e} {:8.1e}\n",
                   Fb.shortSymbol(), ch.kappa, ch.K, piD, ch.K / piD - 1.0);
        REQUIRE(ch.K == Approx(piD).epsilon(5.0e-2));
        ++n_channels;
      }
    }
    REQUIRE(n_channels == 2); // 2p_1/2 and 2p_3/2, kappa_e = kappa_a

    // Seeded homogeneous solves (the K matrix, built by solve_core): the
    // diagonal seed is orthogonalised against phi_a (constraint <a|w> = 0);
    // reciprocity (Kbar symmetry) is not imposed anywhere, so it is a real
    // consistency check of that path
    const auto &kmat = *rpa.kmatrix();
    REQUIRE(kmat.channels.size() == 2);
    REQUIRE(kmat.asymmetry < 0.05);
    for (const auto &ch : kmat.channels) {
      REQUIRE(std::isfinite(rpa.D_phys(wf.core()[ch.i_core], ch.kappa)));
    }
    fmt::print("Kbar asymmetry (diagonal-seeded solves): {:.1e}\n",
               kmat.asymmetry);
  }

  // --- (c) identity limit and q^2 scaling of the response
  {
    // K amplitudes per (2p_1/2, 2p_3/2) channel at each q
    const auto solve_Ks = [&](double q) {
      const auto t0 = t0_at_q(q);
      auto rpa = ExternalField::TDHFcntm(&t0, wf.vHF());
      rpa.eps_target() = 1.0e-6;
      rpa.set_unitarise(false); // standing K amplitudes only
      rpa.solve_core(omega, 60, false);
      std::vector<double> Ks;
      for (const auto &Fb : wf.core()) {
        for (const auto &ch : rpa.open_channels(Fb)) {
          Ks.push_back(ch.K);
        }
      }
      REQUIRE(Ks.size() == 2);
      return Ks;
    };

    const auto K_tiny = solve_Ks(q_tiny);
    const auto K_1 = solve_Ks(q_1);
    const auto K_2 = solve_Ks(q_2);
    const auto K_mod = solve_Ks(q_mod);

    fmt::print("\nt0(qr) identity limit and q^2 scaling, omega = {:.2f} au:\n",
               omega);
    fmt::print("{:>10s} {:>12s} {:>12s}\n", "q", "K(2p-)", "K(2p+)");
    for (const auto &[q, Ks] : {std::pair{q_tiny, &K_tiny},
                                {q_1, &K_1},
                                {q_2, &K_2},
                                {q_mod, &K_mod}}) {
      fmt::print("{:10.2e} {:12.5e} {:12.5e}\n", q, (*Ks)[0], (*Ks)[1]);
    }

    for (std::size_t i = 0; i < 2; ++i) {
      // Identity limit: j0 -> 1, zero response (norm conservation). With
      // the de term missing, |K(q_tiny)| ~ |K(q_mod)| instead:
      REQUIRE(std::abs(K_tiny[i]) < 1.0e-4 * std::abs(K_mod[i]));
      // Quadratic regime: K linear in (j0 - <j0>) ~ q^2:
      const auto expect = (q_1 / q_2) * (q_1 / q_2);
      REQUIRE(K_1[i] / K_2[i] == Approx(expect).epsilon(0.05));
    }
  }
}

//==============================================================================
//! Gauge invariance: exact RRPA amplitudes are identical in length (E1) and
//! velocity (E1v) form -- the decisive end-to-end test of the continuum RRPA
//! machinery (Johnson & Lin, PRA 19, 964 (1979)). Kbar is
//! operator-independent, so the test reduces to the driven amplitudes:
//! D_KS_L = D_KS_V per channel, and the unitarised |A_L| = |A_V|. The bare
//! (HF) amplitudes agree only approximately (nonlocal-exchange gauge
//! ambiguity); the RRPA resummation must close that gap. Requires the
//! consistent occupied-component convention in every channel solve (keep
//! them all; project only the diagonal) -- the old blanket core projection
//! in the open channels broke this at the several-% level.
//! E1v is frequency-dependent: t_+ at +omega, t_- at -omega.
TEST_CASE("cntmRPA: gauge invariance (E1 vs E1v)",
          "[ExternalField][TDHF][cntmrpa][integration]") {

  Wavefunction wf({5000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  const auto E1 = DiracOperator::E1(wf.grid());

  // Smooth region (2s+2p open), and high omega (all shells open):
  for (const auto omega : {4.0, 38.0}) {

    auto E1v = DiracOperator::E1v(wf.alpha(), omega);
    auto E1v_minus = DiracOperator::E1v(wf.alpha(), -omega);

    // Independent instances: each builds its own K matrix
    auto rpa_L = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa_L.eps_target() = 1.0e-7;
    rpa_L.solve_core(omega, 60, false);
    const auto &kmat_L = *rpa_L.kmatrix();

    auto rpa_V = ExternalField::TDHFcntm(&E1v, wf.vHF(), &E1v_minus);
    rpa_V.eps_target() = 1.0e-7;
    rpa_V.solve_core(omega, 60, false);
    const auto &kmat_V = *rpa_V.kmatrix();

    fmt::print("\nGauge test at omega = {:.1f} au ({} open channels):\n"
               "{:>7s} {:>11s} {:>11s} {:>9s}\n",
               omega, kmat_L.channels.size(), "kappa", "|A_L|/pi", "|A_V|/pi",
               "V/L - 1");
    REQUIRE(!kmat_L.channels.empty());
    REQUIRE(kmat_L.channels.size() == kmat_V.channels.size());

    // Kbar is operator-independent: identical matrices (machinery check)
    double dK = 0.0, kmax = 0.0;
    for (std::size_t i = 0; i < kmat_L.channels.size(); ++i) {
      for (std::size_t j = 0; j < kmat_L.channels.size(); ++j) {
        dK = std::max(dK, std::abs(kmat_L.Kbar(i, j) - kmat_V.Kbar(i, j)));
        kmax = std::max(kmax, std::abs(kmat_L.Kbar(i, j)));
      }
    }
    REQUIRE(dK < 1.0e-10 * kmax);

    // Per-channel gauge invariance of the physical amplitudes. Tolerance is
    // the numerical (extraction/TDHF) floor, not physics: observed ~1e-4.
    double sum_L = 0.0, sum_V = 0.0;
    for (const auto &ch : kmat_L.channels) {
      const auto &Fa = wf.core()[ch.i_core];
      const auto A_L = rpa_L.D_phys(Fa, ch.kappa);
      const auto A_V = rpa_V.D_phys(Fa, ch.kappa);
      fmt::print("{:>7d} {:11.6f} {:11.6f} {:9.1e}\n", ch.kappa, A_L, A_V,
                 A_V / A_L - 1.0);
      REQUIRE(A_V == Approx(A_L).epsilon(1.0e-2));
      sum_L += A_L * A_L;
      sum_V += A_V * A_V;
    }
    // Total (cross-section level): tighter
    REQUIRE(sum_V == Approx(sum_L).epsilon(2.0e-3));
  }
}

//==============================================================================
//! The batched dV builder (dV_rhs_all: shared yk screening functions) must
//! reproduce the per-task TDHF::dV_rhs exactly (same angular factors, same
//! radial integrals; only the evaluation order differs). Checked on a
//! mid-convergence state (nonzero X and Y, open and closed channels), for
//! both task types (conj), for an odd rank-1 operator (E1) and an even
//! rank-0 one (temporal t0: diagonal channels).
TEST_CASE("cntmRPA: batched dV_rhs (dV_rhs_all vs dV_rhs)",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  Wavefunction wf({2000, 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");

  const auto E1 = DiracOperator::E1(wf.grid());
  const auto t0 = DiracOperator::Phik(wf.grid(), 0, 1.0 / PhysConst::alpha);

  const std::vector<const DiracOperator::TensorOperator *> ops{&E1, &t0};
  for (const auto *h : ops) {
    auto rpa = ExternalField::TDHFcntm(h, wf.vHF());
    rpa.set_unitarise(false);
    // A few iterations: 2p/2s open at omega = 2.0, 1s closed; X and Y both
    // nonzero and far from converged (the comparison must hold at ANY state)
    rpa.solve_core(2.0, 4, false);

    const auto dv = rpa.dV_rhs_all(true);
    double worst = 0.0;
    int n_checked = 0;
    for (std::size_t ib = 0; ib < wf.core().size(); ++ib) {
      const auto &Fa = wf.core()[ib];
      const auto &Xs = rpa.get_dPsis(Fa, ExternalField::dPsiType::X);
      for (std::size_t be = 0; be < Xs.size(); ++be) {
        const auto kappa_n = Xs[be].kappa();
        for (const bool conj : {false, true}) {
          const auto direct = rpa.dV_rhs(kappa_n, Fa, conj);
          const auto &batch = conj ? dv.Y[ib][be] : dv.X[ib][be];
          const auto diff2 = (direct - batch).norm2();
          const auto ref2 = direct.norm2();
          const auto eps =
            ref2 == 0.0 ? std::sqrt(diff2) : std::sqrt(diff2 / ref2);
          worst = std::max(worst, eps);
          ++n_checked;
        }
      }
    }
    fmt::print("dV_rhs_all vs dV_rhs, {}: {} channels, worst rel diff "
               "{:.1e}\n",
               h->name(), n_checked, worst);
    REQUIRE(n_checked > 0);
    REQUIRE(worst < 1.0e-12);
  }
}

//==============================================================================
//! Marginal-band channel solve: at high omega the radial grid resolves the
//! continuum oscillations only marginally at large r (10-40 points per
//! wavelength). The pair is fine-grid solved through that band
//! (solveContinuum), so the driven solve must be too: it is continued by
//! variation of parameters on the pair (solveContinuumForward). Acceptance:
//! the K = pi*D identity, and grid-independence of the channel amplitudes
//! against a band-free reference grid. With a main-grid-only driven solve
//! both fail at O(1): the phi/pair grid errors no longer cancel in the
//! (c, K) extraction once the pair is band-accurate.
TEST_CASE("cntmRPA: marginal-band channel solve (high omega)",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  // en_+ ~ 240 au: lambda ~ 0.29 au. Coarse grid: ~25 points/wavelength at
  // r_max (band engaged); dense grid: >= 40 everywhere (band-free reference)
  const double omega = 240.0;
  const std::vector<std::size_t> npts{5000, 9000};

  // per grid: {kappa_index label -> K}, and the KpiD deviation
  std::vector<std::vector<double>> Ks(2);
  std::vector<std::vector<std::string>> labs(2);
  std::vector<double> kpid(2);

  for (std::size_t ig = 0; ig < npts.size(); ++ig) {
    Wavefunction wf({npts[ig], 1.0e-6, 40.0, 1.0, "loglinear", -1.0},
                    {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Ne]");
    const auto E1 = DiracOperator::E1(wf.grid());
    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-8;
    rpa.set_unitarise(false); // standing K amplitudes only
    rpa.solve_core(omega, 40, false);
    kpid[ig] = rpa.KpiD_dev();
    for (const auto &Fa : wf.core()) {
      for (const auto &ch : rpa.open_channels(Fa)) {
        Ks[ig].push_back(ch.K);
        labs[ig].push_back(Fa.shortSymbol() + "," +
                           AtomData::kappa_symbol(ch.kappa));
      }
    }
  }

  fmt::print("\nMarginal-band solve, omega = {:.0f} au: KpiD = {:.1e} "
             "(coarse) / {:.1e} (dense)\n",
             omega, kpid[0], kpid[1]);
  REQUIRE(kpid[0] < 0.05);
  REQUIRE(kpid[1] < 0.005);

  // Channel amplitudes grid-independent (compare against the band-free
  // reference; weak channels compared on the dominant-amplitude scale)
  REQUIRE(Ks[0].size() == Ks[1].size());
  double K_max = 0.0;
  for (const auto K : Ks[1]) {
    K_max = std::max(K_max, std::abs(K));
  }
  fmt::print("{:>10s} {:>13s} {:>13s} {:>10s}\n", "channel", "K (coarse)",
             "K (dense)", "rel diff");
  for (std::size_t i = 0; i < Ks[0].size(); ++i) {
    const auto dev = std::abs(Ks[0][i] - Ks[1][i]) / K_max;
    fmt::print("{:>10s} {:13.6e} {:13.6e} {:10.1e}\n", labs[0][i], Ks[0][i],
               Ks[1][i], dev);
    REQUIRE(dev < 0.03);
  }
}

//==============================================================================
//! Barely-open channels (en_+ far below one asymptotic wavelength in the
//! radial box): the asymptotic series does not converge at the box edge, so
//! the pair is built by the outward-extension construction (method B' in
//! solveContinuumIrregular: continue F_reg out on the H-like tail until the
//! series projection converges, seed F_irr there, integrate back in).
//! These channels used to be EXCLUDED (zeroed, warned); now they must
//! solve: SCF converged, K = pi*D, and amplitudes box-independent (the
//! physics lives inside the core region -- the box only ever mattered for
//! the pair construction).
TEST_CASE("cntmRPA: barely-open channels (extension pair)",
          "[ExternalField][TDHF][cntmrpa][unit]") {

  // en_+(2p-) fixed at 0.002 au on each grid (omega from the grid's own
  // 2p- energy): under one wavelength in either box (old exclusion was
  // en_+ < (2pi/r_max)^2/2 ~ 0.012 / 0.005); 2p+ then opens at ~0.0066
  const double en_plus = 2.0e-3;
  const std::vector<std::pair<std::size_t, double>> grids{{4000, 40.0},
                                                          {5400, 60.0}};

  std::vector<std::vector<double>> Ks(2);
  std::vector<std::vector<std::string>> labs(2);
  std::vector<double> kpid(2), eps(2);

  for (std::size_t ig = 0; ig < grids.size(); ++ig) {
    Wavefunction wf(
      {grids[ig].first, 1.0e-6, grids[ig].second, 1.0, "loglinear", -1.0},
      {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Ne]");
    const auto &F2pm = wf.core()[2]; // 2p-
    const auto omega = en_plus - F2pm.en();
    const auto E1 = DiracOperator::E1(wf.grid());
    auto rpa = ExternalField::TDHFcntm(&E1, wf.vHF());
    rpa.eps_target() = 1.0e-8;
    rpa.set_unitarise(false); // standing K amplitudes only
    rpa.solve_core(omega, 60, false);
    kpid[ig] = rpa.KpiD_dev();
    eps[ig] = rpa.last_eps();
    for (const auto &Fa : wf.core()) {
      for (const auto &ch : rpa.open_channels(Fa)) {
        Ks[ig].push_back(ch.K);
        labs[ig].push_back(Fa.shortSymbol() + "," +
                           AtomData::kappa_symbol(ch.kappa));
      }
    }
  }

  fmt::print("\nBarely-open channels, en_+(2p-) = {:.1e} au:\n"
             "  eps = {:.1e} / {:.1e},  KpiD = {:.1e} / {:.1e} "
             "(r_max 40 / 60)\n",
             en_plus, eps[0], eps[1], kpid[0], kpid[1]);
  REQUIRE(eps[0] < 1.0e-6);
  REQUIRE(eps[1] < 1.0e-6);
  REQUIRE(kpid[0] < 0.05);
  REQUIRE(kpid[1] < 0.05);

  // Channels must be present (not excluded), and box-independent
  REQUIRE(Ks[0].size() >= 5);
  REQUIRE(Ks[0].size() == Ks[1].size());
  double K_max = 0.0;
  for (const auto K : Ks[1]) {
    K_max = std::max(K_max, std::abs(K));
  }
  fmt::print("{:>10s} {:>13s} {:>13s} {:>10s}\n", "channel", "K (40 au)",
             "K (60 au)", "rel diff");
  for (std::size_t i = 0; i < Ks[0].size(); ++i) {
    REQUIRE(Ks[0][i] != 0.0);
    const auto dev = std::abs(Ks[0][i] - Ks[1][i]) / K_max;
    fmt::print("{:>10s} {:13.6e} {:13.6e} {:10.1e}\n", labs[0][i], Ks[0][i],
               Ks[1][i], dev);
    REQUIRE(dev < 0.03);
  }
}

//==============================================================================
//! The K matrix is operator-independent: it is built from field-free seeded
//! solves, so it depends only on (core, omega, rank, parity) -- never on the
//! external field's radial form. Checks, with the two E1 gauges (same rank
//! and parity, different radial operators):
//! (a) the K matrices built by independent length and velocity instances
//!     agree element-wise to machine-level precision (identical field-free
//!     solves);
//! (b) re-use: an instance switched from the length to the velocity
//!     operator (set_operator) keeps the length K matrix and channel caches,
//!     and its unitarised velocity amplitudes equal those of the independent
//!     velocity instance -- the pattern the modules rely on (one K matrix
//!     per (omega, rank, parity) serves every operator of the block);
//! (c) switching to an operator of different rank/parity drops the caches.
TEST_CASE("cntmRPA: K matrix operator-independence",
          "[ExternalField][TDHF][cntmrpa][integration]") {

  Wavefunction wf({4000, 1.0e-6, 30.0, 1.0, "loglinear", -1.0},
                  {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");

  const double omega = 1.5;
  const auto E1 = DiracOperator::E1(wf.grid());
  auto E1v = DiracOperator::E1v(wf.alpha(), omega);
  auto E1v_minus = DiracOperator::E1v(wf.alpha(), -omega);

  // Independent instances: each builds its own K matrix
  auto rpa_L = ExternalField::TDHFcntm(&E1, wf.vHF());
  auto rpa_V = ExternalField::TDHFcntm(&E1v, wf.vHF(), &E1v_minus);
  for (auto *rpa : {&rpa_L, &rpa_V}) {
    rpa->eps_target() = 1.0e-7;
    rpa->solve_core(omega, 60, false);
    REQUIRE(rpa->kmatrix().has_value());
  }
  const auto &kmat_L = *rpa_L.kmatrix();
  const auto &kmat_V = *rpa_V.kmatrix();

  // (a) identical K matrices
  REQUIRE(kmat_L.channels.size() == 5);
  REQUIRE(kmat_V.channels.size() == kmat_L.channels.size());
  double worst = 0.0;
  double kmax = 0.0;
  for (auto i = 0ul; i < kmat_L.channels.size(); i++) {
    REQUIRE(kmat_V.channels[i].i_core == kmat_L.channels[i].i_core);
    REQUIRE(kmat_V.channels[i].kappa == kmat_L.channels[i].kappa);
    for (auto j = 0ul; j < kmat_L.channels.size(); j++) {
      kmax = std::max(kmax, std::abs(kmat_L.Kbar(i, j)));
      worst = std::max(worst, std::abs(kmat_L.Kbar(i, j) - kmat_V.Kbar(i, j)));
    }
  }
  fmt::print("\nKbar(L) vs Kbar(V) at omega = {:.2f} au: worst abs diff "
             "{:.1e} (scale {:.2f})\n",
             omega, worst, kmax);
  REQUIRE(worst < 1.0e-10 * kmax);

  // (b) re-use: the length instance switched to the velocity operator
  // keeps its (length-built) K matrix and channel caches; only the driven
  // solve is redone
  auto rpa_VL = rpa_L;
  rpa_VL.set_operator(&E1v, &E1v_minus);
  REQUIRE(rpa_VL.kmatrix().has_value());
  REQUIRE(rpa_VL.channel_list().size() == kmat_L.channels.size());
  rpa_VL.solve_core(omega, 60, false);
  const auto &kmat_VL = *rpa_VL.kmatrix();
  REQUIRE(kmat_VL.channels.size() == kmat_L.channels.size());
  for (auto i = 0ul; i < kmat_L.channels.size(); i++) {
    for (auto j = 0ul; j < kmat_L.channels.size(); j++) {
      // kept, not re-built: bitwise the length one
      REQUIRE(kmat_VL.Kbar(i, j) == kmat_L.Kbar(i, j));
    }
  }
  fmt::print("{:>8s} {:>6s} {:>13s} {:>13s} {:>13s}\n", "shell", "kappa", "L",
             "V (own K)", "V (L's K)");
  for (const auto &ch : kmat_L.channels) {
    const auto &Fa = wf.core()[ch.i_core];
    const auto L_own = rpa_L.D_phys(Fa, ch.kappa);
    const auto V_own = rpa_V.D_phys(Fa, ch.kappa);
    const auto V_reused = rpa_VL.D_phys(Fa, ch.kappa);
    fmt::print("{:>8s} {:>6d} {:13.6e} {:13.6e} {:13.6e}\n", Fa.shortSymbol(),
               ch.kappa, L_own, V_own, V_reused);
    // Same driven solve (same operator, same channel pairs), K matrices
    // equal to 1e-10: machine-level agreement
    REQUIRE(V_reused == Approx(V_own).epsilon(1.0e-8));
  }

  // (c) a different (rank, parity) has different channels: caches dropped
  const auto t0 = DiracOperator::Phik(wf.grid(), 0, 1.0 / PhysConst::alpha);
  auto rpa_t0 = rpa_L;
  rpa_t0.set_operator(&t0);
  REQUIRE(rpa_t0.rank() == 0);
  REQUIRE(rpa_t0.parity() == 1);
  REQUIRE(!rpa_t0.kmatrix().has_value());
  REQUIRE(rpa_t0.channel_list().empty());
}
