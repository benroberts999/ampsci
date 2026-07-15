#include "DiracOperator/include.hpp"
#include "ExternalField/TDHFcntm.hpp"
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

    // New continuum TDHF, at the various levels:
    auto rpa_g0 = ExternalField::TDHFcntm(&E1, wf.vHF());  // fwd solve, no dV
    rpa_g0.solve_core(omega, 0);                           // 0 iterations
    auto rpa_1it = ExternalField::TDHFcntm(&E1, wf.vHF()); // one dV iteration
    rpa_1it.solve_core(omega, 1);
    auto rpa_closed = ExternalField::TDHFcntm(&E1, wf.vHF()); // open chnls off
    rpa_closed.set_suppress_open(true);
    rpa_closed.solve_core(omega, 60);
    auto rpa_full = ExternalField::TDHFcntm(&E1, wf.vHF()); // full dV
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
        const auto D_g0 = D_vn1 + rpa_g0.dV_cntm(Fe, Fb);
        const auto D_1it = D_vn1 + rpa_1it.dV_cntm(Fe, Fb);
        const auto D_closed = D_vn1 + rpa_closed.dV_cntm(Fe, Fb);
        const auto D_full = D_vn1 + rpa_full.dV_cntm(Fe, Fb);

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
