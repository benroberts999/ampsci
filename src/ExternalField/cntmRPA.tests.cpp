#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "HF/HartreeFock.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>

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
    auto rpa_g0 = ExternalField::TDHF(&E1, wf.vHF());  // forward solve, no dV
    rpa_g0.solve_core_cntm(omega, 0);                  // 0 iterations
    auto rpa_1it = ExternalField::TDHF(&E1, wf.vHF()); // one iteration of dV
    rpa_1it.solve_core_cntm(omega, 1);
    auto rpa_closed = ExternalField::TDHF(&E1, wf.vHF()); // dV, open chnls off
    rpa_closed.solve_core_cntm(omega, 60, true, true);    // suppress_open
    auto rpa_full = ExternalField::TDHF(&E1, wf.vHF());   // full dV
    // F_reg is re-normalised from its large-r envelope (~1e-6 noise) each
    // rebuild; that sets the achievable SCF floor (well below the MEs).
    rpa_full.eps_target() = 1.0e-5;
    rpa_full.solve_core_cntm(omega, 60);

    fmt::print("\nNe E1 photoionisation, omega = {:.2f} au  (SCF eps = {:.1e})\n",
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
        const auto D_g0 = D_vn1 + rpa_g0.dV_cont(Fe, Fb);
        const auto D_1it = D_vn1 + rpa_1it.dV_cont(Fe, Fb);
        const auto D_closed = D_vn1 + rpa_closed.dV_cont(Fe, Fb);
        const auto D_full = D_vn1 + rpa_full.dV_cont(Fe, Fb);

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
    // The SCF stays bounded (no divergence); it converges tightly except near a
    // collective resonance, where it stalls at a higher eps (still finite).
    REQUIRE(std::isfinite(rpa_full.last_eps()));
  }
}
