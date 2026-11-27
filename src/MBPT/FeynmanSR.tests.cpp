#include "MBPT/FeynmanSR.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "Modules/SRBreakdown.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <optional>
#include <utility>

//==============================================================================
TEST_CASE("FeynSR: perturbed Greens function", "[FeynSR][MBPT][unit]") {
  // For Hartree-Fock eigenstates, G(e)|v> = |v>/(e - e_v), so
  //   <w| G(e1) h G(e2) |v> = <w||h||v> / [(e1 - e_w)*(e2 - e_v)].
  // Tests the Green's functions, the radial part of the operator, and its
  // angular factor

  // Li: Hartree-Fock core, with 2s and 2p valence states
  Wavefunction wf({1600, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                  {"Li", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-14, false);
  wf.solve_valence("2sp", false);

  // Sub-grid: every 6th point over r = 1e-3 to 30. With g components.
  // Q*Pi*Q is not needed for the Green's functions, so is not formed
  const std::size_t stride = 6;
  const auto i0 = wf.grid().getIndex(1.0e-3);
  const auto num_subgrid_points = (wf.grid().getIndex(30.0) - i0) / stride + 1;
  MBPT::FeynmanOptions options{};
  options.max_l_internal = 3;
  const MBPT::FeynmanSR Fy(wf.vHF(), i0, stride, num_subgrid_points, options, 1,
                           true, false, "", false);

  const auto h = DiracOperator::E1(wf.grid());

  // Complex energies, as on the integration contour
  const std::complex<double> e1{-0.3, 0.2};
  const std::complex<double> e2{-0.5, 0.2};

  for (const auto &w : wf.valence()) {
    for (const auto &v : wf.valence()) {
      if (h.isZero(w.kappa(), v.kappa()))
        continue;

      // dG = G_w(e1) * h * G_v(e2)
      const auto dG =
        Fy.perturbed_greens(h, w.kappa(), v.kappa(), Fy.green(w.kappa(), e1),
                            Fy.green(v.kappa(), e2));

      // <w|dG|v>: the matrix must carry the measure for the inner integral
      const std::complex<double> found{w * (dG.real().drj() * v),
                                       w * (dG.imag().drj() * v)};

      const auto expected = h.reducedME(w, v) / ((e1 - w.en()) * (e2 - v.en()));

      const auto eps = std::abs(found / expected - 1.0);
      REQUIRE(eps < 1.0e-3);
    }
  }
}

//==============================================================================
TEST_CASE("FeynSR: core polarisation in perturbed Greens function",
          "[FeynSR][MBPT][unit]") {
  // With core polarisation, the external field vertex is h + dV. dV is
  // non-local, so is written as a matrix on the sub-grid. As above:
  //   <w| G(e1) dV G(e2) |v> = <w||dV||v> / [(e1 - e_w)*(e2 - e_v)],
  // where the dV part is the difference of dG with and without dV

  // Li: Hartree-Fock core, with 2s and 2p valence states
  Wavefunction wf({1600, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                  {"Li", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-14, false);
  wf.solve_valence("2sp", false);

  // Sub-grid: every 6th point over r = 1e-3 to 30. With g components.
  // Q*Pi*Q is not needed for the Green's functions, so is not formed
  const std::size_t stride = 6;
  const auto i0 = wf.grid().getIndex(1.0e-3);
  const auto num_subgrid_points = (wf.grid().getIndex(30.0) - i0) / stride + 1;
  MBPT::FeynmanOptions options{};
  options.max_l_internal = 3;
  MBPT::FeynmanSR Fy(wf.vHF(), i0, stride, num_subgrid_points, options, 1, true,
                     false, "", false);

  const auto h = DiracOperator::E1(wf.grid());

  // Core polarisation (TDHF), solved at zero frequency
  const auto dV = ExternalField::make_rpa("TDHF", &h, wf.vHF(), false);
  dV->solve_core(0.0, 100, false);

  // dV as a matrix, for each pair of partial waves: indexed by kappa index
  const auto dV_kappa = Fy.dV_matrices(h, dV.get());

  const std::complex<double> e1{-0.3, 0.2};
  const std::complex<double> e2{-0.5, 0.2};

  for (const auto &w : wf.valence()) {
    for (const auto &v : wf.valence()) {
      if (h.isZero(w.kappa(), v.kappa()))
        continue;

      const auto G_w = Fy.green(w.kappa(), e1);
      const auto G_v = Fy.green(v.kappa(), e2);
      const auto &dV_wv = dV_kappa(w.k_index(), v.k_index());

      // G_w * (h + dV) * G_v, and G_w * h * G_v
      const auto dG_h_dV =
        Fy.perturbed_greens(h, w.kappa(), v.kappa(), G_w, G_v, dV_wv);
      const auto dG_h = Fy.perturbed_greens(h, w.kappa(), v.kappa(), G_w, G_v);
      const auto dG_dV = dG_h_dV - dG_h;

      const std::complex<double> found{w * (dG_dV.real().drj() * v),
                                       w * (dG_dV.imag().drj() * v)};

      const auto expected = dV->dV(w, v) / ((e1 - w.en()) * (e2 - v.en()));

      const auto eps = std::abs(found / expected - 1.0);
      REQUIRE(eps < 1.0e-3);
    }
  }
}

//==============================================================================
TEST_CASE("FeynSR: diagrams vs Goldstone", "[FeynSR][MBPT][integration]") {
  // Structure radiation diagrams (a), (b), (c), (a_x) and (c_x), with core
  // polarisation of the external field, against the same diagrams from a
  // sum over a B-spline basis.
  // Small sub-grid and basis: tests the method, not high accuracy

  // Li: Hartree-Fock core, with 2s and 2p valence states, and basis
  Wavefunction wf({1600, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                  {"Li", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-14, false);
  wf.solve_valence("2sp", false);
  wf.formBasis(SplineBasis::Parameters("30spdf", 40, 7, {1.0e-4}, 40.0, "",
                                       SplineBasis::SplineType::Derevianko,
                                       false, false));

  const auto h = DiracOperator::E1(wf.grid());

  // Core polarisation (TDHF), solved at zero frequency
  const auto dV = ExternalField::make_rpa("TDHF", &h, wf.vHF(), false);
  dV->solve_core(0.0, 100, false);

  // Feynman: sub-grid of every 6th point over r = 1e-3 to 30; l <= 3 for
  // the internal lines; no g components. Nothing is written to disk
  const std::size_t stride = 6;
  const auto i0 = wf.grid().getIndex(1.0e-3);
  const auto num_subgrid_points = (wf.grid().getIndex(30.0) - i0) / stride + 1;
  MBPT::FeynmanOptions options{};
  options.max_l_internal = 3;
  options.omre = -std::abs(MBPT::best_omre(wf.core(), wf.valence()));
  MBPT::FeynmanSR Fy(wf.vHF(), i0, stride, num_subgrid_points, options, 1,
                     false, false, "");

  // Sum over basis states: Coulomb integrals, and matrix elements of h + dV
  auto [core, excited] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel(), 1);
  auto basis = core;
  basis.insert(basis.end(), excited.begin(), excited.end());
  Coulomb::YkTable yk;
  yk.extend_angular(DiracSpinor::max_tj(basis));
  yk.calculate(basis);
  const auto h_table = Amplitudes::me_table(basis, &h, dV.get(), nullptr, 0.0);

  const auto &F2s = *wf.getState("2s");
  const auto &F2p = *wf.getState("2p-");

  fmt::print("{:3s} {:3s} {:4s} {:>11s} {:>11s} {:>8s}\n", "w", "v", "",
             "Basis", "Feynman", "eps");

  // Both orderings: (a_x) and (c_x) swap
  for (const auto &[Fw, Fv] : {std::pair{&F2s, &F2p}, std::pair{&F2p, &F2s}}) {
    const auto omega = Fw->en() - Fv->en();

    // Basis: the external states are the basis states matching w and v
    const auto &w = *std::find(excited.cbegin(), excited.cend(), *Fw);
    const auto &v = *std::find(excited.cbegin(), excited.cend(), *Fv);
    const auto [a, c, ax, cx] = Module::SR_a_c_ax_cx(
      w, v, omega, core, excited, yk, h_table, h.rank(), {}, {});
    const auto [b, bx] = Module::SR_b_bx(w, v, omega, core, excited, yk,
                                         h_table, h.rank(), {}, {});

    // Feynman: each returns a matrix; the reduced matrix element is
    // Fw * (dSigma * Fv)
    const auto dSigma_a =
      Fy.dSigma_a(h, Fw->kappa(), Fv->kappa(), Fv->en(), omega, dV.get());
    const auto [dSigma_b, dSigma_c] =
      Fy.dSigma_bc(h, Fw->kappa(), Fv->kappa(), Fv->en(), omega, dV.get());
    const auto [dSigma_ax, dSigma_cx] =
      Fy.dSigma_ax_cx(h, Fw->kappa(), Fv->kappa(), Fv->en(), omega, dV.get());
    const auto a_Feynman = *Fw * (dSigma_a * *Fv);
    const auto b_Feynman = *Fw * (dSigma_b * *Fv);
    const auto c_Feynman = *Fw * (dSigma_c * *Fv);
    const auto ax_Feynman = *Fw * (dSigma_ax * *Fv);
    const auto cx_Feynman = *Fw * (dSigma_cx * *Fv);

    const auto eps_a = std::abs(a_Feynman / a - 1.0);
    const auto eps_b = std::abs(b_Feynman / b - 1.0);
    const auto eps_c = std::abs(c_Feynman / c - 1.0);
    const auto eps_ax = std::abs(ax_Feynman / ax - 1.0);
    const auto eps_cx = std::abs(cx_Feynman / cx - 1.0);

    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n",
               Fw->shortSymbol(), Fv->shortSymbol(), "a", a, a_Feynman, eps_a);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "b",
               b, b_Feynman, eps_b);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "c",
               c, c_Feynman, eps_c);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "ax",
               ax, ax_Feynman, eps_ax);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "cx",
               cx, cx_Feynman, eps_cx);

    REQUIRE(eps_a < 0.01);
    REQUIRE(eps_b < 0.01);
    REQUIRE(eps_c < 0.01);
    REQUIRE(eps_ax < 0.01);
    REQUIRE(eps_cx < 0.01);
  }
}
