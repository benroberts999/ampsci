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
  // Structure radiation diagrams (a), (b), (c), (a_x), (b_x) and (c_x), with
  // core polarisation of the external field, against the same diagrams from
  // a sum over a B-spline basis.
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

  // Loop Green's functions: the same for every matrix element and diagram.
  // Optional: otherwise formed by the first diagram that needs them
  Fy.calculate_gex();

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
    const auto dSigma_bx =
      Fy.dSigma_bx(h, Fw->kappa(), Fv->kappa(), Fv->en(), omega, dV.get());
    const auto bx_Feynman = *Fw * (dSigma_bx * *Fv);

    // The total: all six diagrams from the one function
    const auto dSigma_total =
      Fy.dSigma_total(h, Fw->kappa(), Fv->kappa(), Fv->en(), omega, dV.get());
    const auto total_Feynman = *Fw * (dSigma_total * *Fv);
    const auto sum_Feynman =
      a_Feynman + b_Feynman + c_Feynman + ax_Feynman + bx_Feynman + cx_Feynman;
    const auto total = a + b + c + ax + bx + cx;

    const auto eps_a = std::abs(a_Feynman / a - 1.0);
    const auto eps_b = std::abs(b_Feynman / b - 1.0);
    const auto eps_c = std::abs(c_Feynman / c - 1.0);
    const auto eps_ax = std::abs(ax_Feynman / ax - 1.0);
    const auto eps_bx = std::abs(bx_Feynman / bx - 1.0);
    const auto eps_cx = std::abs(cx_Feynman / cx - 1.0);
    const auto eps_total = std::abs(total_Feynman / total - 1.0);
    const auto eps_sum = std::abs(total_Feynman / sum_Feynman - 1.0);

    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n",
               Fw->shortSymbol(), Fv->shortSymbol(), "a", a, a_Feynman, eps_a);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "b",
               b, b_Feynman, eps_b);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "c",
               c, c_Feynman, eps_c);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "ax",
               ax, ax_Feynman, eps_ax);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "bx",
               bx, bx_Feynman, eps_bx);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "cx",
               cx, cx_Feynman, eps_cx);
    fmt::print("{:3s} {:3s} {:4s} {:+11.4e} {:+11.4e} {:8.1e}\n", "", "", "sum",
               total, total_Feynman, eps_total);
    fmt::print("dSigma_total vs sum of the six diagrams: {:.1e}\n", eps_sum);

    REQUIRE(eps_a < 0.01);
    REQUIRE(eps_b < 0.01);
    REQUIRE(eps_c < 0.01);
    REQUIRE(eps_ax < 0.01);
    // (b_x) converges more slowly with the sub-grid
    REQUIRE(eps_bx < 0.05);
    REQUIRE(eps_cx < 0.01);
    // The diagrams partly cancel in the sum
    REQUIRE(eps_total < 0.05);
    REQUIRE(eps_sum < 1.0e-9);
  }
}

//==============================================================================
TEST_CASE("FeynSR: stored loop Greens functions", "[FeynSR][MBPT]") {
  // The excited Green's functions of the core loops are taken from the
  // stored set when there is one, and are solved as needed otherwise. The
  // diagrams must be the same either way.
  // Static field, and field at the transition frequency

  // Li: Hartree-Fock core, with 2s and 2p valence states
  Wavefunction wf({1600, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                  {"Li", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-14, false);
  wf.solve_valence("2sp", false);

  const auto h = DiracOperator::E1(wf.grid());

  // Small sub-grid; l <= 2 for the internal lines; hole-particle
  // interaction included (it is part of the stored Green's functions).
  // Nothing is written to disk
  const std::size_t stride = 8;
  const auto i0 = wf.grid().getIndex(1.0e-3);
  const auto num_subgrid_points = (wf.grid().getIndex(30.0) - i0) / stride + 1;
  MBPT::FeynmanOptions options{};
  options.max_l_internal = 2;
  options.omre = -std::abs(MBPT::best_omre(wf.core(), wf.valence()));
  options.hole_particle = MBPT::HoleParticle::include;

  options.store_gex = true;
  MBPT::FeynmanSR Fy_stored(wf.vHF(), i0, stride, num_subgrid_points, options,
                            1, false, false, "");
  options.store_gex = false;
  MBPT::FeynmanSR Fy_solved(wf.vHF(), i0, stride, num_subgrid_points, options,
                            1, false, false, "");

  const auto &Fw = *wf.getState("2p-");
  const auto &Fv = *wf.getState("2s");

  // Forms the stored set (no effect for Fy_solved: store_gex is false)
  Fy_stored.calculate_gex();
  Fy_solved.calculate_gex();

  fmt::print("{:4s} {:>9s} {:>12s} {:>12s} {:>8s}\n", "", "omega", "solved",
             "stored", "eps");

  for (const auto omega : {0.0, Fw.en() - Fv.en()}) {

    const auto [b_solved, c_solved] =
      Fy_solved.dSigma_bc(h, Fw.kappa(), Fv.kappa(), Fv.en(), omega);
    const auto [ax_solved, cx_solved] =
      Fy_solved.dSigma_ax_cx(h, Fw.kappa(), Fv.kappa(), Fv.en(), omega);
    const auto bx_solved =
      Fy_solved.dSigma_bx(h, Fw.kappa(), Fv.kappa(), Fv.en(), omega);

    const auto [b_stored, c_stored] =
      Fy_stored.dSigma_bc(h, Fw.kappa(), Fv.kappa(), Fv.en(), omega);
    const auto [ax_stored, cx_stored] =
      Fy_stored.dSigma_ax_cx(h, Fw.kappa(), Fv.kappa(), Fv.en(), omega);
    const auto bx_stored =
      Fy_stored.dSigma_bx(h, Fw.kappa(), Fv.kappa(), Fv.en(), omega);

    const auto b_1 = Fw * (b_solved * Fv);
    const auto b_2 = Fw * (b_stored * Fv);
    const auto c_1 = Fw * (c_solved * Fv);
    const auto c_2 = Fw * (c_stored * Fv);
    const auto ax_1 = Fw * (ax_solved * Fv);
    const auto ax_2 = Fw * (ax_stored * Fv);
    const auto cx_1 = Fw * (cx_solved * Fv);
    const auto cx_2 = Fw * (cx_stored * Fv);
    const auto bx_1 = Fw * (bx_solved * Fv);
    const auto bx_2 = Fw * (bx_stored * Fv);

    const auto eps_b = std::abs(b_2 / b_1 - 1.0);
    const auto eps_c = std::abs(c_2 / c_1 - 1.0);
    const auto eps_ax = std::abs(ax_2 / ax_1 - 1.0);
    const auto eps_cx = std::abs(cx_2 / cx_1 - 1.0);
    const auto eps_bx = std::abs(bx_2 / bx_1 - 1.0);

    fmt::print("{:4s} {:9.6f} {:+12.5e} {:+12.5e} {:8.1e}\n", "b", omega, b_1,
               b_2, eps_b);
    fmt::print("{:4s} {:9.6f} {:+12.5e} {:+12.5e} {:8.1e}\n", "c", omega, c_1,
               c_2, eps_c);
    fmt::print("{:4s} {:9.6f} {:+12.5e} {:+12.5e} {:8.1e}\n", "ax", omega, ax_1,
               ax_2, eps_ax);
    fmt::print("{:4s} {:9.6f} {:+12.5e} {:+12.5e} {:8.1e}\n", "cx", omega, cx_1,
               cx_2, eps_cx);
    fmt::print("{:4s} {:9.6f} {:+12.5e} {:+12.5e} {:8.1e}\n", "bx", omega, bx_1,
               bx_2, eps_bx);

    REQUIRE(eps_b < 1.0e-9);
    REQUIRE(eps_c < 1.0e-9);
    REQUIRE(eps_ax < 1.0e-9);
    REQUIRE(eps_cx < 1.0e-9);
    REQUIRE(eps_bx < 1.0e-9);
  }
}
