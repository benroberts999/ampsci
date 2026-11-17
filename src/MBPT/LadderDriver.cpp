#include "MBPT/LadderDriver.hpp"
#include "Angular/include.hpp"
#include "CI/CI_Integrals.hpp"
#include "Coulomb/include.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/Ladder.hpp"
#include "MBPT/LadderExternalLine.hpp"
#include "MBPT/LadderIteration.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace MBPT {

//==============================================================================
// Driver for ladder diagram calculation: Ladder{} input block
void ladder(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"min_n_core", "lowest core n to include [1]"},
     {"basis", "Basis string specifying which states to include; must be a "
               "subset of full "
               "basis), e.g. '30spdf'. Default: include entire excited basis."},
     {"max_k", "maximum k to include in Qk. Put -1 to include all. [8]"},
     {"include_L4", "Inlcude 4th Ladder diagram [true]"},
     {"include_P", "Include core polarisation in the ladder [false]:\n"
                   " - with its exchange and hole-particle parts\n"
                   " - no Sigma_L: prints valence energy corrections,\n"
                   "   without the part in the Feynman direct diagram\n"
                   "   (Correlations with screening and hole_particle)\n"
                   " - use min_n_core = Correlations n_min_core\n"
                   " - needs max_it > 0"},
     {"include_R", "Include the core relaxation [false]:\n"
                   " - core single excitations in the ladder\n"
                   " - adds the interaction with the change of the\n"
                   "   core density, <v|dV^R|v>, to the valence\n"
                   "   energy corrections\n"
                   " - needs include_P"},
     {"sigma_method",
      "Method for the ladder correlation potential Sigma_L. Options:\n"
      "   - basis  : projection Sigma_L = sum_x L_x|x><x| onto the excited "
      "states x of kappa_v, with the ladder iterated for each x at the "
      "valence energy (uses max_it, damp and eps_target)\n"
      "   - ratio  : no projection; rescale each Sigma(2) term by L/Q "
      "(fast)\n"
      "   - ratio2 : exact acting on the valence state; ratio only for the "
      "change of the valence orbital in the Brueckner equation (fast)\n"
      "  [basis]"},
     {"mirror_terms",
      "Add the mirror terms to Sigma_L: the mirror images of the ladder "
      "diagrams that the ladder equation does not generate (lowest order "
      "only). Their off-diagonal part is projected onto every basis state "
      "of kappa_v, core states included: complete up to the basis (the "
      "positive-energy spline states of kappa_v in the cavity). The Coulomb "
      "integrals of those states are calculated each run, not stored in the "
      "Qk file [false]"},
     {"each_valence",
      "If true, iterate L and calculate Sigma_L for each valence state "
      "separately. If false, only for the lowest valence state of each kappa "
      "(Correlations applies it to all states of that kappa). [false]"},
     {"Qk_file", "Filename for storing Qk Coulomb integrals. By default, is "
                 "<Identity>.qk. If 'false' will not read or write."},
     {"Lk_file", "Filename for storing Lk ladder integrals. By default, is "
                 "<Identity>_l4.lk (<Identity>.lk if include_L4 is false; _P "
                 "added before .lk with include_P, _PR with include_R). Do "
                 "not use one file for runs with and without include_P or "
                 "include_R. If 'false' will not read or write."},
     {"sl_file", "Filename for storing the ladder correlation potential, "
                 "Sigma_L. By default, is <Identity>.sl.abf. If 'false' will "
                 "not write. Read in via Correlations{ladder=...;}"},
     {"from_scratch", "If true, don't read existing Qk/Lk files (still "
                      "writes). [false]"},
     {"max_it", "Max # iterations. If zero, will simply read ladder diagrams "
                "in (any _new_ ladder diagrams will be calculated,). Must be > "
                "0 with include_P. [32]"},
     {"damp",
      "Damping factor for iterations, [0,1). 0 means no damping. [0.0]"},
     {"eps_target", "Target for convergance [1.0e-5]"},
     {"rmin", "minimum radius of Sigma_L sub-grid [1.0e-4]"},
     {"rmax", "maximum radius of Sigma_L sub-grid [30.0]"},
     {"stride", "Only calculate Sigma_L every <stride> points. Default such "
                "that there are 150 points between (1e-4, 30)"},
     {"include_G", "Inlcude lower g-part into Sigma_L [false]"},
     {"check_symmetry",
      "Print the asymmetry of Sigma_L, before it is symmetrised [false]"}});
  // If we are just requesting 'help', don't run:
  if (input.has_option("help")) {
    return;
  }

  std::cout << "\n";
  IO::print_line('*', 50);
  IO::print_line('*', 50);
  std::cout << "Ladder Diagrams\n\n";

  // Get input
  const auto min_n_core = input.get("min_n_core", 1);
  const auto basis_str = input.get<std::string>("basis", "");
  const auto max_k = input.get("max_k", 8);
  const auto include_L4 = input.get("include_L4", true);
  const auto include_P = input.get("include_P", false);
  const auto include_R = input.get("include_R", false);
  // Diagrams of the ladder equation: L1-L4 (L4 optional), and the three
  // pieces of P with include_P
  const MBPT::LadderDiagrams diagrams{
    true, true, true, include_L4, include_P, include_P, include_P};

  // Method for Sigma_L: basis projection, or ratio (rescale Sigma(2) terms by
  // L/Q; Dzuba)
  using namespace std::string_literals;
  const auto method = parseSigmaLMethod(input.get("sigma_method", "basis"s));

  const auto mirror_terms = input.get("mirror_terms", false);

  const auto each_valence = input.get("each_valence", false);

  const auto max_it = input.get("max_it", 32);
  const auto a_damp = input.get("damp", 0.0);
  const auto eps_target = input.get("eps_target", 1.0e-5);

  // Sub-grid for Sigma_L matrix (same defaults as Correlations):
  const auto sig_r0 = input.get("rmin", 1.0e-4);
  const auto sig_rmax = input.get("rmax", 30.0);
  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto stride =
      int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;
    return (stride <= 2) ? 2 : stride;
  }();
  const auto sig_stride = std::size_t(input.get("stride", default_stride));
  const auto include_G = input.get("include_G", false);

  const auto check_symmetry = input.get("check_symmetry", false);

  // Sort basis into core/excited/valence
  const auto core_excited =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel(), min_n_core);
  const auto &holes = core_excited.first;
  const auto excited =
    basis_str.empty() ?
      core_excited.second :
      DiracSpinor::subset(core_excited.second, basis_str, false);

  // nb: use _basis_ version of valence states in iterations.
  // Take only the lowest valence state of each kappa (Unless each_valence=true)
  std::vector<DiracSpinor> valence;
  const auto max_ki = DiracSpinor::max_kindex(wf.valence());
  for (std::size_t ki = 0; ki <= max_ki; ++ki) {
    const auto kappa = Angular::kindex_to_kappa(ki);
    for (const auto &Fv : wf.valence()) {
      if (Fv.kappa() != kappa)
        continue;
      const auto pFv = std::find(wf.basis().cbegin(), wf.basis().cend(), Fv);
      if (pFv == wf.basis().cend()) {
        std::cout << "Warning: Basis missing valence state: " << Fv << "\n";
      } else {
        valence.push_back(*pFv);
      }
      if (!each_valence)
        break;
    }
  }

  std::cout << "basis        = " << DiracSpinor::state_config(excited) << "\n";
  std::cout << "min_n (core) = " << min_n_core << "\n";
  std::cout << std::boolalpha;
  std::cout << "include_L4   = " << include_L4 << "\n";
  std::cout << "include_P    = " << include_P << "\n";
  std::cout << "include_R    = " << include_R << "\n";
  std::cout << "sigma_method = " << parseSigmaLMethod(method) << "\n";
  std::cout << "mirror_terms = " << mirror_terms << "\n";
  std::cout << "each_valence = " << each_valence << "\n";
  std::cout << "include_G    = " << include_G << "\n";
  std::cout << "max_k        = " << max_k << "\n";
  std::cout << "max_it       = " << max_it << "\n";
  std::cout << "damp         = " << a_damp << "\n";
  std::cout << "eps_target   = " << eps_target << "\n";

  // x_F is not stored: it is iterated in every run, from the lowest order
  if (include_P && max_it < 1) {
    std::cout << "\nFAIL: include_P requires max_it > 0\n";
    return;
  }
  // The energy with R is formed only on the include_P path
  if (include_R && !include_P) {
    std::cout << "\nFAIL: include_R requires include_P\n";
    return;
  }

  // in/out file names (default based on atom identity)
  const auto ident = wf.identity();
  const auto Qk_file = input.get("Qk_file", ident + ".qk.abf"s);
  const auto lk4 = include_L4 ? "_l4" : ""s;
  const auto lkP = include_R ? "_PR" : include_P ? "_P" : ""s;
  const auto Lk_file =
    input.get<std::string>("Lk_file", ident + lk4 + lkP + ".lk.abf"s);
  const auto sl_file = input.get<std::string>("sl_file", ident + ".sl.abf"s);
  const auto from_scratch = input.get("from_scratch", false);

  // Create the "total" basis, which has core+excited, but only those states
  // actually included (i.e., [n_min, n_max]). This is used to calculate Qk.
  // Reduces size of Qk; should also reduce lookup time.
  const auto both = qip::merge(holes, excited);
  // the "i" orbitals: core + valence
  const auto core_and_val = qip::merge(holes, valence);

  // Fill Yk table (used to fill Qk table)
  const Coulomb::YkTable yk(both);
  // Store reference to Yk's 6J table (save typing):
  const auto &sjt = yk.SixJ();

  // Form the Qk table.
  Coulomb::QkTable qk;
  std::cout << "\nFill (or read) Qk table:\n";
  if (!from_scratch) {
    qk.read(Qk_file);
  }
  const auto count_init = qk.count();
  qk.fill(both, yk, max_k);
  const auto count_after = qk.count();
  if (count_init != count_after) {
    fmt::print("Calculated {} new Qk integrals\n", count_after - count_init);
    qk.write(Qk_file);
  } else {
    std::cout << "No new Qk integrals required\n";
  }

  std::cout << "\nMBPT(2) energy shifts, using HF vs. spline legs" << std::endl;
  fmt::print("{:<5s} {:>13s} {:>13s} {:>13s}   {}\n", "state", "de(HF)",
             "de(basis)", "de(Qk)", "eps");
  const auto dec_1 = MBPT::de_core(yk, yk, wf.core(), excited);
  const auto dec_2 = MBPT::de_core(yk, yk, holes, excited);
  const auto dec_3 = MBPT::de_core(qk, qk, holes, excited);
  const auto dec_eps = std::abs(2.0 * (dec_1 - dec_2) / (dec_1 + dec_2));
  fmt::print("{:<5s} {:+13.7f} {:+13.7f} {:+13.7f}   {:.0e}\n", "core", dec_1,
             dec_2, dec_3, dec_eps);
  for (const auto &sv : valence) {
    const auto pv = wf.getState(sv.n(), sv.kappa());
    assert(pv != nullptr); //valence is subset of wf.valence()
    const auto de_1 = MBPT::de_valence(*pv, yk, yk, holes, excited);
    const auto de_2 = MBPT::de_valence(sv, yk, yk, holes, excited);
    const auto de_3 = MBPT::de_valence(sv, qk, qk, holes, excited);
    const auto de_eps = std::abs(2.0 * (de_1 - de_2) / (de_1 + de_2));
    fmt::print("{:<5s} {:+13.7f} {:+13.7f} {:+13.7f}   {:.0e}\n",
               sv.shortSymbol(), de_1, de_2, de_3, de_eps);
  }
  std::cout << "\n";

  // Fill Lk table:
  Coulomb::LkTable lk;
  Coulomb::LkTable lk_next;
  // Each run will continue from where last left off
  const bool did_read_Lk = !from_scratch ? lk.read(Lk_file) : false;
  std::cout << (did_read_Lk ? "\nRe-starting using existing ladder diagrams\n" :
                              "\nCalculating ladder diagrams from scratch\n");

  //----------------------------------------------------------------------------
  // Iterate Lk equations

  // "First iteration": fill_Lk_mnib keeps existing (read) entries and adds only
  // newly-required diagrams (or fills from scratch if not read). Subsequent
  // iterations below re-iterate the existing table.
  std::cout << "\nInitial fill of Lk table: core + valence\n" << std::flush;
  const auto n_lk_initial = lk.count();
  MBPT::fill_Lk_mnib(&lk, qk, excited, holes, core_and_val, diagrams, max_k,
                     true);
  const auto n_lk_after = lk.count();
  if (n_lk_initial != n_lk_after) {
    fmt::print("Calculated {} new Lk integrals\n", n_lk_after - n_lk_initial);
    lk.summary();
    lk.write(Lk_file);
  } else {
    std::cout << "No new Lk integrals required\n";
  }
  lk_next = lk;

  // convert to inverse cm
  const auto icm = PhysConst::Hartree_invcm;

  // initial corrections
  double de_c0 = MBPT::de_core(qk, lk, holes, excited);
  std::vector<double> de_0;
  std::cout << "\nInitial Ladder Corrections:\n";
  fmt::print("de_l({:4}): {:+12.9f} au  {:+12.5f} cm^-1\n", "core", de_c0,
             de_c0 * icm);
  for (const auto &v : valence) {
    const auto de_v = MBPT::de_valence(v, qk, lk, holes, excited);
    de_0.push_back(de_v);
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n", v.shortSymbol(),
               de_v, de_v * icm);
  }

  // With P the change per sweep is not monotone: converged only after this
  // many consecutive sweeps below eps_target
  const int required_sweeps_below_eps = include_P ? 2 : 1;

  // Core single excitations rho_ma (include_R): rows excited, columns core.
  // Not stored: solved for the x the sweeps start from (read, or lowest
  // order), then updated with x in each core sweep
  LinAlg::Matrix<double> rho_ma;
  if (include_R) {
    std::cout << "\nCore single excitations for the starting x:\n";
    const auto source = MBPT::core_singles_source(qk, excited, holes, &lk);
    rho_ma = MBPT::core_singles(source, qk, excited, holes, nullptr);
    int steps = 0;
    double eps_rho = 1.0;
    for (int it = 1; it <= max_it && eps_rho >= eps_target; ++it) {
      auto rho_next = MBPT::core_singles(source, qk, excited, holes, &rho_ma);
      eps_rho = MBPT::core_singles_change(rho_next, rho_ma);
      rho_ma = std::move(rho_next);
      steps = it;
    }
    fmt::print("rho_ma: {} steps, change {:.1e}\n", steps, eps_rho);
  }
  const auto *const rho_ma_or_null = include_R ? &rho_ma : nullptr;

  // Iterate for core states:
  int core_sweeps_below_eps = 0;
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "\ncore it:" << it << "\n";
    MBPT::update_Lk_mnib(&lk_next, qk, excited, holes, holes, diagrams, &lk,
                         a_damp, true, rho_ma_or_null);
    // for the damping, don't swap; lk and lk_prev must match before next iter
    lk = lk_next;
    lk.write(Lk_file);
    const auto de_c = MBPT::de_core(qk, lk, holes, excited);
    const auto eps_c = std::abs(2.0 * (de_c - de_c0) / (de_c + de_c0));
    de_c0 = de_c;
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n", "core",
               de_c, de_c * icm, eps_c);
    // The core single excitations from the new x, one step
    double eps_rho = 0.0;
    if (include_R) {
      auto rho_next =
        MBPT::core_singles(MBPT::core_singles_source(qk, excited, holes, &lk),
                           qk, excited, holes, &rho_ma);
      eps_rho = MBPT::core_singles_change(rho_next, rho_ma);
      rho_ma = std::move(rho_next);
      fmt::print("rho_ma: change {:.1e}\n", eps_rho);
    }
    core_sweeps_below_eps =
      (std::max(eps_c, eps_rho) < eps_target) ? core_sweeps_below_eps + 1 : 0;
    if (core_sweeps_below_eps == required_sweeps_below_eps)
      break;
  }

  // Print energy corrections after core convergence; refresh de_0 for valence
  if (max_it > 0) {
    std::cout << "\nAfter core convergence:\n";
    fmt::print("de_l({:4}): {:+12.9f} au  {:+12.5f} cm^-1\n", "core", de_c0,
               de_c0 * icm);
    for (std::size_t i = 0; i < valence.size(); ++i) {
      const auto de_v = MBPT::de_valence(valence[i], qk, lk, holes, excited);
      de_0[i] = de_v;
      fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n",
                 valence[i].shortSymbol(), de_v, de_v * icm);
    }
  }

  // Iterate for each required valence state; drop converged states each iter
  std::vector<DiracSpinor> unconverged_valence = valence;
  std::vector<int> valence_sweeps_below_eps(valence.size(), 0);
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "\nvalence it:" << it << "\n";
    MBPT::update_Lk_mnib(&lk_next, qk, excited, holes, unconverged_valence,
                         diagrams, &lk, a_damp, true, rho_ma_or_null);
    lk = lk_next;
    lk.write(Lk_file);
    for (std::size_t i = 0; i < valence.size(); ++i) {
      const auto &v = valence[i];
      const auto de_v = MBPT::de_valence(v, qk, lk, holes, excited);
      const auto eps_v = std::abs(2.0 * (de_v - de_0[i]) / (de_v + de_0[i]));
      fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n",
                 v.shortSymbol(), de_v, de_v * icm, eps_v);
      auto it2 =
        std::find(unconverged_valence.begin(), unconverged_valence.end(), v);
      if (it2 != unconverged_valence.end()) {
        de_0[i] = de_v;
        valence_sweeps_below_eps[i] =
          (eps_v < eps_target) ? valence_sweeps_below_eps[i] + 1 : 0;
        if (valence_sweeps_below_eps[i] == required_sweeps_below_eps) {
          // if converged: remove from list (speed up convergance)
          unconverged_valence.erase(it2);
        }
      }
    }
    if (unconverged_valence.empty())
      break;
  }
  std::cout << "\n";

  // With P: the valence energy correction de = E_W[x] - E_Q[x_F]. x_F is the
  // part of x that the Feynman direct diagram contains: the polarisation loop
  // (piece 1 of P) with the hole-particle interaction, x_F = P_F[g + x_F]. Its
  // direct (Q) closing uses the closing amplitude x' [calculate_Lk_vnab_F()].
  // x_F has no core single excitations: with R, the same x_F, and de has
  // <v|dV^R|v> added. Sigma_L is not formed with P: only the valence energy
  // corrections.
  if (include_P) {
    MBPT::LadderDiagrams diagrams_F{false, false, false, false};
    diagrams_F.P1 = true;
    diagrams_F.hole_particle = true;

    std::cout << "\nInitial fill of x_F table: core + valence\n" << std::flush;
    Coulomb::LkTable lk_F;
    MBPT::fill_Lk_mnib(&lk_F, qk, excited, holes, core_and_val, diagrams_F,
                       max_k, true);
    Coulomb::NkTable L_vnab_F;
    MBPT::add_Lk_vnab_F_entries(&L_vnab_F, excited, holes, valence,
                                lk_F.max_k());
    MBPT::calculate_Lk_vnab_F(&L_vnab_F, qk, excited, holes, valence,
                              diagrams_F, nullptr, nullptr, 0.0);
    auto lk_F_next = lk_F;
    auto L_vnab_F_next = L_vnab_F;

    std::vector<double> E_Q_F;
    std::cout << "\nInitial x_F corrections:\n";
    for (const auto &v : valence) {
      const auto E_Q_v =
        MBPT::de_valence_F(v, qk, lk_F, L_vnab_F, holes, excited);
      E_Q_F.push_back(E_Q_v);
      fmt::print("E_Q({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n", v.shortSymbol(),
                 E_Q_v, E_Q_v * icm);
    }

    // x_F (core and valence entries) and x' together, each from the previous
    // sweep: the core entries depend on no others, so the fixed point is the
    // same as converging the core first
    int sweeps_below_eps = 0;
    for (int it = 1; it <= max_it; ++it) {
      std::cout << "\nx_F it:" << it << "\n";
      MBPT::calculate_Lk_vnab_F(&L_vnab_F_next, qk, excited, holes, valence,
                                diagrams_F, &lk_F, &L_vnab_F, a_damp);
      MBPT::update_Lk_mnib(&lk_F_next, qk, excited, holes, core_and_val,
                           diagrams_F, &lk_F, a_damp, true);
      // for the damping, don't swap: tables and their copies must match
      lk_F = lk_F_next;
      L_vnab_F = L_vnab_F_next;
      double max_eps = 0.0;
      for (std::size_t i = 0; i < valence.size(); ++i) {
        const auto &v = valence[i];
        const auto E_Q_v =
          MBPT::de_valence_F(v, qk, lk_F, L_vnab_F, holes, excited);
        const auto eps_v =
          std::abs(2.0 * (E_Q_v - E_Q_F[i]) / (E_Q_v + E_Q_F[i]));
        E_Q_F[i] = E_Q_v;
        max_eps = std::max(max_eps, eps_v);
        fmt::print("E_Q({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n",
                   v.shortSymbol(), E_Q_v, E_Q_v * icm, eps_v);
      }
      sweeps_below_eps = (max_eps < eps_target) ? sweeps_below_eps + 1 : 0;
      if (sweeps_below_eps == required_sweeps_below_eps)
        break;
    }
    // x_F starts from the lowest order in every run, so may be unconverged
    // even when x was read in converged
    if (sweeps_below_eps < required_sweeps_below_eps) {
      fmt::print("\nWarning: x_F not converged in max_it = {} sweeps\n",
                 max_it);
    }

    // With R: the interaction with the change of the core density,
    // <v|dV^R|v>, from rho_ma (all orders, with its mirror image) and the
    // occupation blocks (third order); and its third-order value, with
    // rho_ma at lowest order
    std::vector<double> dV_R(valence.size(), 0.0);
    if (include_R) {
      const auto density = MBPT::relaxation_density(qk, holes, excited);
      const auto rho_ma_lowest = MBPT::core_singles(
        MBPT::core_singles_source(qk, excited, holes, nullptr), qk, excited,
        holes, nullptr);
      std::cout << "\nCore relaxation, <v|dV^R|v> (cm^-1): rho_ma part, "
                   "occupation part (third order), total; third order\n";
      fmt::print("{:4s}  {:11}  {:11}  {:11}  {:11}\n", "", "rho_ma",
                 "occupation", "dV^R", "dV^R(3)");
      for (std::size_t i = 0; i < valence.size(); ++i) {
        const auto &v = valence[i];
        const auto [rho_part, occupation_part] =
          MBPT::de_valence_relaxation(v, qk, rho_ma, density, holes, excited);
        const auto [rho_part_3, occupation_part_3] =
          MBPT::de_valence_relaxation(v, qk, rho_ma_lowest, density, holes,
                                      excited);
        dV_R[i] = rho_part + occupation_part;
        fmt::print("{:4s}  {:11.7}  {:11.7}  {:11.7}  {:11.7}\n",
                   v.shortSymbol(), rho_part * icm, occupation_part * icm,
                   dV_R[i] * icm, (rho_part_3 + occupation_part_3) * icm);
      }
      std::cout << "\nValence energy corrections, de = E_W[x] - E_Q[x_F] + "
                   "<v|dV^R|v> (cm^-1):\n";
      fmt::print("{:4s}  {:11}  {:11}  {:11}  {:11}\n", "", "E_W[x]",
                 "E_Q[x_F]", "dV^R", "de");
    } else {
      std::cout << "\nValence energy corrections, de = E_W[x] - E_Q[x_F] "
                   "(cm^-1):\n";
      fmt::print("{:4s}  {:11}  {:11}  {:11}\n", "", "E_W[x]", "E_Q[x_F]",
                 "de");
    }
    for (std::size_t i = 0; i < valence.size(); ++i) {
      const auto &v = valence[i];
      const auto E_W_x = MBPT::de_valence_w(v, qk, lk, holes, excited, &sjt);
      if (include_R) {
        fmt::print("{:4s}  {:11.7}  {:11.7}  {:11.7}  {:11.7}\n",
                   v.shortSymbol(), E_W_x * icm, E_Q_F[i] * icm, dV_R[i] * icm,
                   (E_W_x - E_Q_F[i] + dV_R[i]) * icm);
      } else {
        fmt::print("{:4s}  {:11.7}  {:11.7}  {:11.7}\n", v.shortSymbol(),
                   E_W_x * icm, E_Q_F[i] * icm, (E_W_x - E_Q_F[i]) * icm);
      }
    }
    std::cout << "\nSigma_L is not formed with include_P: no sl file written\n";
    return;
  }

  //----------------------------------------------------------------------------

  // Build all Sigma_L first (parallelisable), then print
  std::vector<MBPT::GMatrix> SigL_v;
  fmt::print("\nCalculating Sigma_L matrix (method: {}):\n",
             parseSigmaLMethod(method));
  fmt::print("Sigma_L sub-grid: r0={:.1e}, rmax={:.1f}, stride={}\n", sig_r0,
             sig_rmax, sig_stride);
  {
    IO::ChronoTimer t("", true);
    for (const auto &v : valence) {
      std::cout << v << "\n";
      if (method == SigmaLMethod::basis) {
        const auto ladder_v =
          MBPT::external_line_ladder(v, holes, excited, qk, lk, diagrams,
                                     max_it, eps_target, a_damp, true);
        SigL_v.push_back(MBPT::Sigma_ladder(v, ladder_v, holes, excited, sig_r0,
                                            sig_rmax, sig_stride, include_G));
      } else if (method == SigmaLMethod::ratio2) {
        SigL_v.push_back(MBPT::Sigma_ladder_ratio2(
          v, holes, excited, qk, lk, sig_r0, sig_rmax, sig_stride, include_G));
      } else {
        SigL_v.push_back(MBPT::Sigma_ladder(v, holes, excited, qk, &lk, sig_r0,
                                            sig_rmax, sig_stride, include_G));
      }
    }
  }

  // Mirror images of the ladder diagrams that the ladder equation does not
  // generate
  if (mirror_terms) {
    // Kept in memory only: the table is then not a full fill, so it is not
    // written to the Qk file [see fill_Qk_projection_states()]
    std::cout << "\nQk integrals of the projection states (mirror terms):\n";
    {
      IO::ChronoTimer t("", true);
      const auto num_added =
        MBPT::fill_Qk_projection_states(&qk, both, wf.basis(), valence, max_k);
      fmt::print("Calculated {} new Qk integrals\n", num_added);
    }
    std::cout << "\nMirror terms (lowest order), added to Sigma_L:\n";
    IO::ChronoTimer t("", true);
    for (std::size_t i = 0; i < valence.size(); ++i) {
      const auto &v = valence[i];
      // Projection states x: every basis state of kappa_v
      const auto x_orbitals = MBPT::projection_states(v, wf.basis());
      const auto Sigma_mirror = MBPT::Sigma_ladder_mirror(
        v, x_orbitals, holes, excited, qk, lk.max_k(), include_L4, sig_r0,
        sig_rmax, sig_stride, include_G);
      const auto de_mirror = v * (Sigma_mirror * v);
      fmt::print("de({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n", v.shortSymbol(),
                 de_mirror, de_mirror * icm);
      SigL_v[i] += Sigma_mirror;
    }
  }

  std::cout << "\nEnergy corrections:\n";
  fmt::print("{:4s}  {:13}  {:11}  {:10}  {:10}  {:10}\n", "", "HF", "de(2)",
             "de(l) [Q*L]", "de(l) [W*L]", "<v|Sig_L|v>");
  const auto Ec_HF = wf.coreEnergyHF();
  const auto de2_c = MBPT::de_core(qk, qk, holes, excited);
  const auto del_c = MBPT::de_core(qk, lk, holes, excited);
  const auto del_c_2 = MBPT::de_core_w(qk, lk, holes, excited, &sjt);
  fmt::print("{:4s}  {:13.7}  {:11.7}  {:10.7}  {:10.7}\n", "core", Ec_HF * icm,
             de2_c * icm, del_c * icm, del_c_2 * icm);

  for (std::size_t i = 0; i < valence.size(); ++i) {
    const auto &v = valence[i];
    const auto e0 = v.en();
    const auto de2 = MBPT::de_valence(v, qk, qk, holes, excited);
    const auto del = MBPT::de_valence(v, qk, lk, holes, excited);
    // Alternative: antisymmetrise the first (Coulomb) integral via W=Q+P,
    // rather than the second (ladder) integral via lk.P. Should match 'del'.
    const auto del2 = MBPT::de_valence_w(v, qk, lk, holes, excited, &sjt);
    const auto deS = v * (SigL_v[i] * v);
    fmt::print("{:4s}  {:13.7}  {:11.7}  {:10.7}  {:10.7}  {:10.7}\n",
               v.shortSymbol(), e0 * icm, de2 * icm, del * icm, del2 * icm,
               deS * icm);
  }

  // Hermiticity check of Sigma_L. The exact correlation potential, at fixed
  // energy, is Hermitian (real-symmetric kernel). Sigma_L as built (|W><L|)
  // is not, however exactly L is known: the ladder equation generates only
  // one member of each mirror (time-reversed) pair of diagrams in which a
  // core pair is excited before the valence electron interacts [closings
  // (a) of L2, (a+b) of L4, (c+d) of L1, (c) of L3]. The transpose of such a
  // diagram is its mirror image, which is not generated. The asymmetry here
  // measures those diagrams. With mirror_terms, the third-order mirror images
  // are added above, so the asymmetry is from higher orders, and from the
  // exchange closings of L2 and L3 (their mirror images are hole-particle
  // diagrams). The stored matrix is K_ij*dr_j (drj_in_place), so compare
  // w_i*S_ij with w_j*S_ji. Uses only the matrix as built: no new integrals
  // are evaluated.
  if (check_symmetry) {
    std::cout << "\nSigma_L hermiticity: eps = max|K - K^T| / max|K|\n";
    for (std::size_t iv = 0; iv < valence.size(); ++iv) {
      const auto &S = SigL_v[iv];
      const auto w = [&](std::size_t i) {
        return wf.grid().drdu(S.i0() + i * S.stride()) * wf.grid().du() *
               double(S.stride());
      };
      double max_K = 0.0, asym = 0.0;
      double max_g = 0.0, asym_g = 0.0;
      for (std::size_t i = 0; i < S.size(); ++i) {
        for (std::size_t j = 0; j < S.size(); ++j) {
          const auto Kij = w(i) * S.ff(i, j);
          const auto Kji = w(j) * S.ff(j, i);
          max_K = std::max(max_K, std::abs(Kij));
          asym = std::max(asym, std::abs(Kij - Kji));
          if (S.includes_g()) {
            // gg block symmetric; fg block the transpose of gf
            const auto gij = w(i) * S.gg(i, j);
            const auto gji = w(j) * S.gg(j, i);
            const auto xij = w(i) * S.fg(i, j);
            const auto xji = w(j) * S.gf(j, i);
            max_g = std::max({max_g, std::abs(gij), std::abs(xij)});
            asym_g =
              std::max({asym_g, std::abs(gij - gji), std::abs(xij - xji)});
          }
        }
      }
      fmt::print("{:4s}  ff: eps = {:.1e}", valence[iv].shortSymbol(),
                 max_K > 0.0 ? asym / max_K : 0.0);
      if (S.includes_g()) {
        fmt::print("   g: eps = {:.1e}", max_g > 0.0 ? asym_g / max_g : 0.0);
      }
      std::cout << "\n";
    }
  }

  // Symmetrise Sigma_L, S -> (S + S^T)/2: the correlation potential must be
  // Hermitian. Each diagram whose mirror image is not generated (see above)
  // is replaced by the average of itself and its mirror image, i.e., the
  // complete pair at half weight: consistent with the energy, which holds one
  // member of the pair. With mirror_terms, both members of each third-order
  // pair are already in Sigma_L at full weight, and so in the energy:
  // symmetrising leaves them unchanged. Expectation values are unchanged.
  // The stored matrix is K_ij*dr_j, so symmetrise the kernel:
  // S_ij -> (S_ij + (w_j/w_i) S_ji)/2. nb: the hermiticity check above
  // reports the raw (pre-symmetrisation) asymmetry.
  for (auto &S : SigL_v) {
    const auto w = [&](std::size_t i) {
      return wf.grid().drdu(S.i0() + i * S.stride()) * wf.grid().du() *
             double(S.stride());
    };
    for (std::size_t i = 0; i < S.size(); ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        const auto rji = w(j) / w(i);
        const auto ff_ij = S.ff(i, j);
        const auto ff_ji = S.ff(j, i);
        S.ff(i, j) = 0.5 * (ff_ij + rji * ff_ji);
        S.ff(j, i) = 0.5 * (ff_ji + ff_ij / rji);
        if (S.includes_g()) {
          const auto gg_ij = S.gg(i, j);
          const auto gg_ji = S.gg(j, i);
          S.gg(i, j) = 0.5 * (gg_ij + rji * gg_ji);
          S.gg(j, i) = 0.5 * (gg_ji + gg_ij / rji);
        }
      }
    }
    if (S.includes_g()) {
      // Off-diagonal spinor blocks pair across the diagonal:
      // K_fg(i,j) = K_gf(j,i). Each (i,j) handles one fg/gf pair.
      for (std::size_t i = 0; i < S.size(); ++i) {
        for (std::size_t j = 0; j < S.size(); ++j) {
          const auto rji = w(j) / w(i);
          const auto fg_ij = S.fg(i, j);
          const auto gf_ji = S.gf(j, i);
          S.fg(i, j) = 0.5 * (fg_ij + rji * gf_ji);
          S.gf(j, i) = 0.5 * (gf_ji + fg_ij / rji);
        }
      }
    }
  }

  // Write the Sigma_L matrices to disk
  std::vector<MBPT::SigmaLData> SLs;
  for (std::size_t i = 0; i < valence.size(); ++i) {
    const auto &v = valence[i];
    SLs.push_back({v.kappa(), v.n(), v.en(), std::move(SigL_v[i])});
  }
  MBPT::write_SigmaL(sl_file, SLs, wf.grid());
}

} // namespace MBPT
