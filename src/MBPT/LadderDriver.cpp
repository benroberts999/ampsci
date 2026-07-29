#include "MBPT/LadderDriver.hpp"
#include "Angular/include.hpp"
#include "CI/CI_Integrals.hpp"
#include "Coulomb/include.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/Ladder.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace MBPT {

namespace {
// Helper, defined below.
void check_L_symmetry(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &valence,
                      const Coulomb::QkTable &qk, bool include_L4,
                      const Angular::SixJTable &sj,
                      const Coulomb::LkTable *const lk = nullptr);
} // namespace

//==============================================================================
// Driver for ladder diagram calculation: Ladder{} input block
void ladder(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"min_n_core", "lowest core n to include [1]"},
     {"basis", "Basis string specifying which states to include; must be a "
               "subset of full "
               "basis), e.g. '30spdf'. Default: include entire excited basis."},
     {"max_k", "maximum k to include in Qk. Put -1 to include all. [8]"},
     {"include_L4", "Inlcude 4th Ladder diagram [false]"},
     {"projection",
      "Projection basis {|i>} for the ladder correlation potential Sigma_L = "
      "sum_i L|i><i|. Options:\n"
      "   - single : just the single HF |v> eigenstate\n"
      "   - ladder : states in the ladder basis (see basis option)\n"
      "   - full   : entire basis; requires extending Qk (slow)\n"
      "   - ratio  : no projection; rescale each Sigma(2) term by L/Q "
      "(fast)\n"
      "   - direct : no projection; open the external line exactly (ladder "
      "vertex)\n"
      "  [ratio]"},
     {"each_valence",
      "If true, iterate L and calculate Sigma_L for each valence state "
      "separately. If false, only for the lowest valence state of each kappa "
      "(Correlations applies it to all states of that kappa). [false]"},
     {"Qk_file", "Filename for storing Qk Coulomb integrals. By default, is "
                 "<Identity>.qk. If 'false' will not read or write."},
     {"Lk_file", "Filename for storing Lk ladder integrals. By default, is "
                 "<Identity>.lk. If 'false' will not read or write."},
     {"sl_file", "Filename for storing the ladder correlation potential, "
                 "Sigma_L. By default, is <Identity>.sl. If 'false' will not "
                 "write. Read in via Correlations{ladder_file=...;}"},
     {"from_scratch", "If true, don't read existing Qk/Lk files (still "
                      "writes). [false]"},
     {"max_it", "Max # iterations. If zero, will simply read ladder diagrams "
                "in (any _new_ ladder diagrams will be calculated,). [15]"},
     {"damp",
      "Damping factor for iterations, [0,1). 0 means no damping. [0.0]"},
     {"eps_target", "Target for convergance [1.0e-5]"},
     {"rmin", "minimum radius of Sigma_L sub-grid [1.0e-4]"},
     {"rmax", "maximum radius of Sigma_L sub-grid [30.0]"},
     {"stride", "Only calculate Sigma_L every <stride> points. Default such "
                "that there are 150 points between (1e-4, 30)"},
     {"include_G", "Inlcude lower g-part into Sigma_L [false]"},
     {"check_symmetry",
      "Run the (slow) numerical symmetry test of the ladder integrals at the "
      "end [false]"}});
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
  const auto include_L4 = input.get("include_L4", false);

  // Method for Sigma_L: projection (single |v>, ladder basis, full basis),
  // Dzuba (no projection: rescale Sigma(2) terms by L/Q), or direct (no
  // projection: open the external line exactly)
  using namespace std::string_literals;
  const auto method = parseSigmaLMethod(input.get("projection", "ratio"s));

  const auto each_valence = input.get("each_valence", false);

  const auto max_it = input.get("max_it", 15);
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
  const auto en_core = wf.FermiLevel();
  const auto holes = qip::select_if(wf.basis(), [&](const auto &Fn) {
    return Fn.en() < en_core && Fn.n() >= min_n_core;
  });
  const auto excited =
    CI::basis_subset(wf.basis(), basis_str, wf.coreConfiguration());

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
  std::cout << "projection   = " << parseSigmaLMethod(method) << "\n";
  std::cout << "each_valence = " << each_valence << "\n";
  std::cout << "include_G    = " << include_G << "\n";
  std::cout << "max_k        = " << max_k << "\n";
  std::cout << "max_it       = " << max_it << "\n";
  std::cout << "damp         = " << a_damp << "\n";
  std::cout << "eps_target   = " << eps_target << "\n";

  // in/out file names (default based on atom identity)
  const auto ident = wf.identity();
  const auto Qk_file = input.get("Qk_file", ident + ".qk.abf"s);
  const auto lk4 = include_L4 ? "_l4" : ""s;
  const auto Lk_file =
    input.get<std::string>("Lk_file", ident + lk4 + ".lk.abf"s);
  const auto sl_ext =
    ".sl"s + (include_L4 ? "4" : "") + (include_G ? "g" : "") + ".abf"s;
  const auto sl_file = input.get<std::string>("sl_file", ident + sl_ext);
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
  MBPT::fill_Lk_mnib(&lk, qk, excited, holes, core_and_val, include_L4, sjt,
                     max_k, true);
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

  // Iterate for core states:
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "\ncore it:" << it << "\n";
    MBPT::update_Lk_mnib(&lk_next, qk, excited, holes, holes, include_L4, sjt,
                         &lk, a_damp, true);
    // for the damping, don't swap; lk and lk_prev must match before next iter
    lk = lk_next;
    lk.write(Lk_file);
    const auto de_c = MBPT::de_core(qk, lk, holes, excited);
    const auto eps_c = std::abs(2.0 * (de_c - de_c0) / (de_c + de_c0));
    de_c0 = de_c;
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n", "core",
               de_c, de_c * icm, eps_c);
    if (eps_c < eps_target)
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
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "\nvalence it:" << it << "\n";
    MBPT::update_Lk_mnib(&lk_next, qk, excited, holes, unconverged_valence,
                         include_L4, sjt, &lk, a_damp, true);
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
        if (eps_v < eps_target) {
          // if converged: remove from list (speed up convergance)
          unconverged_valence.erase(it2);
        }
      }
    }
    if (unconverged_valence.empty())
      break;
  }
  std::cout << "\n";

  //----------------------------------------------------------------------------

  // Extend Qk table for the projection states. With projection=full,
  // Sigma_ladder projects the ladder onto the full basis of each kappa_v, so
  // Lkmnij() needs Qk integrals involving those (high-n) projection states -
  // absent if Qk was filled only over the [n_min,n_max] subset 'both'.
  // Only basis states with a valence kappa are used in the projection, so
  // extend over both + those (not the entire basis). Not required for
  // single-state or ladder-basis projection (those states already in 'both').
  if (method == SigmaLMethod::full) {
    std::cout
      << "\nExtending Qk table for projection states (for Sigma_ladder):\n"
      << std::flush;

    // Basis states with a valence kappa, not already in 'both':
    const auto proj_states = qip::select_if(wf.basis(), [&](const auto &Fi) {
      const auto vkappa =
        std::any_of(valence.cbegin(), valence.cend(),
                    [&](const auto &v) { return v.kappa() == Fi.kappa(); });
      return vkappa && std::find(both.cbegin(), both.cend(), Fi) == both.cend();
    });

    std::cout << "Adding: " << DiracSpinor::state_config(proj_states) << "\n";
    const auto ext_basis = qip::merge(both, proj_states);
    const Coulomb::YkTable yk_ext(ext_basis);
    qk.fill(ext_basis, yk_ext, max_k);
    qk.write(Qk_file);
  }

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
      if (method == SigmaLMethod::direct) {
        SigL_v.push_back(MBPT::Sigma_ladder_direct(
          v, holes, excited, qk, yk, &lk, sjt, include_L4, sig_r0, sig_rmax,
          sig_stride, include_G));
        continue;
      }
      // Projection basis for Sigma_ladder; empty => ratio method
      const std::vector<DiracSpinor> proj_v{v};
      const std::vector<DiracSpinor> empty{};
      const auto &proj = method == SigmaLMethod::single ? proj_v :
                         method == SigmaLMethod::full   ? wf.basis() :
                         method == SigmaLMethod::ratio  ? empty :
                                                          excited;
      SigL_v.push_back(MBPT::Sigma_ladder(v, holes, excited, proj, qk, &lk, sjt,
                                          include_L4, sig_r0, sig_rmax,
                                          sig_stride, include_G));
    }
  }

  std::cout << "\nEnergy corrections:\n";
  fmt::print("{:4s}  {:13}  {:11}  {:10}  {:10}  {:10}\n", "", "HF", "de(2)",
             "de(l) [Q*L]", "de(l) [W*L]", "<v|Sig_L|v>");
  const auto Ec_HF = wf.coreEnergyHF();
  const auto de2_c = MBPT::de_core(qk, qk, holes, excited);
  const auto del_c = MBPT::de_core(qk, lk, holes, excited);
  const auto del_c_2 = MBPT::de_core(lk, qk, holes, excited);
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

  // Write the Sigma_L matrices to disk
  std::vector<MBPT::SigmaLData> SLs;
  for (std::size_t i = 0; i < valence.size(); ++i) {
    const auto &v = valence[i];
    SLs.push_back({v.kappa(), v.n(), v.en(), std::move(SigL_v[i])});
  }
  MBPT::write_SigmaL(sl_file, SLs, wf.grid());

  // Use this to check the symmetries
  if (check_symmetry)
    check_L_symmetry(holes, excited, valence, qk, include_L4, sjt, &lk);
}

//==============================================================================
namespace {
void check_L_symmetry(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &valence,
                      const Coulomb::QkTable &qk, bool include_L4,
                      const Angular::SixJTable &sj,
                      const Coulomb::LkTable *const lk) {
  std::cout << "\ncheck_L_symmetry\n";
  std::random_device rd;
  std::mt19937 gen(rd());
  if (excited.empty() || core.empty() || valence.empty())
    return;
  std::uniform_int_distribution<std::size_t> e_index(0, excited.size() - 1);
  std::uniform_int_distribution<std::size_t> c_index(0, core.size() - 1);
  std::uniform_int_distribution<std::size_t> v_index(0, valence.size() - 1);

  const int num_to_test = 150000;

  // worst-case trackers; init to -1 so first element always sets a result
  double max1 = -1.0, max2 = -1.0, max3 = -1.0;
  std::string label1, label2, label3;
  std::string info1, info2, info3;
  int count = 0;

  for (int tries = 0; tries < num_to_test; ++tries) {
    const auto &m = excited[e_index(gen)];
    const auto &n = excited[e_index(gen)];

    // half the time: core, other, valence
    const auto &a = tries % 2 == 0 ? core[c_index(gen)] : valence[v_index(gen)];

    const auto &b = core[c_index(gen)];
    auto sym = [](const auto &x) { return x.shortSymbol(); };

    // Full ladder k range: both parities (see k_minmax_L)
    const auto [k0, kI] = k_minmax_L(m, n, a, b);
    for (int k = k0; k <= kI; ++k) {
      ++count;
      auto Qkmnab = qk.Q(k, m, n, a, b);
      auto lkmnab =
        MBPT::Lkmnij(k, m, n, a, b, qk, core, excited, include_L4, sj);
      auto lknmba =
        MBPT::Lkmnij(k, n, m, b, a, qk, core, excited, include_L4, sj);

      auto lkmnab_tab = lk ? lk->Q(k, m, n, a, b) : 0.0;
      auto lknmba_tab = lk ? lk->Q(k, n, m, b, a) : 0.0;

      auto lbl = "L^" + std::to_string(k) + "_(" + sym(m) + sym(n) + sym(a) +
                 sym(b) + ")";

      // 1: Lk_mnab vs Lk_nmba
      auto d1 = std::abs(lkmnab - lknmba);
      if (d1 > max1) {
        max1 = d1;
        label1 = lbl;
        info1 = fmt::format("Qk     = {:11.4e}\n"
                            "Lk_mnab= {:12.5e}\n"
                            "Lk_nmba= {:12.5e}\n"
                            "del    = {:8.1e}",
                            Qkmnab, lkmnab, lknmba, d1);
      }

      if (lk) {
        // 2: Lk_mnab(T) vs Lk_nmba(T)
        auto d2 = std::abs(lkmnab_tab - lknmba_tab);
        if (d2 > max2) {
          max2 = d2;
          label2 = lbl;
          info2 = fmt::format("Qk        = {:11.4e}\n"
                              "Lk_mnab(T)= {:12.5e}\n"
                              "Lk_nmba(T)= {:12.5e}\n"
                              "del       = {:8.1e}",
                              Qkmnab, lkmnab_tab, lknmba_tab, d2);
        }

        // 3: Lk_mnab vs Lk_mnab(T)
        auto d3 = std::abs(lkmnab - lkmnab_tab);
        if (d3 > max3) {
          max3 = d3;
          label3 = lbl;
          info3 = fmt::format("Qk        = {:11.4e}\n"
                              "Lk_mnab   = {:12.5e}\n"
                              "Lk_mnab(T)= {:12.5e}\n"
                              "del       = {:8.1e}",
                              Qkmnab, lkmnab, lkmnab_tab, d3);
        }
      }
    }
  }

  std::cout << "Checked " << count << " Lk elements (worst case shown)\n\n";

  std::cout << "Lk_mnab vs Lk_nmba:\n  " << label1 << "\n" << info1 << "\n\n";

  if (lk) {
    std::cout << "Lk_mnab(T) vs Lk_nmba(T) [should be trivially the same]:\n  "
              << label2 << "\n"
              << info2 << "\n\n";
    std::cout << "Lk_mnab vs Lk_mnab(T) [direct vs table; differs if >1 "
                 "iteration]:\n  "
              << label3 << "\n"
              << info3 << "\n";
  } else {
    std::cout << "(no table -- skipping table comparisons)\n";
  }
}
} // namespace

} // namespace MBPT
