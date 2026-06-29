#include "MBPT/Ladder.hpp"
#include "Angular/include.hpp"
#include "CI/CI_Integrals.hpp"
#include "Coulomb/include.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include <numeric>
#include <random>

namespace Module {

// Declare, register, then define below.
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_ladder{
  "ladder", "Calculates ladder diagrams and energy corrections", &ladder};
} // namespace

// Helper, defined below.
void check_L_symmetry(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &valence,
                      const Coulomb::QkTable &qk, bool include_L4,
                      const Angular::SixJTable &sj,
                      const Coulomb::LkTable *const lk = nullptr);

//==============================================================================
//==============================================================================
// Module for testing ladder diagram implementation
void ladder(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"min_n_core", "lowest core n to include [1]"},
     {"basis", "Basis string specifying which states to include; must be a "
               "subset of full "
               "basis), e.g. '30spdf'. Default: include entire excited basis."},
     {"max_k", "maximum k to include in Qk. Put -1 to include all. [8]"},
     {"include_L4", "Inlcude 4th Ladder diagram [false]"},
     {"full_basis",
      "Construct the ladder correlation potential Sigma_L = sum_i L|i><i| "
      "using the full basis for {|i>}, otherwise just single HF |v> "
      "eigenstate. Using the full basis usually requires extending Qk. "
      "[false]"},
     {"Qk_file", "Filename for storing Qk Coulomb integrals. By default, is "
                 "<Identity>.qk. If 'false' will not read or write."},
     {"Lk_file", "Filename for storing Lk ladder integrals. By default, is "
                 "<Identity>.lk. If 'false' will not read or write."},
     {"from_scratch", "If true, don't read existing Qk/Lk files (still "
                      "writes). [false]"},
     {"max_it", "Max # iterations. If zero, will simply read ladder diagrams "
                "in (and then run a symmetry test). [15]"},
     {"damp",
      "Damping factor for iterations, [0,1). 0 means no damping. [0.0]"},
     {"eps_target", "Target for convergance [1.0e-4]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Get input + print to screen
  const auto min_n_core = input.get("min_n_core", 1);
  const auto basis_str = input.get<std::string>("basis", "");
  const auto max_k = input.get("max_k", 8);
  const auto include_L4 = input.get("include_L4", false);
  const auto full_basis = input.get("full_basis", false);

  const auto max_it = input.get("max_it", 15);
  const auto a_damp = input.get("damp", 0.0);
  const auto eps_target = input.get("eps_target", 1.0e-5);

  // Sort basis into core/excited/valence
  const auto en_core = wf.FermiLevel();
  const auto holes = qip::select_if(wf.basis(), [&](const auto &Fn) {
    return Fn.en() < en_core && Fn.n() >= min_n_core;
  });
  const auto excited =
    CI::basis_subset(wf.basis(), basis_str, wf.coreConfiguration());

  // nb: use _basis_ version of valence states in iterations
  std::vector<DiracSpinor> valence;
  for (const auto &Fv : wf.valence()) {
    const auto pFv = std::find(wf.basis().cbegin(), wf.basis().cend(), Fv);
    if (pFv == wf.basis().cend()) {
      std::cout << "Warning: Basis missing valence state: " << Fv << "\n";
      continue;
    }
    // const auto orth = *pFv * Fv - 1.0;
    // const auto eps = (pFv->en() - Fv.en()) / Fv.en();
    // fmt::print("{:3s} orth = {:.1e}, dE/E = {:.1e}\n", Fv.shortSymbol(), orth,
    //            eps);
    valence.push_back(*pFv);
  }

  std::cout << "\n";
  std::cout << "basis        = " << DiracSpinor::state_config(excited) << "\n";
  std::cout << "min_n (core) = " << min_n_core << "\n";
  std::cout << std::boolalpha;
  std::cout << "include_L4   = " << include_L4 << "\n";
  std::cout << "full_basis   = " << full_basis << "\n";
  std::cout << "max_k        = " << max_k << "\n";
  std::cout << "max_it       = " << max_it << "\n";
  std::cout << "damp         = " << a_damp << "\n";
  std::cout << "eps_target   = " << eps_target << "\n";

  // in/out file names (default based on atom identity)
  using namespace std::string_literals;
  const auto ident = wf.identity();
  const auto Qk_file = input.get("Qk_file", ident + ".qk.abf"s);
  const auto lk4 = include_L4 ? "_l4" : ""s;
  const auto Lk_file =
    input.get<std::string>("Lk_file", ident + lk4 + ".lk.abf"s);
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
  if (!from_scratch)
    qk.read(Qk_file);
  qk.fill(both, yk, max_k);
  qk.write(Qk_file);

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
  const bool did_read_Lk = !from_scratch ? lk_next.read(Lk_file) : false;
  std::cout << (did_read_Lk ? "\nRe-starting using existing ladder diagrams\n" :
                              "\nCalculating ladder diagrams from scratch\n");

  //----------------------------------------------------------------------------
  // Iterate Lk equations

  // "First iteration": fill_Lk_mnib keeps existing (read) entries and adds only
  // newly-required diagrams (or fills from scratch if not read). Subsequent
  // iterations below re-iterate the existing table.
  std::cout << "\nInitial fill of Lk table: core + valence\n" << std::flush;
  MBPT::fill_Lk_mnib(&lk_next, qk, excited, holes, core_and_val, include_L4,
                     sjt, max_k, true);
  lk_next.write(Lk_file);
  lk = lk_next;

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
  std::cout << "\nAfter core convergence:\n";
  fmt::print("de_l({:4}): {:+12.9f} au  {:+12.5f} cm^-1\n", "core", de_c0,
             de_c0 * icm);
  for (std::size_t i = 0; i < valence.size(); ++i) {
    const auto de_v = MBPT::de_valence(valence[i], qk, lk, holes, excited);
    de_0[i] = de_v;
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n",
               valence[i].shortSymbol(), de_v, de_v * icm);
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

  // return;

  //----------------------------------------------------------------------------

  // Extend Qk table to the FULL basis. With full_basis, Sigma_ladder projects
  // the ladder onto the full basis of each kappa_v, so Lkmnij() needs Qk
  // integrals involving those (high-n) projection states - which are absent if
  // Qk was filled only over the [n_min,n_max] subset 'both'. Not always
  // strictly required (only when projection states fall outside 'both'), but we
  // just refill over the full basis to be safe. Skipped for single-state |v>.
  if (full_basis) {
    std::cout << "\nExtending Qk table to full basis (for Sigma_ladder):\n"
              << std::flush;
    const Coulomb::YkTable yk_full(wf.basis());
    qk.fill(wf.basis(), yk_full, max_k);
    qk.write(Qk_file);
  }

  // Build all Sigma_L first (parallelisable), then print
  std::vector<MBPT::GMatrix> SigL_v;
  fmt::print("\nCalculating Sigma_L matrix (using {} projection):\n",
             full_basis ? "full basis" : "single state");
  for (const auto &v : valence) {
    std::cout << v << "\n";
    const std::vector<DiracSpinor> proj_v{v};
    const auto &proj = full_basis ? wf.basis() : proj_v;
    SigL_v.push_back(MBPT::Sigma_ladder(v.kappa(), v.en(), holes, excited, proj,
                                        qk, &lk, sjt, include_L4));
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

  // Use this to check the symmetries
  if (max_it == 0)
    check_L_symmetry(holes, excited, valence, qk, include_L4, sjt, &lk);
}

//==============================================================================
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

    const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, a, b);
    for (int k = k0; k <= kI; k += 2) {
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

} // namespace Module
