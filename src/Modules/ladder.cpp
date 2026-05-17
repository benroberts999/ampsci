#include "MBPT/Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
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
     {"max_n_excited",
      "maximum excited n to include. Default is to include entire basis."},
     {"max_l",
      "maximum excited l to include. Default is to include entire basis."},
     {"max_k", "maximum k to include in Qk. Put -1 to include all. [8]"},
     {"include_L4", "Inlcude 4th Ladder diagram [false]"},
     {"read", "true/false. If false, will not attemp to read Qk or Lk file "
              "from disk, and will start calculation from scratch. "
              "May still overwrite existing file. [true]"},
     {"max_it", "Max # iterations. If zero, will simply read ladder diagrams "
                "in (and then run a symmetry test). [15]"},
     {"damp",
      "Damping factor for iterations, [0,1). 0 means no damping. [0.35]"},
     {"eps_target", "Target for convergance [1.0e-5]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Get input + print to screen
  const auto min_n_core = input.get("min_n_core", 1);
  const auto max_n = input.get("max_n_excited", 999);
  const auto max_l = input.get("max_l", 999);
  const auto max_k = input.get("max_k", 8);
  const auto include_L4 = input.get("include_L4", false);

  const auto readQ = input.get("read", true);

  const auto max_it = input.get("max_it", 15);
  const auto a_damp = input.get("damp", 0.35);
  const auto eps_target = input.get("eps_target", 1.0e-5);

  std::cout << "min_n (core)    = " << min_n_core << "\n";
  std::cout << "max_n (excited) = " << max_n << "\n";
  std::cout << "max_l (excited) = " << max_l << "\n";
  std::cout << std::boolalpha;
  std::cout << "include_L4      = " << include_L4 << "\n";
  std::cout << "max_k           = " << max_k << "\n";
  std::cout << "max_it          = " << max_it << "\n";
  std::cout << "damp            = " << a_damp << "\n";
  std::cout << "eps_target      = " << eps_target << "\n";

  // Sort basis into core/excited/valcne
  std::vector<DiracSpinor> core, excited, valence;
  const auto en_core = wf.FermiLevel();
  for (const auto &Fn : wf.basis()) {
    if (Fn.en() < en_core && Fn.n() >= min_n_core) {
      core.push_back(Fn);
    } else if (Fn.en() > en_core && Fn.n() <= max_n && Fn.l() <= max_l) {
      excited.push_back(Fn);
    }
  }

  // nb: use _basis_ version of valence states in iterations
  std::cout << "\nValence vs. basis orthog:\n";
  for (const auto &Fv : wf.valence()) {
    const auto pFv = std::find(wf.basis().cbegin(), wf.basis().cend(), Fv);
    if (pFv == wf.basis().cend()) {
      std::cout << "Warning: Basis missing valence state: " << Fv << "\n";
      continue;
    }
    const auto orth = *pFv * Fv - 1.0;
    const auto eps = (pFv->en() - Fv.en()) / Fv.en();
    fmt::print("{:3s} orth = {:.1e}, dE/E = {:.1e}\n", Fv.shortSymbol(), orth,
               eps);
    valence.push_back(*pFv);
  }

  // Create the "total" basis, which has core+excited, but only those states
  // actually included (i.e., [n_min, n_max]). This is used to calculate Qk.
  // Reduces size of Qk; should also reduce lookup time.
  const auto both = qip::merge(core, excited);
  // the "i" orbitals: core + valence
  const auto core_and_val = qip::merge(core, valence);

  // in/out file names (default based on basis)
  using namespace std::string_literals;
  const auto ident = wf.identity();
  const auto basis_string = DiracSpinor::state_config(excited);
  const auto Qfname = ident + ".qk.abf"s;
  const auto lk4 = include_L4 ? "_l4" : ""s;
  const auto Lfname = ident + "_" + basis_string + lk4 + ".lk.abf"s;

  // Fill Yk table (used to fill Qk table)
  const Coulomb::YkTable yk(both);
  // Store reference to Yk's 6J table (save typing):
  const auto &sjt = yk.SixJ();

  // Form the Qk table.
  Coulomb::QkTable qk;
  std::cout << "\nFill (or read) Qk table:\n";
  const auto read_Qk = readQ ? qk.read(Qfname) : false;
  // nb: will not extend Qk integrals. Maybe not best choice?
  if (!read_Qk) {
    qk.fill(both, yk, max_k);
    qk.write(Qfname);
  } else {
    qk.summary();
    std::cout
      << "\nNote: no new Qk integrals will be calculated, even if basis "
         "has changed!\n";
  }

  std::cout << "\nMBPT(2) energy shifts, using HF vs. spline legs" << std::endl;
  fmt::print("{:<5s} {:>13s} {:>13s} {:>13s}   {}\n", "state", "de(HF)",
             "de(basis)", "de(Qk)", "eps");
  const auto dec_1 = MBPT::de_core(yk, yk, wf.core(), excited);
  const auto dec_2 = MBPT::de_core(yk, yk, core, excited);
  const auto dec_3 = MBPT::de_core(qk, qk, core, excited);
  const auto dec_eps = std::abs(2.0 * (dec_1 - dec_2) / (dec_1 + dec_2));
  fmt::print("{:<5s} {:+13.7f} {:+13.7f} {:+13.7f}   {:.0e}\n", "core", dec_1,
             dec_2, dec_3, dec_eps);
  for (const auto &sv : valence) {
    const auto pv = wf.getState(sv.n(), sv.kappa());
    assert(pv != nullptr); //valence is subset of wf.valence()
    const auto de_1 = MBPT::de_valence(*pv, yk, yk, core, excited);
    const auto de_2 = MBPT::de_valence(sv, yk, yk, core, excited);
    const auto de_3 = MBPT::de_valence(sv, qk, qk, core, excited);
    const auto de_eps = std::abs(2.0 * (de_1 - de_2) / (de_1 + de_2));
    fmt::print("{:<5s} {:+13.7f} {:+13.7f} {:+13.7f}   {:.0e}\n",
               sv.shortSymbol(), de_1, de_2, de_3, de_eps);
  }
  std::cout << "\n";

  // Fill Lk table:
  Coulomb::LkTable lk;
  Coulomb::LkTable lk_next;
  // Each run will continue from where last left off
  const bool read_Lk = readQ ? lk.read(Lfname) : false;
  std::cout << (read_Lk ? "\nRe-starting using existing ladder diagrams\n" :
                          "\nCalculating ladder diagrams from scratch\n");

  //----------------------------------------------------------------------------
  // Iterate Lk equations

  if (read_Lk) {
    lk.summary();
    // must have correct dimension!
    lk_next = lk;
  } else {
    std::cout << "\nInitial fill of Lk table: core + valence\n" << std::flush;
    MBPT::fill_Lk_mnib(&lk_next, qk, excited, core, core_and_val, include_L4,
                       sjt, true);
    lk_next.write(Lfname);
    lk = lk_next;
  }

  // convert to inverse cm
  const auto icm = PhysConst::Hartree_invcm;

  // initial corrections
  double de_c0 = MBPT::de_core(qk, lk, core, excited);
  std::vector<double> de_0;
  std::cout << "\nInitial Ladder Corrections:\n";
  fmt::print("de_l({:4}): {:+12.9f} au  {:+12.5f} cm^-1\n", "core", de_c0,
             de_c0 * icm);
  for (const auto &v : valence) {
    const auto de_v = MBPT::de_valence(v, qk, lk, core, excited);
    de_0.push_back(de_v);
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n", v.shortSymbol(),
               de_v, de_v * icm);
  }

  bool core_converged = false;
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "\nit:" << it << "\n";

    // nb: if core or particular valence state has converged
    // then we _could_ speed this up by skipping...
    // Also: in converge core first, valence probably faster
    // Would result in fewer integrals to calculate
    // Achieve by re-introducing {i_orbs}....
    MBPT::update_Lk_mnib(&lk_next, qk, excited, core, {}, include_L4, sjt, &lk,
                         a_damp, true, {});

    // std::swap(lk, lk_next); // lk never empty
    // for the damping, don't want to swap
    // lk and lk_prev must be the same before next iteration begins
    lk = lk_next;
    lk.write(Lfname);
    std::cout << "\n";

    // check convergance (core):
    const auto de_c = MBPT::de_core(qk, lk, core, excited);
    const auto eps_c = std::abs(2.0 * (de_c - de_c0) / (de_c + de_c0));
    de_c0 = de_c;
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n", "core",
               de_c, de_c * icm, eps_c);
    if (eps_c < eps_target)
      core_converged = true;

    // check convergance (valence):
    double eps = 0.0;
    for (std::size_t i = 0; i < valence.size(); ++i) {
      const auto &v = valence[i];
      const auto de_v = MBPT::de_valence(v, qk, lk, core, excited);
      const auto eps_v = std::abs(2.0 * (de_v - de_0[i]) / (de_v + de_0[i]));
      de_0[i] = de_v;
      fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n",
                 v.shortSymbol(), de_v, de_v * icm, eps_v);
      if (eps_v > eps)
        eps = eps_v;
    }
    // keep going until core and all valence have converged!
    if (eps < eps_target && core_converged)
      break;
  }
  std::cout << "\n";

  //----------------------------------------------------------------------------

  std::cout << "Energy corrections:\n";
  fmt::print("{:4s}  {:13}  {:11}  {:10}  {:10}\n", "", "HF", "de(2)",
             "de(l) [Q*L]", "de(l) [L*Q]");
  const auto Ec_HF = wf.coreEnergyHF();
  const auto de2_c = MBPT::de_core(qk, qk, core, excited);
  const auto del_c = MBPT::de_core(qk, lk, core, excited);
  const auto del_c_2 = MBPT::de_core(lk, qk, core, excited);
  fmt::print("{:4s}  {:13.7}  {:11.7}  {:10.7}  {:10.7}\n", "core", Ec_HF * icm,
             de2_c * icm, del_c * icm, del_c_2 * icm);

  for (const auto &v : valence) {

    const auto e0 = v.en();
    const auto de2 = MBPT::de_valence(v, qk, qk, core, excited);
    const auto del = MBPT::de_valence(v, qk, lk, core, excited);
    const auto del2 = MBPT::de_valence(v, lk, qk, core, excited);
    fmt::print("{:4s}  {:13.7}  {:11.7}  {:10.7}  {:10.7}\n", v.shortSymbol(),
               e0 * icm, de2 * icm, del * icm, del2 * icm);
  }

  // Use this to check the symmetries
  // old code that hasn't been tested in a while
  // Doing odd things?
  if (max_it == 0)
    check_L_symmetry(core, excited, valence, qk, include_L4, sjt, &lk);
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

  double max_del = 0.0;
  std::string s1, s2;
  const int num_to_test = 150000;
  for (int tries = 0; tries < num_to_test; ++tries) {
    const auto &m = excited[e_index(gen)];
    const auto &n = excited[e_index(gen)];

    // half the time: core, other, valence
    const auto &a = tries % 2 == 0 ? core[c_index(gen)] : valence[v_index(gen)];

    const auto &b = core[c_index(gen)];
    auto sym = [](const auto &x) { return x.shortSymbol(); };

    const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, a, b);
    for (int k = k0; k <= kI; k += 2) {
      auto gkmnab = qk.Q(k, m, n, a, b);
      auto lkmnab =
        MBPT::Lkmnij(k, m, n, a, b, qk, core, excited, include_L4, sj);
      auto lknmba =
        MBPT::Lkmnij(k, n, m, b, a, qk, core, excited, include_L4, sj);

      auto lkmnab_tab = lk ? lk->Q(k, m, n, a, b) : 0.0;
      auto lknmba_tab = lk ? lk->Q(k, n, m, b, a) : 0.0;

      auto del1 = std::abs(lkmnab - lknmba);
      auto del_tab = std::abs(lkmnab_tab - lknmba_tab);
      auto del2 = lk ? std::abs(lkmnab_tab - lkmnab) : 0.0;
      auto del = del1 + del_tab + del2;

      if (del > max_del) {
        max_del = del;
        s1 = qip::fstring("g         =%11.4e\n"
                          "Lk_mnab   =%12.5e\n"
                          "Lk_nmba   =%12.5e\n"
                          "Lk_mnab(T)=%12.5e\n"
                          "Lk_nmba(T)=%12.5e\n"
                          "del       =%8.1e",
                          gkmnab, lkmnab, lknmba, lkmnab_tab, lknmba_tab, del);
        s2 = "L^" + std::to_string(k) + "_(" + sym(m) + sym(n) + sym(a) +
             sym(b) + ")";
      }
    }
  }
  std::cout << s1 << "\n";
  std::cout << s2 << "\n";
  std::cout
    << "nb: Lk_mnab is calc'd directly, Lk_mnab(T) is from table - so "
       "won't be the same if more than 1 iteration has been perormed.\n";
}

} // namespace Module
