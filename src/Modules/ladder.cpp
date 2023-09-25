#include "Modules/ladder.hpp"
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/Ladder.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include <numeric>
#include <random>

namespace Module {

//==============================================================================
//==============================================================================
// Module for testing ladder diagram implementation
void ladder(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nLadder Module:\n\n";

  input.check(
      {{"min", "lowest core n to include [0]"},
       {"max", "maximum excited n to include [99]"},
       {"max_l", "maximum excited l to include [99]"},
       {"max_k", "maximum k to include in Qk [99]"},
       {"include_L4", "Inlcude 4th Ladder diagram [false]"},
       {"fk", "List of doubles. Effective screening factors. Used to "
              "calculate Lk. []"},
       {"eta", "List of doubles. Effective hp factors. Only used to "
               "print energy shift. []"},
       {"Qfile", "filename to read/write Qk integrals"},
       {"form_Q", "Form or read Qk? (if have lk already, dont' need!) [true]"},
       {"Lfile", "filename to read/write Qk integrals"},
       {"progbar", "Print progress bar? [true]"},
       {"max_it", "Max # iterations [15]"},
       {"eps_target", "Target for convergance [1.0e-4]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Get input + print to screen
  const auto min_n = input.get("min", 0);
  const auto max_n = input.get("max", 99);
  const auto max_l = input.get("max_l", 99);
  const auto max_k = input.get("max_k", 99);
  const auto include_L4 = input.get("include_L4", false);
  const auto fk = input.get("fk", std::vector<double>{});
  const auto etak = input.get("eta", std::vector<double>{});
  const auto max_it = input.get("max_it", 15);
  const auto eps_target = input.get("eps_target", 1.0e-4);
  const auto print_progbar = input.get("progbar", true);
  std::cout << "min_n (core)    = " << min_n << "\n";
  std::cout << "max_n (excited) = " << max_n << "\n";
  std::cout << "max_l (excited) = " << max_l << "\n";
  std::cout << std::boolalpha;
  std::cout << "include_L4      = " << include_L4 << "\n";
  std::cout << "max_k           = " << max_k << "\n";
  std::cout << "max_it          = " << max_it << "\n";
  std::cout << "eps_target      = " << eps_target << "\n";

  if (fk.empty()) {
    std::cout << "No screening.\n";
  } else {
    std::cout << "Using screening factors, fk = ";
    std::for_each(fk.cbegin(), fk.cend(),
                  [](auto x) { std::cout << x << ", "; });
    std::cout << "1.0\n";
  }
  if (!etak.empty()) {
    std::cout << "Using eta (hp) factors, eta_k = ";
    std::for_each(etak.cbegin(), etak.cend(),
                  [](auto x) { std::cout << x << ", "; });
    std::cout << "1.0 (only in de, does not affect Lk)\n";
  }

  // Sort basis into core/excited/valcne
  std::vector<DiracSpinor> core, excited, valence;
  const auto en_core = wf.FermiLevel();
  for (const auto &Fn : wf.basis()) {
    if (Fn.en() < en_core && Fn.n() >= min_n) {
      core.push_back(Fn);
    } else if (Fn.en() > en_core && Fn.n() <= max_n && Fn.l() <= max_l) {
      excited.push_back(Fn);
    }
  }
  for (const auto &Fv : wf.valence()) {
    // nb: use _basis_ version of valence states (only for de, makes no diff)
    const auto pFv = std::find(wf.basis().cbegin(), wf.basis().cend(), Fv);
    valence.push_back(*pFv);
  }

  // in/out file names (default based on basis)
  using namespace std::string_literals;
  const auto name = wf.identity() + DiracSpinor::state_config(excited);
  const auto Qfname = input.get("Qfile", name + ".qk"s);
  const auto Lfname = input.get("Lfile", name + ".lk"s);

  // Create the "total" basis, which has core+excited, but only those states
  // actually included (i.e., [n_min, n_max]). This is used to calculate Qk.
  // Reduces size of Qk; should also reduce lookup time.
  // Use immediately invoked lambda
  const std::vector<DiracSpinor> both = [&]() {
    std::vector<DiracSpinor> t_basis = core;
    t_basis.insert(t_basis.end(), excited.cbegin(), excited.cend());
    return t_basis;
  }();

  // Fill Yk table (used to fill Qk table)
  const Coulomb::YkTable yk(both);
  // Store reference to Yk's 6J table (save typing):
  const auto &sjt = yk.SixJ();

  // Calculate initial energy corrections (mainly for checks):
  std::cout << "\nCore/Valence MBPT(2) shifts, using Yk table" << std::endl;
  std::cout << "Core: " << MBPT::de_core(yk, yk, core, excited) << "\n";
  //
  std::cout << "Valence, using wf.valence()" << std::endl;
  for (const auto &v : wf.valence()) {
    std::cout << v.symbol() << " " << MBPT::de_valence(v, yk, yk, core, excited)
              << "\n";
  }
  //
  std::cout << "Valence, using basis" << std::endl;
  for (const auto &v : valence) {
    std::cout << v.symbol() << " " << MBPT::de_valence(v, yk, yk, core, excited)
              << "\n";
  }
  if (!fk.empty() || !etak.empty()) {
    std::cout << "Valence, using basis + fk + eta_k" << std::endl;
    for (const auto &v : valence) {
      std::cout << v.symbol() << " "
                << MBPT::de_valence(v, yk, yk, core, excited, fk, etak) << "\n";
    }
  }

  // Form the Qk table.
  // If we already have L, and just want to print de, don't need to do this
  const auto formQ = input.get("form_Q", true);
  Coulomb::QkTable qk;
  if (formQ) {
    std::cout << "\nFill Qk table:\n";
    const auto ok = qk.read(Qfname);
    if (!ok) {
      qk.fill(both, yk, max_k);
      qk.write(Qfname);
    }
  }

  std::cout << "\nCore/Valence MBPT(2) shifts, using Qk table" << std::endl;
  std::cout << "Core: " << MBPT::de_core(qk, qk, core, excited) << "\n";
  for (const auto &v : valence) {
    std::cout << v.symbol() << " " << MBPT::de_valence(v, qk, qk, core, excited)
              << "\n";
  }
  std::cout << "\n";

  // Fill Lk table:
  Coulomb::LkTable lk;
  Coulomb::LkTable lk_next;
  // Each run will continue from where last left off
  const bool read_lad = lk.read(Lfname);
  std::cout << (read_lad ? "Re-starting using existing ladder diagrams\n" :
                           "Calculating ladder diagrams from scratch\n");

  //----------------------------------------------------------------------------
  // Iterate Lk equations
  std::cout << "\nFilling Lk table: core + valence & iterate\n" << std::flush;

  // initial corrections (will be zero if this is first run)
  std::vector<double> de_0(valence.size());
  for (std::size_t i = 0; i < valence.size(); ++i) {
    de_0[i] = MBPT::de_valence(valence[i], qk, lk, core, excited);
  }
  double de_c0 = MBPT::de_core(qk, lk, core, excited);

  bool core_converged = false;
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "it:" << it << "\n";
    {
      IO::ChronoTimer t("Fill Lk");
      if (!core_converged) {
        // Don't update core terms if core energy shift converged?
        // include screening for core parts?
        MBPT::fill_Lk_mnib(&lk_next, qk, excited, core, core, include_L4, sjt,
                           &lk, print_progbar, fk);
      }
      // in theory: each valence may converge differently, don't need to re-run
      // for valence states which already converged...
      MBPT::fill_Lk_mnib(&lk_next, qk, excited, core, valence, include_L4, sjt,
                         &lk, print_progbar, fk);
    }
    lk = lk_next; // XXX use swap or similar?
    lk.write(Lfname);

    // check convergance (core):
    const auto de_c = MBPT::de_core(qk, lk, core, excited);
    const auto eps_c = std::abs((de_c - de_c0) / de_c);
    de_c0 = de_c;
    std::cout << "de_l(core): ";
    printf("%10.7f %.1e\n", de_c, eps_c);
    if (eps_c < eps_target)
      core_converged = true;

    // check convergance (valence):
    double eps = 0.0;
    for (std::size_t i = 0; i < valence.size(); ++i) {
      // XXX include eta here?
      const auto de_v =
          MBPT::de_valence(valence[i], qk, lk, core, excited, fk, etak);
      const auto eps_v = std::abs((de_v - de_0[i]) / de_v);
      de_0[i] = de_v;
      std::cout << "de_l(" << valence[i].shortSymbol() << ") : ";
      printf("%10.7f %.1e\n", de_v, eps_v);
      if (eps_v > eps)
        eps = eps_v;
    }
    // keep going until core and all valence have converged!
    if (eps < eps_target && core_converged)
      break;
  }
  lk.summary();
  //----------------------------------------------------------------------------

  std::cout << "\nEnergy corrections:\n";
  std::cout << "       de(2)/au   de(A)/au   de(l)/au    de(2)/cm^-1 "
               "  de(A)/cm^-1   de(l)/cm^-1\n";
  for (const auto &v : valence) {

    const auto de2 = MBPT::de_valence(v, yk, yk, core, excited);
    const auto dea = MBPT::de_valence(v, yk, yk, core, excited, fk, etak) - de2;
    const auto del = MBPT::de_valence(v, yk, lk, core, excited, fk, etak);

    printf("%4s: %11.8f %11.8f %11.8f   %9.3f %9.3f %9.3f\n",
           v.shortSymbol().c_str(), de2, dea, del,
           de2 * PhysConst::Hartree_invcm, dea * PhysConst::Hartree_invcm,
           del * PhysConst::Hartree_invcm);
  }

  // check_L_symmetry(core, excited, valence, qk, include_L4, yk.SixJ(), &lk);

  // // Calculate Sigma_l, compare de_l to <v|Sigma_l|v>
  // bool include_G = false;
  // const auto sigp = MBPT::Sigma_params{MBPT::Method::Goldstone,
  //                                      min_n,
  //                                      include_G,
  //                                      0,
  //                                      false,
  //                                      false,
  //                                      0.0,
  //                                      0.0,
  //                                      1.0,
  //                                      false,
  //                                      false,
  //                                      "",
  //                                      fk,
  //                                      etak};

  // MBPT::GoldstoneSigma Sigma(wf.vHF(), wf.basis(), sigp,
  //                            MBPT::rgrid_params{1.0e-3, 30.0, 6}, "na");

  // std::cout << "\nEnergy corrections, using Sigma:\n";
  // for (const auto &v : valence) {
  //   const auto sig_l = Sigma.Sigma_l(v, yk, lk, core, excited);

  //   const auto SlFv = Sigma.Sigmal_Fv(v, yk, lk, core, excited);

  //   std::cout << v.symbol() << " "
  //             << v * Sigma.act_G_Fv(sig_l, v) * PhysConst::Hartree_invcm << " "
  //             << v * SlFv * PhysConst::Hartree_invcm << "\n";
  // }
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
