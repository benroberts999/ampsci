#include "Modules/laddertests.hpp"
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

//******************************************************************************
//******************************************************************************
// Module for testing ladder diagram implementation
void ladder(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nLadder Module:\n\n";

  input.check({{"min", "lowest core n to include"},
               {"max", "maximum excited n to include"},
               {"max_l", "maximum excited l to include"},
               {"Qfile", "filename to read/write Qk integrals"},
               {"Lfile", "filename to read/write Qk integrals"},
               {"progbar", "Print progress bar? [true]"},
               {"max_it", "Max # iterations"},
               {"eps_target", "Target for convergance [1.0e-3]"}});

  // Example for retrieving input options:
  const auto min_n = input.get("min", 0);
  const auto max_n = input.get("max", 99);
  const auto max_l = input.get("max_l", 99);
  const auto max_it = input.get("max_it", 15);
  const auto eps_target = input.get("eps_target", 1.0e-3);
  std::cout << "min_n (core)    = " << min_n << "\n";
  std::cout << "max_n (excited) = " << max_n << "\n";
  std::cout << "max_l (excited) = " << max_l << "\n";
  std::cout << "max_it          = " << max_it << "\n";
  std::cout << "eps_target      = " << eps_target << "\n";

  const auto print_progbar = input.get("progbar", true);

  using namespace std::string_literals;
  const auto Qfname = input.get("Qfile", "qk.out"s);
  const auto Lfname = input.get("Lfile", "lk.out"s);

  // Sort basis into core/excited/valcne
  std::vector<DiracSpinor> core, excited, valence;
  const auto en_core = wf.en_coreval_gap();
  for (const auto &Fn : wf.basis) {
    if (Fn.en() < en_core && Fn.n >= min_n) {
      core.push_back(Fn);
    } else if (Fn.en() > en_core && Fn.n <= max_n && Fn.l() <= max_l) {
      excited.push_back(Fn);
    }
  }
  for (const auto &Fv : wf.valence) {
    // valence.push_back(Fv);
    // nb: use basis version of valence states
    const auto pFv = std::find(wf.basis.cbegin(), wf.basis.cend(), Fv);
    valence.push_back(*pFv);
  }

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

  // Calculate initial energy corrections:
  std::cout << "\nCore/Valence MBPT(2) shifts, using Yk table" << std::endl;
  std::cout << "Core: " << MBPT::de_core(yk, yk, core, excited) << "\n";
  //
  std::cout << "Valence, using wf.valence" << std::endl;
  for (const auto &v : wf.valence) {
    std::cout << v.symbol() << " " << MBPT::de_valence(v, yk, yk, core, excited)
              << "\n";
  }
  //
  std::cout << "Valence, using basis" << std::endl;
  for (const auto &v : valence) {
    std::cout << v.symbol() << " " << MBPT::de_valence(v, yk, yk, core, excited)
              << "\n";
  }

  std::cout << "\nFill Qk table:\n";
  Coulomb::QkTable qk;
  const auto ok = qk.read(Qfname);
  if (!ok) {
    qk.fill(both, yk);
    qk.write(Qfname);
  }

  std::cout << "\nCore/Valence MBPT(2) shifts, using Qk table" << std::endl;
  std::cout << "Core: " << MBPT::de_core(qk, qk, core, excited) << "\n";
  for (const auto &v : valence) {
    std::cout << v.symbol() << " " << MBPT::de_valence(v, qk, qk, core, excited)
              << "\n";
  }

  // Fill Lk table:
  Coulomb::LkTable lk;
  Coulomb::LkTable lk_next;
  // Each run will continue from where last left off
  const bool read_lad = lk.read(Lfname);
  std::cout << (read_lad ? "Read in existing ladder diagrams\n" :
                           "Calculating ladder diagrams from scratch\n");

  //----------------------------------------------------------------------------
  // Iterate Lk equations
  std::cout << "\nFilling Lk table: core + valence & iterate\n" << std::flush;
  double de_0 = MBPT::de_valence(valence.front(), qk, lk, core, excited);
  double de_c0 = MBPT::de_core(qk, lk, core, excited);
  bool core_converged = false;
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "it:" << it << "\n";
    {
      if (!core_converged) {
        // Don't update core terms if core energy shift converged?
        IO::ChronoTimer t("Fill Lk(core)");
        MBPT::calculate_Lk_mnib(&lk_next, qk, excited, core, core, &sjt, &lk,
                                print_progbar);
      }
      IO::ChronoTimer t("Fill Lk(valence)");
      MBPT::calculate_Lk_mnib(&lk_next, qk, excited, core, valence, &sjt, &lk,
                              print_progbar);
    }
    lk = lk_next; // XXX use swap or similar?
    lk.write(Lfname);

    // check convergance (core):
    const auto de_c = MBPT::de_core(qk, lk, core, excited);
    const auto eps_c = std::abs((de_c - de_c0) / de_c);
    de_c0 = de_c;
    std::cout << "de_l(core): " << de_c << " " << eps_c << "\n";
    if (eps_c < eps_target)
      core_converged = true;

    // check convergance (valence):
    const auto de_i = MBPT::de_valence(valence.front(), qk, lk, core, excited);
    const auto eps = std::abs((de_i - de_0) / de_i);
    de_0 = de_i;
    std::cout << "de_l(" << valence.front().shortSymbol() << "): " << de_i
              << " " << eps << "\n";
    if (eps < eps_target)
      break;
  }
  lk.summary();
  //----------------------------------------------------------------------------

  std::cout << "\nEnergy corrections:\n";
  std::cout << "       de(2)/au   de(l)/au    de(2)/cm^-1 "
               "  de(l)/cm^-1\n";
  for (const auto &v : valence) {

    auto de2 = MBPT::de_valence(v, qk, qk, core, excited);
    auto del = MBPT::de_valence(v, qk, lk, core, excited);

    printf("%4s: %11.8f %11.8f   %9.3f %9.3f\n", v.shortSymbol().c_str(), de2,
           del, de2 * PhysConst::Hartree_invcm, del * PhysConst::Hartree_invcm);
  }

  check_L_symmetry(core, excited, valence, qk, yk.SixJ(), &lk);
}

//******************************************************************************
void check_L_symmetry(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &valence,
                      const Coulomb::CoulombTable &qk,
                      const Angular::SixJTable &sj,
                      const Coulomb::CoulombTable *const lk) {

  std::cout << "\ncheck_L_symmetry\n";
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<std::size_t> e_index(0, excited.size() - 1);
  std::uniform_int_distribution<std::size_t> c_index(0, core.size() - 1);
  std::uniform_int_distribution<std::size_t> v_index(0, valence.size() - 1);

  double max_del = 0.0;
  std::string s1, s2;
  const int num_to_test = 15000;
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
      auto lkmnab = MBPT::Lkmnab(k, m, n, a, b, qk, core, excited, &sj);
      auto lknmba = MBPT::Lkmnab(k, n, m, b, a, qk, core, excited, &sj);

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
}

} // namespace Module
