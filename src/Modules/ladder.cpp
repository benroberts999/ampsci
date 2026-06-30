#include "MBPT/Ladder.hpp"
#include "Angular/include.hpp"
#include "CI/CI_Integrals.hpp"
#include "Coulomb/include.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace Module {

// Declare, register, then define below.
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_ladder{"ladder",
                        "Checks ladder-diagram symmetries (reads Lk/Qk files; "
                        "see the Ladder{} block to compute them)",
                        &ladder};
} // namespace

// Helper, defined below.
void check_L_symmetry(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &valence,
                      const Coulomb::QkTable &qk, bool include_L4,
                      const Angular::SixJTable &sj,
                      const Coulomb::LkTable *const lk = nullptr);

//==============================================================================
// Module for testing/validating the ladder-diagram implementation.
// Production runs use the Ladder{} main block; this module only checks the
// numerical symmetries of the ladder integrals.
void ladder(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"min_n_core", "lowest core n to include [1]"},
     {"basis", "Basis string specifying which states to include; must be a "
               "subset of full "
               "basis), e.g. '30spdf'. Default: include entire excited basis."},
     {"max_k", "maximum k to include in Qk. Put -1 to include all. [8]"},
     {"include_L4", "Inlcude 4th Ladder diagram [false]"},
     {"Qk_file", "Filename for storing Qk Coulomb integrals. By default, is "
                 "<Identity>.qk. If 'false' will not read or write."},
     {"Lk_file", "Filename for storing Lk ladder integrals. By default, is "
                 "<Identity>.lk. If 'false' will not read or write."}});
  if (input.has_option("help")) {
    return;
  }

  const auto min_n_core = input.get("min_n_core", 1);
  const auto basis_str = input.get<std::string>("basis", "");
  const auto max_k = input.get("max_k", 8);
  const auto include_L4 = input.get("include_L4", false);

  using namespace std::string_literals;
  const auto ident = wf.identity();
  const auto Qk_file = input.get("Qk_file", ident + ".qk.abf"s);
  const auto lk4 = include_L4 ? "_l4"s : ""s;
  const auto Lk_file =
    input.get<std::string>("Lk_file", ident + lk4 + ".lk.abf"s);

  // Sort basis into core (holes) and excited sets
  const auto en_core = wf.FermiLevel();
  const auto holes = qip::select_if(wf.basis(), [&](const auto &Fn) {
    return Fn.en() < en_core && Fn.n() >= min_n_core;
  });
  const auto excited =
    CI::basis_subset(wf.basis(), basis_str, wf.coreConfiguration());

  // basis-version valence states
  std::vector<DiracSpinor> valence;
  for (const auto &Fv : wf.valence()) {
    const auto pFv = std::find(wf.basis().cbegin(), wf.basis().cend(), Fv);
    if (pFv != wf.basis().cend())
      valence.push_back(*pFv);
  }

  // Yk table provides the 6j table (and fills Qk if not read from file)
  const auto both = qip::merge(holes, excited);
  const Coulomb::YkTable yk(both);
  const auto &sjt = yk.SixJ();

  // Read Qk (rebuild if absent) and Lk
  Coulomb::QkTable qk;
  Coulomb::LkTable lk;
  const bool got_qk = qk.read(Qk_file);
  const bool got_lk = lk.read(Lk_file);
  if (!got_qk) {
    std::cout << "No Qk file (" << Qk_file
              << "); building from basis subset.\n";
    qk.fill(both, yk, max_k);
  }
  if (!got_lk) {
    std::cout << "No Lk file (" << Lk_file
              << "); checking direct (untabulated) symmetries only.\n";
  }

  check_L_symmetry(holes, excited, valence, qk, include_L4, sjt,
                   got_lk ? &lk : nullptr);
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
