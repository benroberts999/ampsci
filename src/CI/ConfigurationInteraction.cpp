#include "ConfigurationInteraction.hpp"
#include "Angular/Angular.hpp"
#include "CI_Integrals.hpp"
#include "CSF.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/Matrix.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/Sigma2.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include "qip/Vector.hpp"
#include <array>
#include <fstream>
#include <vector>

namespace CI {
std::vector<PsiJPi> configuration_interaction(const IO::InputBlock &input,
                                              const Wavefunction &wf) {

  // Check input options:
  input.check(
      {{"ci_basis",
        "Basis used for CI expansion; must be a sub-set of full ampsci basis "
        "[default: 10spdf]"},
       {"J", "List of total angular momentum J for CI solutions (comma "
             "separated). Must be integers (two-electron only). [0]"},
       {"J+", "As above, but for EVEN CSFs only (takes precedence over J)."},
       {"J-", "As above, but for ODD CSFs (takes precedence over J)."},
       {"num_solutions", "Number of CI solutions to find (for each J/pi) [5]"},
       {"sigma1", "Include one-body MBPT correlations? [false]"},
       {"sigma2", "Include two-body MBPT correlations? [false]"},
       {"cis2_basis",
        "The subset of ci_basis for which the two-body MBPT corrections are "
        "calculated. Must be a subset of ci_basis. If existing sk file has "
        "more integrals, they will be used. [default: Nspdf, where N is "
        "maximum n for core + 3]"},
       {"s1_basis",
        "Basis used for the one-body MBPT diagrams (Sigma^1). These are the "
        "most important, so in general the default (all basis states) should "
        "be used. Must be a subset of full ampsci basis. [default: full "
        "basis]\n"
        " - Note: if CorrelationPotential is available, it will be used "
        "instead of calculating the Sigma_1 integrals"},
       {"s2_basis",
        "Basis used for internal lines of the two-body MBPT diagrams "
        "(Sigma^2). Must be a subset of s1_basis. [default: s1_basis]"},
       {"n_min_core", "Minimum n for core to be included in MBPT [1]"},
       {"max_k",
        "Maximum k (multipolarity) to include when calculating new "
        "Coulomb integrals. Higher k often contribute negligably. Note: if qk "
        "file already has higher-k terms, they will be included. Set negative "
        "(or very large) to include all k. [6]"},
       {"qk_file",
        "Filename for storing two-body Coulomb integrals. By default, is "
        "At.qk, where At is atomic symbol."},
       {"sk_file",
        "Filename for storing two-body Sigma_2 integrals. By default, is "
        "At_n_b_k.sk, where At is atomic symbol, n is n_min_core, b is "
        "cis2_basis, k is max_k."},
       {"exclude_wrong_parity_box",
        "Excludes the Sigma_2 box corrections that "
        "have 'wrong' parity when calculating Sigma2 matrix elements. Note: If "
        "existing sk file already has these, they will be included [false]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return {};
  }

  //----------------------------------------------------------------------------
  // Single-particle basis:
  std::cout << "\nConstruct single-particle basis:\n";

  // Determine the sub-set of basis to use in CI:
  const auto basis_string = input.get("ci_basis", std::string{"10spdf"});

  // Select from wf.basis() [MBPT basis], those which match input 'basis_string'
  // exclude those in coreConfiguration
  const std::vector<DiracSpinor> ci_sp_basis =
      CI::basis_subset(wf.basis(), basis_string, wf.coreConfiguration());

  // Print info re: basis to screen:
  std::cout << "\nUsing " << DiracSpinor::state_config(ci_sp_basis) << " = "
            << ci_sp_basis.size() << " orbitals in CI expansion\n";

  //----------------------------------------------------------------------------

  // Determine different basis subsets

  // check if including MBPT corrections
  const auto include_Sigma1 = input.get("sigma1", false);
  const auto include_Sigma2 = input.get("sigma2", false);
  const auto max_k_Coulomb = input.get("max_k", 6);
  const auto exclude_wrong_parity_box =
      input.get("exclude_wrong_parity_box", false);
  const auto include_MBPT = include_Sigma1 || include_Sigma2;

  // s1 and s2 MBPT basis
  const auto s1_basis_string = input.get("s1_basis");
  const auto &s1_basis = s1_basis_string ?
                             CI::basis_subset(wf.basis(), *s1_basis_string) :
                             wf.basis();
  const auto s2_basis_string = input.get("s2_basis");
  const auto &s2_basis = s2_basis_string ?
                             CI::basis_subset(wf.basis(), *s2_basis_string) :
                             s1_basis;

  // Ensure s2_basis is subset of s1_basis
  assert(s2_basis.size() <= s1_basis.size() &&
         "s2_basis must be a subset of s1_basis");

  // Split basis' into core/excited (for MBPT evaluations)
  const auto n_min_core = input.get("n_min_core", 1);
  const auto [core_s1, excited_s1] =
      MBPT::split_basis(s1_basis, wf.FermiLevel(), n_min_core);
  const auto [core_s2, excited_s2] =
      MBPT::split_basis(s2_basis, wf.FermiLevel(), n_min_core);

  // S2 corrections are included only for this subset of the CI basis:
  const auto Ncore = DiracSpinor::max_n(wf.core()) + 3;
  const auto cis2_basis_string =
      input.get("cis2_basis", std::to_string(Ncore) + "spdf");
  const auto &cis2_basis = CI::basis_subset(ci_sp_basis, cis2_basis_string);

  //----------------------------------------------------------------------------

  if (include_MBPT) {
    std::cout << "\nIncluding MBPT: "
              << (include_Sigma1 && include_Sigma2 ? "Σ_1 + Σ_2" :
                  include_Sigma1                   ? "Σ_1" :
                                                     "Σ_2")
              << "\n";
    std::cout << "Including core excitations from n ≥ " << n_min_core << "\n";
    if (max_k_Coulomb >= 0 && max_k_Coulomb < 50) {
      std::cout << "Including k ≤ " << max_k_Coulomb
                << " in Coulomb integrals (unless already calculated)\n";
    }
    if (include_Sigma1) {
      if (wf.Sigma()) {
        std::cout << "With existing Correlation Potential for Σ_1:\n";
        wf.Sigma()->print_info();
      } else {
        std::cout << "With basis for Σ_1: "
                  << DiracSpinor::state_config(s1_basis) << "\n";
      }
    }
    if (include_Sigma2) {
      std::cout << "With basis for Σ_2: " << DiracSpinor::state_config(s2_basis)
                << "\n";

      std::cout << "Including Σ_2 correction to Coulomb integrals up to: "
                << DiracSpinor::state_config(cis2_basis) << "\n";
      if (exclude_wrong_parity_box) {
        std::cout
            << "Excluding the Σ_2 diagrams that have the 'wrong' parity\n";
        std::cout << "(Unless they were already calculated in sk file)\n";
      }
    }
    std::cout << "\n";
  }

  //----------------------------------------------------------------------------

  // Lookup table; stores all qk's
  Coulomb::QkTable qk;
  {
    std::cout << "Calculate two-body Coulomb integrals: Q^k_abcd\n";

    const auto qk_filename = input.get("qk_file", wf.atomicSymbol() + ".qk");

    // Try to read from disk (may already have calculated Qk)
    qk.read(qk_filename);
    const auto existing = qk.count();
    {

      // Try to limit number of Coulomb integrals we calculate
      // use whole basis (these are used inside Sigma_2)
      // If not including MBPT, only need to caculate smaller set of integrals

      // First, calculate the integrals between ci basis states:
      {
        std::cout << "For: " << DiracSpinor::state_config(ci_sp_basis) << "\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(ci_sp_basis);
        qk.fill(ci_sp_basis, yk, max_k_Coulomb, false);
      }

      // Selection function for which Qk's to calculate.
      // For Sigma, we only need those with 1 or 2 core electrons
      // i.e., no Q_vwxy, Q_vabc, or Q_abcd
      // Note: we *do* need Q_vwxy for the CI part (but with smaller basis)
      const auto select_Q_sigma =
          [eF = wf.FermiLevel()](int, const DiracSpinor &s,
                                 const DiracSpinor &t, const DiracSpinor &u,
                                 const DiracSpinor &v) {
            // Only calculate Coulomb integrals with 1 or 2 electrons in the core
            auto num = MBPT::number_below_Fermi(s, t, u, v, eF);
            return num == 1 || num == 2;
          };

      // Then, add those required for Sigma_1 (unless we have matrix!)
      if (include_Sigma1 && !wf.Sigma()) {
        const auto temp_basis = qip::merge(core_s1, excited_s1);
        std::cout << "and: " << DiracSpinor::state_config(temp_basis) << "\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(temp_basis);
        qk.fill_if(temp_basis, yk, select_Q_sigma, max_k_Coulomb, false);
      }

      // Then, add those required for Sigma_2 (unless we already did Sigma_1)
      if (include_Sigma2 && !(include_Sigma1 && !wf.Sigma())) {
        const auto temp_basis = qip::merge(core_s2, excited_s2);
        std::cout << "and: " << DiracSpinor::state_config(temp_basis) << "\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(temp_basis);
        qk.fill_if(temp_basis, yk, select_Q_sigma, max_k_Coulomb, false);
      }

      // print summary
      qk.summary();

      // If we calculated new integrals, write to disk
      const auto total = qk.count();
      assert(total >= existing);
      const auto new_integrals = total - existing;
      std::cout << "Calculated " << new_integrals << " new Coulomb integrals\n";
      if (new_integrals > 0) {
        qk.write(qk_filename);
      }
    }
    std::cout << "\n" << std::flush;
  }

  //----------------------------------------------------------------------------

  // Create lookup table for one-particle matrix elements, h1
  const auto h1 =
      wf.Sigma() ?
          CI::calculate_h1_table(ci_sp_basis, *wf.Sigma(), include_Sigma1) :
          CI::calculate_h1_table(ci_sp_basis, core_s1, excited_s1, qk,
                                 include_Sigma1);

  //----------------------------------------------------------------------------
  // Calculate MBPT corrections to two-body Coulomb integrals

  // Here, write basis info into filename, since these are _internal_ lines!
  const auto Sk_filename = input.get(
      "sk_file", wf.atomicSymbol() + "_" + std::to_string(n_min_core) + "_" +
                     DiracSpinor::state_config(excited_s2) +
                     (max_k_Coulomb >= 0 && max_k_Coulomb < 50 ?
                          "_" + std::to_string(max_k_Coulomb) :
                          "") +
                     ".sk");

  Coulomb::LkTable Sk;
  if (include_Sigma2) {
    std::cout << "Calculate two-body MBPT integrals: Σ^k_abcd\n";

    std::cout << "For: " << DiracSpinor::state_config(cis2_basis) << ", using "
              << DiracSpinor::state_config(excited_s2) << "\n";

    Sk = CI::calculate_Sk(Sk_filename, cis2_basis, core_s2, excited_s2, qk,
                          max_k_Coulomb, exclude_wrong_parity_box);
    std::cout << "\n" << std::flush;
  }

  //----------------------------------------------------------------------------
  const auto J_list = input.get("J", std::vector<int>{0});
  const auto J_even_list = input.get("J+", J_list);
  const auto J_odd_list = input.get("J-", J_list);
  const auto num_solutions = input.get("num_solutions", 5);
  const auto ci_input = input.get("ci_input", std::string{""});

  std::vector<PsiJPi> levels;

  // even parity:
  for (const auto J : J_even_list) {
    const auto t_levels = run_CI(ci_sp_basis, int(std::round(2 * J)), +1,
                                 num_solutions, h1, qk, Sk, include_Sigma2);
    levels.push_back(t_levels);
  }

  // odd parity:
  for (const auto J : J_odd_list) {
    const auto t_levels = run_CI(ci_sp_basis, int(std::round(2 * J)), -1,
                                 num_solutions, h1, qk, Sk, include_Sigma2);
    levels.push_back(t_levels);
  }

  // const auto sort_ouput = input.get("sort", false);
  const auto compare_energy = [](const auto &a, const auto &b) {
    return a.energy(0) < b.energy(0);
  };

  // Find minimum (ground-state) energy (for level comparison)
  // Assumes levels for each J/Pi are sorted (they always are)
  const auto e0 =
      !levels.empty() ?
          std::min_element(levels.cbegin(), levels.cend(), compare_energy)
              ->energy(0) :
          0.0;

  std::cout << "\nLevel Summary:\n\n";
  std::cout << "config Term J  π  #   Energy(au)  Energy(/cm)   Level(/cm) "
               " gJ\n";
  for (const auto &Psis : levels) {

    std::string pi = Psis.parity() == 1 ? "" : "°";
    const auto twoj = Psis.twoJ();

    for (std::size_t i = 0; i < Psis.num_solutions(); ++i) {

      const auto [config, gJ, L, twoS] = Psis.info(i);

      fmt::print("{:<6s} {}{:<2s} {:>2} {:+2} {:2}  {:+11.8f}  {:+11.2f}  "
                 "{:11.2f}",
                 config, int(std::round(twoS + 1)),
                 AtomData::L_symbol((int)std::round(L)) + pi, twoj / 2,
                 Psis.parity(), i, Psis.energy(i),
                 Psis.energy(i) * PhysConst::Hartree_invcm,
                 (Psis.energy(i) - e0) * PhysConst::Hartree_invcm);
      if (gJ != 0.0) {
        fmt::print("  {:.4f}", gJ);
      }
      std::cout << "\n";
    }
  }

  return levels;
}

//==============================================================================
//==============================================================================
//==============================================================================
PsiJPi run_CI(const std::vector<DiracSpinor> &ci_sp_basis, int twoJ, int parity,
              int num_solutions, const Coulomb::meTable<double> &h1,
              const Coulomb::QkTable &qk, const Coulomb::LkTable &Sk,
              bool include_Sigma2) {

  auto printJ = [](int twoj) {
    return twoj % 2 == 0 ? std::to_string(twoj / 2) :
                           std::to_string(twoj) + "/2";
  };
  auto printPi = [](int pi) { return pi > 0 ? "even" : "odd"; };

  fmt::print("Run CI for J={}, {} parity\n", printJ(twoJ), printPi(parity));
  std::cout << std::flush;

  PsiJPi psi{twoJ, parity, ci_sp_basis};

  if (twoJ < 0) {
    std::cout << "Fail: twoJ must >=0\n";
    return psi;
  }
  if (twoJ % 2 != 0) {
    std::cout << "Fail: twoJ must be even for two-electron CSF\n";
    return psi;
  }

  const auto N_CSFs = psi.CSFs().size();
  std::cout << "Total CSFs: " << psi.CSFs().size() << "\n";
  std::cout << std::flush;

  //----------------------------------------------------------------------------

  // Construct the CI matrix:
  const auto Hci = include_Sigma2 ? CI::construct_Hci(psi, h1, qk, &Sk) :
                                    CI::construct_Hci(psi, h1, qk);

  //----------------------------------------------------------------------------

  fmt::print("Find first {} solutions\n", num_solutions);
  {
    IO::ChronoTimer t2("Diagonalise");
    psi.solve(Hci, num_solutions);
  }
  const auto E0 = psi.energy(0);

  // For calculating g-factors
  // (Use non-rel formula? Or relativistic M1?)
  DiracOperator::M1nr m1{};
  DiracOperator::M1 m1_rel{ci_sp_basis.front().grid(), PhysConst::alpha, 0.0};

  // const auto m1_tab = ExternalField::me_table(ci_sp_basis, &m1_rel);
  const auto m1_tab = ExternalField::me_table(ci_sp_basis, &m1_rel);

  int l1{-1}, l2{-1};

  for (std::size_t i = 0; i < N_CSFs && int(i) < num_solutions; ++i) {

    fmt::print(
        "{:<2} {} {:+1}  {:+11.8f} au  {:+11.2f} cm^-1  {:11.2f} cm^-1\n", i,
        0.5 * twoJ, parity, psi.energy(i),
        psi.energy(i) * PhysConst::Hartree_invcm,
        (psi.energy(i) - E0) * PhysConst::Hartree_invcm);

    const double minimum_percentage = 1.0; // min % to print
    std::size_t max_j = 0;
    double max_cj = 0.0;
    for (std::size_t j = 0ul; j < N_CSFs; ++j) {
      const auto cj = 100.0 * std::pow(psi.coef(i, j), 2);
      if (cj > max_cj) {
        max_cj = cj;
        max_j = j;
        l1 = Angular::nkindex_to_l(psi.CSF(j).state(0));
        l2 = Angular::nkindex_to_l(psi.CSF(j).state(1));
      }
      if (cj > minimum_percentage) {
        fmt::print("   {:<6s} {:5.3f}%\n", psi.CSF(j).config(true), cj);
      }
    }

    // g_J <JJz|J|JJz> = <JJz|L + 2*S|JJz>
    // take J=Jz, <JJz|J|JJz> = J
    // then: g_J = <JJ|L + 2*S|JJ> / J
    // And: <JJ|L + 2*S|JJ> = 3js * <A||L+2S||A> (W.E. Theorem)
    // const auto m1AA_NR = CI::ReducedME(psi.coefs(i), psi.CSFs(), twoJ, &m1);
    const auto m1AA_R =
        CI::ReducedME(psi.coefs(i), psi.CSFs(), twoJ, psi.coefs(i), psi.CSFs(),
                      twoJ, m1_tab, m1.rank(), m1.parity());
    const auto tjs = Angular::threej_2(twoJ, twoJ, 2, twoJ, -twoJ, 0);

    // Calculate g-factors, for line identification. Only defined for J!=0
    // const double gJnr = twoJ != 0 ? tjs * m1AA_NR / (0.5 * twoJ) : 0.0;
    // const double gJ = twoJ != 0 ? tjs * m1AA_R / (0.5 * twoJ) : 0.0;

    const double gJ = twoJ != 0 ? tjs * m1AA_R / (0.5 * twoJ) : 0.0;
    const double gJnr = gJ;

    // Determine Term Symbol, from g-factor
    // Use non-relativistic M1 operator for closest match
    const auto [S, L] = CI::Term_S_L(l1, l2, twoJ, gJnr);

    std::cout << "   --------------\n";
    if (twoJ != 0) {
      std::cout << "   gJ = " << gJ << "\n";
    }

    const auto tSp1 = 2.0 * S + 1.0;
    const auto config = psi.CSF(max_j).config();
    const auto pm = parity == 1 ? "" : "°";

    fmt::print("   {:<6s} {}^{}{}_{}\n", config, int(std::round(tSp1)),
               AtomData::L_symbol((int)std::round(L)), pm, twoJ / 2);

    psi.update_config_info(i, {config, gJ, double(L), 2.0 * S});

    std::cout << "\n";
  }

  return psi;
}

} // namespace CI