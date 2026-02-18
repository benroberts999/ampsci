#include "ConfigurationInteraction.hpp"
#include "Angular/include.hpp"
#include "CI_Integrals.hpp"
#include "CSF.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
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
#include <iostream>
#include <vector>

namespace CI {
std::vector<PsiJPi> configuration_interaction(const IO::InputBlock &input,
                                              const Wavefunction &wf) {

  // Check input options:
  input.check(
    {{"ci_basis",
      "Basis used for CI expansion; must be a sub-set of full ampsci basis "
      "[default: 20spdf]"},
     {"J", "List of total angular momentum J for CI solutions (comma "
           "separated). Must be integers (two-electron only). []"},
     {"J+", "As above, but for EVEN CSFs only (takes precedence over J)."},
     {"J-", "As above, but for ODD CSFs (takes precedence over J)."},
     {"num_solutions", "Number of CI solutions to find (for each J/pi) [5]"},
     {"sigma1", "Include one-body MBPT correlations? [false]"},
     {"sigma2", "Include two-body MBPT correlations? [false]"},
     {"Brueckner", "Use Brueckner (spectrum) states for CI basis? Must have "
                   "Spectrum and sigma1. [false]"},
     {"cis2_basis",
      "The subset of ci_basis for which the two-body MBPT corrections are "
      "calculated. Must be a subset of ci_basis. If existing sk file has "
      "more integrals, they will be used. [default: Nspdf, where N is "
      "maximum n for core + 3]"},
     {"Breit2", "Include two-body Breit? Default is true if Breit included in "
                "HF. Ignored if Breit not included in HF. [true]"},
     {"Breit_basis",
      "Subset of ci_basis used to include two-body Breit "
      "corrections into CI matrix. Large basis is slow, uses "
      "huge memory, and makes small contribution. [default: Nspdf, where N is "
      "maximum n for core + 6]"},
     {"s1_basis",
      "Usually should be left as default. Basis used for the one-body MBPT "
      "diagrams (Sigma^1) internal lines. These are the "
      "most important, so in general the default (all basis states) should "
      "be used. Must be a subset of full ampsci basis. [default: full "
      "basis]\n"
      " - Note: if CorrelationPotential is available, it will be used "
      "instead of calculating the Sigma_1 integrals"},
     {"s2_basis",
      "Usually should be left blank. Basis used for internal lines of the "
      "two-body MBPT diagrams "
      "(Sigma^2) internal lines. Must be a subset of s1_basis. [default: "
      "s1_basis]"},
     {"n_min_core", "Minimum n for core to be included in MBPT [1]"},
     {"max_k",
      "Maximum k (multipolarity) to include when calculating new "
      "Coulomb integrals. Higher k often contribute negligably. Note: if qk "
      "file already has higher-k terms, they will be included. Set negative "
      "(or very large) to include all k. [8]"},
     {"qk_file",
      "Filename for storing two-body Coulomb integrals. By default, is "
      "~ At.qk, where At is atomic symbol + 'identity'."},
     {"sk_file",
      "Filename for storing two-body Sigma_2 integrals. By default, is "
      "At_n_b_k.sk, where At is atomic symbol, n is n_min_core, b is "
      "cis2_basis, k is max_k."},
     {"bk_file",
      "Filename for storing two-body Breit integrals. By default, is "
      "~ At.bk, where At is atomic symbol + 'identity'."},
     {"no_new_integrals",
      "Usually false. If set to true, ampsci will not calculate any new "
      "Coulomb or Sigma_2 integrals, even if they are implied by the above "
      "settings. This saves time when we know all required integrals already "
      "exist, since the code doesn't need to check. [true]"},
     {"exclude_wrong_parity_box",
      "Excludes the Sigma_2 box corrections that "
      "have 'wrong' parity when calculating Sigma2 matrix elements. Note: If "
      "existing sk file already has these, they will be included [false]"},
     {"sort_output", "Sort output by energy? Default is to sort by J and Pi "
                     "first. [false]"},
     {"parallel_ci", "Run CI in parallel (solve each J/Pi in parallel). "
                     "Faster, uses slightly more memory [true]"}});

  // construct first, for RVO
  std::vector<PsiJPi> levels;

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return levels;
  }

  //----------------------------------------------------------------------------
  // Single-particle basis:
  std::cout << "\nConstruct single-particle basis:\n";

  // Determine the sub-set of basis to use in CI:
  const auto basis_string = input.get("ci_basis", std::string{"20spdf"});
  // options to include MBPT
  const auto include_Sigma1 = input.get("sigma1", false);
  const auto include_Sigma2 = input.get("sigma2", false);

  // Use Breuckner states for MBPT
  // nb: currently also use these Bruckner states in place of the "core"
  // I think this is the best option, due to orthogonality
  // These are not eigenstates of the V^HF potential used in the core
  // However, we don't explicitely use this for <v|h1|w>
  // We assume eignestates, so <v|h1|w> = E_w \delta_vw
  // (When not using Brueckner there is also <v|Sigma1|w>, which is not diag)
  // Perhaps a slight inconsistancy when including RPA into matrix elements?
  const auto Brueckner_raw = input.get("Brueckner", false);
  const auto Brueckner =
    Brueckner_raw && include_Sigma1 && !wf.spectrum().empty();

  const auto &t_basis = Brueckner ? wf.spectrum() : wf.basis();

  // maximum n present in core: used for default basis
  const auto N_max_core = DiracSpinor::max_n(wf.core());

  // Select from basis those which match input 'basis_string'
  // exclude those in coreConfiguration
  const std::vector<DiracSpinor> ci_sp_basis =
    CI::basis_subset(t_basis, basis_string, wf.coreConfiguration());

  // Print info re: basis to screen:
  std::cout << "\nUsing " << DiracSpinor::state_config(ci_sp_basis) << " = "
            << ci_sp_basis.size() << " orbitals in CI expansion\n";

  if (Brueckner) {
    std::cout
      << "CI + MBPT + Brueckner method: using Brueckner states for CI basis\n";
  }
  if (Brueckner_raw && !Brueckner) {
    fmt2::warning();
    std::cout << "Requested Brueckner method, but conditions (Spectrum and "
                 "Sigma1) not met\n";
  }

  // aditional output string for br method:
  using namespace std::string_literals;
  const auto br_string = Brueckner ? "_bru"s : "";

  //----------------------------------------------------------------------------

  // Determine different basis subsets

  // Details of MBPT
  const auto max_k_Coulomb = input.get("max_k", 8);
  const auto exclude_wrong_parity_box =
    input.get("exclude_wrong_parity_box", false);
  const auto include_MBPT = include_Sigma1 || include_Sigma2;

  // s1 and s2 MBPT basis
  const auto s1_basis_string = input.get("s1_basis");
  const auto &s1_basis =
    s1_basis_string ? CI::basis_subset(t_basis, *s1_basis_string) : t_basis;
  const auto s2_basis_string = input.get("s2_basis");
  const auto &s2_basis =
    s2_basis_string ? CI::basis_subset(t_basis, *s2_basis_string) : s1_basis;

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
  const auto cis2_basis_string =
    input.get("cis2_basis", std::to_string(N_max_core + 3) + "spdf");
  const auto cis2_basis = CI::basis_subset(ci_sp_basis, cis2_basis_string);

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
  // Often, we don't need to calculate new integrals.
  // It takes time to check if we need to, so faster to skip if we already
  // know all integrals exist.
  const auto no_new_integralsQ = input.get("no_new_integrals", false);
  {
    std::cout << "Calculate two-body Coulomb integrals: Q^k_abcd\n";
    std::cout << std::flush;

    const auto qk_filename =
      input.get("qk_file", wf.identity() + br_string + ".qk.abf");

    // Try to read from disk (may already have calculated Qk)
    qk.read(qk_filename);
    const auto existing = qk.count();

    if (!no_new_integralsQ) {
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
        [eF = wf.FermiLevel()](int, const DiracSpinor &s, const DiracSpinor &t,
                               const DiracSpinor &u, const DiracSpinor &v) {
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
      std::cout << std::flush;

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
  // nb: if Brueckner, Sigma1 already accounted for!
  std::cout << "Calculate one-body integrals.\n";
  std::cout << std::flush;
  const auto h1 =
    Brueckner ?
      CI::calculate_h1_table(ci_sp_basis, {}, {}, {}, false) :
    wf.Sigma() ?
      CI::calculate_h1_table(ci_sp_basis, *wf.Sigma(), include_Sigma1) :
      CI::calculate_h1_table(ci_sp_basis, core_s1, excited_s1, qk,
                             include_Sigma1);

  //----------------------------------------------------------------------------
  // Breit and QED

  if (wf.vHF()->Vrad()) {
    std::cout << "Including QED via HF\n";
  }
  if (wf.vHF()->vBreit()) {
    std::cout << "Including one-body Breit via HF\n";
  }
  std::cout << std::flush;

  // Creat Breit table (only for CI, not MBPT part)
  Coulomb::WkTable Bk;
  const auto Breit2 = input.get("Breit2", true);
  if (wf.vHF()->vBreit() && Breit2) {

    // use a subset of basis for Breit?
    const auto Breit_basis_string =
      input.get("Breit_basis", std::to_string(N_max_core + 6) + "spdf");
    const auto Breit_basis = CI::basis_subset(ci_sp_basis, Breit_basis_string);

    std::cout
      << "\nCalculate + include two-body Breit integrals for CI: B^k_abcd\n";
    std::cout << "For: " << DiracSpinor::state_config(Breit_basis) << "\n";
    std::cout << std::flush;

    const auto bk_filename =
      input.get("bk_file", wf.identity() + br_string + ".bk");

    Bk = CI::calculate_Bk(bk_filename, wf.vHF()->vBreit(), Breit_basis,
                          max_k_Coulomb, no_new_integralsQ);
  }

  //----------------------------------------------------------------------------
  // Calculate MBPT corrections to two-body Coulomb integrals

  Coulomb::LkTable Sk;
  if (include_Sigma2) {

    // Here, write basis info into filename, since these are _internal_ lines!
    const auto Sk_filename =
      input.get("sk_file", wf.identity() + "_" + std::to_string(n_min_core) +
                             "_" + DiracSpinor::state_config(excited_s2) +
                             (max_k_Coulomb >= 0 && max_k_Coulomb < 50 ?
                                "_" + std::to_string(max_k_Coulomb) :
                                "") +
                             br_string + ".sk.abf");

    std::cout << "\nCalculate two-body MBPT integrals: Σ^k_abcd\n";

    std::cout << "For: " << DiracSpinor::state_config(cis2_basis) << ", using "
              << DiracSpinor::state_config(excited_s2) << "\n";
    std::cout << std::flush;

    Sk = CI::calculate_Sk(Sk_filename, cis2_basis, core_s2, excited_s2, qk,
                          max_k_Coulomb, exclude_wrong_parity_box,
                          no_new_integralsQ);
  }

  //----------------------------------------------------------------------------
  const auto J_list = input.get("J", std::vector<int>{});
  const auto J_even_list = input.get("J+", J_list);
  const auto J_odd_list = input.get("J-", J_list);
  const auto num_solutions = input.get("num_solutions", 5);
  const auto all_below = input.get<double>("all_below");
  const auto ci_input = input.get("ci_input", std::string{""});
  const auto sort_output = input.get("sort_output", false);
  const auto parallel_ci = input.get("parallel_ci", true);

  const auto n_Js = J_even_list.size() + J_odd_list.size();

  levels.resize(n_Js);
  std::vector<std::ostringstream> os(n_Js); // for parallel output

  fmt::print("Running CI for {} J/pi's {}\n", n_Js,
             parallel_ci ? "in parallel" : "");
  std::cout << std::flush;

  {
    IO::ChronoTimer t2("CI Eigenvalues");
#pragma omp parallel for if (parallel_ci)
    for (std::size_t i = 0; i < n_Js; ++i) {
      const auto [twoj, pi] =
        i < J_even_list.size() ?
          std::pair{2 * J_even_list.at(i), +1} :
          std::pair{2 * J_odd_list.at(i - J_even_list.size()), -1};

      auto &output_stream = parallel_ci ? os.at(i) : std::cout;
      levels.at(i) = run_CI(ci_sp_basis, twoj, pi, num_solutions, all_below, h1,
                            qk, Bk, Sk, include_Sigma2, output_stream);
    }

    // If doing in parallel, output detailed output at end
    if (parallel_ci) {
      for (const auto &out : os)
        std::cout << out.str();
      std::cout << std::flush;
    }
  }

  //----------------------------------------------------------------------------

  // Find minimum (ground-state) energy (for level comparison)
  // Assumes levels for each J/Pi are sorted (they always are)
  double e0 = 0.0;
  for (const auto &lvl : levels) {
    if (lvl.num_solutions() == 0)
      continue;
    const auto e = lvl.energy(0);
    if (e < e0) {
      e0 = e;
    }
  }

  // This is just for screen output:
  // Sort output in pair {energy, output_string}, so we can optionally sort
  std::vector<std::pair<double, std::string>> E_output;
  for (const auto &Psi_Jpi : levels) {

    for (std::size_t i = 0; i < Psi_Jpi.num_solutions(); ++i) {

      const auto [config, pc, gJ, L, twoS] = Psi_Jpi.info(i);
      const auto iL = (int)std::round(L);
      const auto itwoS = (int)std::round(twoS);

      auto out_string = fmt::format(
        "{:<2} {:+2} {:>2}  {:<6s} {:2.0f}  {:<3s}  {:+12.8f}  {:+12.2f} "
        "{:12.2f}",
        Psi_Jpi.twoJ() / 2, Psi_Jpi.parity(), i, config, pc * 100.0,
        Term_Symbol(iL, itwoS, Psi_Jpi.parity()), Psi_Jpi.energy(i),
        Psi_Jpi.energy(i) * PhysConst::Hartree_invcm,
        (Psi_Jpi.energy(i) - e0) * PhysConst::Hartree_invcm);
      if (gJ != 0.0) {
        out_string += fmt::format("  {:.4f}", gJ);
      }
      E_output.emplace_back(Psi_Jpi.energy(i), out_string);
    }
    std::cout << std::flush;
  }

  // optionally, sort by energy
  if (sort_output) {
    std::sort(E_output.begin(), E_output.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
  }

  std::cout << "\nLevel Summary:\n\n";
  if (std::abs(wf.dalpha2()) > 1.0e-5) {
    std::cout << "d(α^2) = " << wf.dalpha2() << "\n";
  }
  std::cout
    << "J   π  #  conf.  %   Term   Energy(au)   Energy(/cm)   Level(/cm) "
       " gJ\n";
  for (const auto &[E, output] : E_output) {
    std::cout << output << "\n";
  }
  std::cout << std::flush;

  return levels;
}

//==============================================================================
//==============================================================================
//==============================================================================
PsiJPi run_CI(const std::vector<DiracSpinor> &ci_sp_basis, int twoJ, int parity,
              int num_solutions, std::optional<double> all_below,
              const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
              const Coulomb::WkTable &Bk, const Coulomb::LkTable &Sk,
              bool include_Sigma2, std::ostream &outstream) {

  auto printJ = [](int twoj) {
    return twoj % 2 == 0 ? std::to_string(twoj / 2) :
                           std::to_string(twoj) + "/2";
  };
  auto printPi = [](int pi) { return pi > 0 ? "even" : "odd"; };

  fmt::print(outstream, "Run CI for J={}, {} parity\n", printJ(twoJ),
             printPi(parity));
  outstream << std::flush;

  PsiJPi psi{twoJ, parity, ci_sp_basis};

  if (twoJ < 0) {
    outstream << "Fail: twoJ must >=0\n";
    return psi;
  }
  if (twoJ % 2 != 0) {
    outstream << "Fail: twoJ must be even for two-electron CSF\n";
    return psi;
  }

  const auto N_CSFs = psi.CSFs().size();
  outstream << "Total CSFs: " << psi.CSFs().size() << "\n";
  outstream << std::flush;

  //----------------------------------------------------------------------------

  // Construct the CI matrix:
  const auto br_ptr = !Bk.emptyQ() ? &Bk : nullptr;
  const auto s2_ptr = include_Sigma2 ? &Sk : nullptr;
  const auto Hci = CI::construct_Hci(psi, h1, qk, br_ptr, s2_ptr);

  //----------------------------------------------------------------------------

  if (all_below) {
    fmt::print(outstream, "Finding all solutions below {} cm^-1\n", *all_below);
  } else if (num_solutions > 0) {
    fmt::print(outstream, "Find first {} solutions\n", num_solutions);
  } else {
    fmt::print(outstream, "Finding all solutions\n", num_solutions);
  }

  {
    IO::ChronoTimer t2("");
    psi.solve(Hci, num_solutions, all_below);
    outstream << psi.num_solutions() << " eigenvalues: T = " << t2.reading_str()
              << "\n\n";
  }
  const auto E0 = psi.num_solutions() > 0 ? psi.energy(0) : 0.0;

  // For calculating g-factors
  DiracOperator::M1 m1{ci_sp_basis.front().grid(), PhysConst::alpha, 0.0};
  // only actually need to do this once..
  const auto m1_tab = ExternalField::me_table(ci_sp_basis, &m1);

  // Print details of each solution, unless we find all:
  const auto print_details = all_below || num_solutions > 0;
  const double minimum_percentage = 5.0; // min % to print

  // XXX nb: sometimes get's non-rel config wrong! (not a big issue)
  for (std::size_t i = 0; i < N_CSFs && i < psi.num_solutions(); ++i) {

    const auto pi = parity == 1 ? '+' : '-';
    if (print_details)
      fmt::print(outstream,
                 "{} {} {:<2}  {:+11.8f} au  {:+11.2f} cm^-1  {:11.2f} cm^-1\n",
                 twoJ / 2, pi, i, psi.energy(i),
                 psi.energy(i) * PhysConst::Hartree_invcm,
                 (psi.energy(i) - E0) * PhysConst::Hartree_invcm);

    // l's of the leading configuration (for gJ)
    int l1{-1}, l2{-1};
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
      if (cj > minimum_percentage && print_details) {
        fmt::print(outstream, "   {:<6s} {:5.3f}%\n", psi.CSF(j).config(true),
                   cj);
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
    const double gJ = twoJ != 0 ? tjs * m1AA_R / (0.5 * twoJ) : 0.0;

    // Determine Term Symbol, from g-factor
    const auto [S, L] = CI::Term_S_L(l1, l2, twoJ, gJ);

    if (print_details) {
      outstream << "   --------------\n";
      if (twoJ != 0) {
        outstream << "   gJ = " << gJ << "\n";
      }
    }

    // maximum relativistic config, into non-relativistic notation:
    const auto config = psi.CSF(max_j).config(false);

    // Percentage of that non-relativistic config*
    double pc = 0.0;
    for (std::size_t j = 0; j < N_CSFs; ++j) {
      if (psi.CSF(j).config(false) == config) {
        pc += psi.coef(i, j) * psi.coef(i, j);
      }
    }
    // * technically, might not be maximum non-rel config.. realistically, fine

    if (print_details) {
      fmt::print(outstream, "   {:<6s} {}\n", config,
                 Term_Symbol(twoJ, L, 2 * S, parity));
      outstream << "\n";
    }

    psi.update_config_info(i, {config, pc, gJ, 1.0 * L, 2.0 * S});
  }

  return psi;
}

} // namespace CI