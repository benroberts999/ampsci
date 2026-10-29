#include "ConfigurationInteraction.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Angular/include.hpp"
#include "CI_Integrals.hpp"
#include "CSF.hpp"
#include "ConfigurationSpace.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/Matrix.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/Sigma2.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

namespace CI {
Solutions configuration_interaction(const IO::InputBlock &input,
                                    const Wavefunction &wf) {

  // Check input options:
  input.check({
    {"ci_basis",
     "Basis used for CI expansion; must be a sub-set of full ampsci basis. "
     "With three valence electrons the number of CSFs grows quickly: a "
     "smaller basis (e.g., 12spdf) is usually required [default: 20spdf]"},
    {"J", "List of total angular momentum J for CI solutions (comma "
          "separated). Integers for two electrons; half-integers, as "
          "fractions (1/2, 3/2) or decimals (0.5, 1.5), for one or three. []"},
    {"J+", "As above, but for EVEN CSFs only (takes precedence over J)."},
    {"J-", "As above, but for ODD CSFs (takes precedence over J)."},
    {"num_solutions", "Number of CI solutions to find (for each J/pi) [5]"},
    {"all_below_cm",
     "Find all CI solutions for energies below this threshold, "
     "in inverse cm. Note that this is the total energy, not the "
     "excitation energy. If set, num_solutions is ignored."},
    {"sigma1", "Include one-body MBPT correlations? [false]"},
    {"sigma2", "Include two-body MBPT correlations? [false]"},
    {"Brueckner", "Use Brueckner (spectrum) states for CI basis? Must have "
                  "Correlations, Spectrum, and sigma1. [false]"},
    {"cis2_basis",
     "Pair of basis strings, e.g. '8spdf, 20spdf'. Sigma_2 is calculated for "
     "<vw|Q|xy> if (vw or xy) is in the first, and (vw and xy) are in second. "
     "All diagonal entries (xy = vw) up to second are included. "
     "Default: '(n+2)spdf..., 20spdf...' (each L in ci_basis)."},
    {"fk",
     "Effective QPQ ~ fk Q screening factors for the Sigma_2 Coulomb lines. "
     "Give an explicit list to use those values; set 'false' for no "
     "screening (or 'fk = 1;'). "
     "If blank (default), takes from correlation potential if it has them."},
    {"Feynman_screen",
     "Use Feynman method for Sigma_2 screening diagrams. Requires "
     "Correlations{Feynman = true;} [false]"},
    {"denominators",
     "Energy denominators in Sigma_2 [DFK]:\n"
     " - DFK: (Dzuba Flambaum Kozlov) target legs as Fermi, others as RS\n"
     " - BW: Brillouin-Wigner E0 (iterated; see BW_each_level).\n"
     " - RS: actual single-particle energies\n"
     " - Fermi: lowest excited energy of each kappa\n"
     " - Fermi0: lowest excited energy of all kappas\n"},
    {"sigma1_dSdE", "Sigma_1 at the level energy E0, via dSigma/dE [false]\n"
                    " - E0: lowest level of J/pi, or each (BW_each_level)\n"
                    " - faster with Correlations{ derivative = true; }\n"
                    " - For testing; usually should be left off"},
    {"BW_each_level", "E0 for each level, not just lowest of J/pi [false]\n"
                      "- for BW and sigma1_dSdE\n"
                      "- levels not orthogonal; one more solve per level"},
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
     "(or very large) to include all k. [10]"},
    //////
    {"num_electrons", "Number of valence electrons: 1, 2 or 3 [2]"},
    {"leading",
     "Leading configurations, from which the CI configurations are formed by "
     "excitations;\n"
     " - shells separated by '.', e.g., 'leading = 6s2.6p, 6s.6p2;'.\n"
     " - a hole (see 'holes') as a negative number, e.g., '6s-1.6p2'.\n"
     " - Default: the lowest few configurations of each parity\n"
     " - NB: Always check for truncated CI\n"
     " - Not required for one or two valence electrons (full CI)"},
    {"excitations",
     "Excitations from the leading configurations, formed stepwise: 'k' or "
     "'k:basis', e.g., '1, 2, 3:6spdf'. Step k moves one electron into its "
     "basis (default: ci_basis); the largest k is the number of steps, so "
     "'3:6spdf' is the same as '1, 2, 3:6spdf'; '0' gives the leading "
     "configurations only. Two full steps include every two-electron "
     "configuration. [1, 2]"},
    {"davidson", "Use iterative Davidson method for CI solutions "
                 "[default: true for 3 electrons, false otherwise]"},
    {"holes",
     "Core shells that may hold a hole, e.g., 'holes = 6s;' (both j). Adds "
     "CSFs with two electrons and one hole, e.g., 6p^2 6s^-1, to the CI "
     "space (one valence electron only). [none]"},
    {"ci_hole_basis", "Basis for the two electrons of the hole CSFs; subset "
                      "of ci_basis. [default: ci_basis]"},
    /////
    // {"twoJ", "List of 2J (integers) for CI solutions; takes precedence over "
    //          "J. Works for either number of electrons."},
    // {"twoJ+", "As above, but for EVEN CSFs only (takes precedence over J+)."},
    // {"twoJ-", "As above, but for ODD CSFs (takes precedence over J-)."},
    {"qk_file",
     "Filename for storing two-body Coulomb integrals. By default, is "
     "At.qk.abf, where At is atomic symbol + 'identity'. Set to 'false' to "
     "disable read/write."},
    {"sk_file",
     "Filename for storing two-body Sigma_2 integrals. By default, is "
     "At_b_hash.sk.abf, where At is atomic symbol + 'identity', b is the s2 "
     "(internal) basis, and hash encodes other settings that changes the "
     "integrals. Set to 'false' to disable read/write."},
    {"dsk_file",
     "As sk_file, but for the second Sigma_2 table used when denominators = "
     "BW. Default is At_b_hash.dsk.abf. nb: it holds "
     "S^k at the shifted E0; the derivative is formed from it after reading."},
    {"ds1_file",
     "Filename for storing the Sigma_1 derivative-correction tables "
     "(sigma1_dSdE). Default: At_b_hash.ds1.abf, where b is the s1 "
     "(internal) basis, and hash encodes other settings that change the "
     "tables. Not used if the Correlations block provides dSigma/dE (option "
     "derivative). Set to 'false' to disable read/write."},
    {"bk_file", "Filename for storing two-body Breit integrals. By default, is "
                "At.bk.abf, where At is atomic symbol + 'identity'. Set to "
                "'false' to disable read/write."},
    {"ci_file",
     "Filename for storing CI solutions. Default: At_b_hash.ci.abf, where At "
     "is atomic symbol + 'identity', b is ci_basis (prefixed with 1el_ or "
     "3el_ for one or three electrons), and hash encodes other settings that "
     "change the CI Hamiltonian. "
     "Set to 'false' to disable read/write."},
    {"read_only",
     "If true, will *only* read in the existing CI solutions, and will not "
     "calculate anything new, even if new states are requested. You must set "
     "the ci_basis (not read in) and the J/pi solutions you want; only these "
     "will be read in. Use if you know CI has already been completed. [false]"},
    {"no_new_integrals",
     "Usually false. If set to true, ampsci will not calculate any new "
     "Coulomb or Sigma_2 integrals, even if they are implied by the above "
     "settings. This saves time when we know all required integrals already "
     "exist, since the code doesn't need to check. [false]"},
    {"sort_output", "Sort output by energy? Default is to sort by J and Pi "
                    "first. [false]"},
    {"print_details", "Condition to print details of each CI solution "
                      "(otherwise just prints summary) [true]"},
    {"all_order_screen",
     "With Feynman_screen: screening to all orders; else second order only "
     "[true]"},
    {"hp_screen",
     "With Feynman_screen: include the hole-particle interaction [true]"},
    {"extrapolate_sigma2",
     "Extrapolate Sigma_2 to diagrams outside cis2_basis, using average "
     "correction ratios: S^k ~ h_k*Q^k, where h_k = <S^k/Q^k> is averaged over "
     "the calculated Sigma_2 integrals. Note: These are stored in the Sk "
     "table, but NOT written to disk [true]"},
    {"exclude_wrong_parity_box",
     "Excludes the Sigma_2 box corrections that "
     "have 'wrong' parity when calculating Sigma2 matrix elements (i.e., that "
     "don't match the Q^k selection rules). Note: If "
     "existing sk file already has these, they will be included [false]"},
  });

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return {};
  }
  //----------------------------------------------------------------------------

  // construct first, for RVO
  std::vector<PsiJPi> levels;

  //----------------------------------------------------------------------------
  // Single-particle basis:
  std::cout << "\nConstruct single-particle basis:\n";

  // Determine the sub-set of basis to use in CI:
  const auto basis_string = input.get("ci_basis", std::string{"20spdf"});

  // Number of valence electrons (in the CI space)
  const auto num_electrons = input.get("num_electrons", 2);
  if (num_electrons < 1 || num_electrons > 3) {
    fmt2::warning();
    fmt::print(": CI is implemented for one, two, or three valence electrons, "
               "not {}. Nothing to do\n",
               num_electrons);
    return {};
  }
  fmt::print("CI for {} valence electron{}\n", num_electrons,
             num_electrons == 1 ? "" : "s");
  // options to include MBPT
  const auto include_Sigma1 = input.get("sigma1", false);
  // Sigma_2 acts between pairs of valence electrons: none for one electron
  const auto include_Sigma2 = input.get("sigma2", false) && num_electrons > 1;
  if (input.get("sigma2", false) && !include_Sigma2) {
    std::cout << "nb: one valence electron: no Sigma_2 (no valence pairs)\n";
  }

  // Use Breuckner states for MBPT
  // nb: currently also use these Bruckner states in place of the "core"
  // I think this is the best option, due to orthogonality
  // These are not eigenstates of the V^HF potential used in the core
  // However, we don't explicitely use this for <v|h1|w>
  // We assume eignestates, so <v|h1|w> = E_w \delta_vw
  // (When not using Brueckner there is also <v|Sigma1|w>, which is not diag)
  // Perhaps a slight inconsistancy when including RPA into matrix elements?
  const auto Brueckner_raw = input.get("Brueckner", false);
  // nb: require wf.Sigma(): without it, the Spectrum is formed without
  // correlations, so the "Brueckner" states would just be regular basis states
  // (and Sigma1 would then be missed entirely)
  const auto Brueckner =
    Brueckner_raw && include_Sigma1 && wf.Sigma() && !wf.spectrum().empty();

  const auto &t_basis = Brueckner ? wf.spectrum() : wf.basis();

  // maximum n present in core: used for default basis
  const auto N_max_core = DiracSpinor::max_n(wf.core());

  // Select from basis those which match input 'basis_string'
  // exclude those in coreConfiguration
  // Negative-energy states must always excluded from CI
  const std::vector<DiracSpinor> ci_sp_basis = qip::select_if(
    CI::basis_subset(t_basis, basis_string, wf.coreConfiguration()),
    [](const auto &Fn) { return !Fn.negativeEnergyStateQ(); });

  // Print info re: basis to screen:
  std::cout << "\nUsing " << DiracSpinor::state_config(ci_sp_basis) << " = "
            << ci_sp_basis.size() << " orbitals in CI expansion\n";

  // Hole CSFs (two electrons + one hole): the core shells that may hold the
  // hole, as orbitals of the basis, and the basis of the two electrons
  std::vector<DiracSpinor> hole_orbitals;
  std::vector<DiracSpinor> ci_hole_basis;
  {
    const auto core_shells = AtomData::core_parser(wf.coreConfiguration());
    for (const auto &shell_string :
         input.get("holes", std::vector<std::string>{})) {
      const auto configuration = CI::parse_configuration(shell_string);
      if (configuration.size() != 1) {
        fmt2::warning();
        fmt::print(": hole shell '{}' is not a single shell: skipping\n",
                   shell_string);
        continue;
      }
      const auto &shell = configuration.front();
      const auto in_core = std::any_of(
        core_shells.cbegin(), core_shells.cend(),
        [&shell](const auto &c) { return c.n == shell.n && c.l == shell.l; });
      if (!in_core) {
        fmt2::warning();
        fmt::print(": hole shell '{}' is not a core shell: skipping\n",
                   shell_string);
        continue;
      }
      for (const auto &Fa : t_basis) {
        if (Fa.n() == shell.n && Fa.l() == shell.l) {
          hole_orbitals.push_back(Fa);
        }
      }
    }
    if (!hole_orbitals.empty() && num_electrons != 1) {
      fmt2::warning();
      std::cout << ": hole CSFs are available for one valence electron only: "
                   "not including them\n";
      hole_orbitals.clear();
    }
    if (!hole_orbitals.empty()) {
      ci_hole_basis =
        CI::basis_subset(ci_sp_basis, input.get("ci_hole_basis", basis_string));
      std::cout << "With hole CSFs: holes in "
                << DiracSpinor::state_config(hole_orbitals)
                << ", two electrons in "
                << DiracSpinor::state_config(ci_hole_basis) << "\n";
    }
  }

  // Configurations in the CI expansion: leading + excitations (including
  // those with a hole, if any)
  const ConfigurationSpace config_space{
    ci_sp_basis,
    num_electrons,
    input.get("leading", std::vector<std::string>{}),
    input.get("excitations", std::vector<std::string>{"1", "2"}),
    hole_orbitals,
    ci_hole_basis};
  std::cout << config_space.summary();

  if (Brueckner) {
    std::cout
      << "CI + MBPT + Brueckner method: using Brueckner states for CI basis\n";
  }
  if (Brueckner_raw && !Brueckner) {
    fmt2::warning();
    std::cout << ": Requested Brueckner method, but conditions not met:\n";
    if (!include_Sigma1) {
      std::cout << " - sigma1 option is false\n";
    }
    if (wf.spectrum().empty()) {
      std::cout << " - no Spectrum\n";
    }
    if (!wf.Sigma()) {
      std::cout << " - no correlation potential: Spectrum was formed without "
                   "Sigma\n";
    }
    std::cout << "Not using Brueckner method: using regular basis\n";
  }

  // aditional output string for br method:
  using namespace std::string_literals;
  const auto br_string = Brueckner ? "_bru"s : "";

  const auto read_only = input.get("read_only", false);
  if (read_only) {
    std::cout << "\nReading existing CI solutions and integrals - no new CI "
                 "performed\n";
  }

  //----------------------------------------------------------------------------

  // Determine different basis subsets

  // Details of MBPT
  const auto max_k_Coulomb = input.get("max_k", 10);
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
    DiracSpinor::split_by_energy(s1_basis, wf.FermiLevel(), n_min_core);
  const auto [core_s2, excited_s2] =
    DiracSpinor::split_by_energy(s2_basis, wf.FermiLevel(), n_min_core);

  // S2 corrections are included for <vw|Q|xy> with the pair vw, or xy, in
  // cis2_basis (the first subset), or diagonal (xy in the same nl shells as
  // vw), and all four orbitals in cis2_legs_basis (the second; all of the CI
  // basis if only the first is given). Defaults: n up to n_core + 2 for the
  // pair, n up to 20 for the legs, both with every l of the CI basis
  std::string ci_l_symbols;
  for (int l = 0; l <= DiracSpinor::max_l(ci_sp_basis); ++l) {
    ci_l_symbols += AtomData::l_symbol(l);
  }
  const auto cis2_basis_strings = input.get(
    "cis2_basis",
    std::vector<std::string>{std::to_string(N_max_core + 2) + ci_l_symbols,
                             "20" + ci_l_symbols});
  const auto cis2_basis =
    CI::basis_subset(ci_sp_basis, cis2_basis_strings.front());
  const auto cis2_legs_basis =
    cis2_basis_strings.size() > 1 ?
      CI::basis_subset(ci_sp_basis, cis2_basis_strings.at(1)) :
      ci_sp_basis;

  // Screening factors for Sigma_2 (QPQ ~ fk*Q):
  // - explicit list: used as given
  // - fk = false; : no screening (fk = 1)
  // - blank: if the correlation potential stores fk (Feynman screening),
  //   use the average of the lowest s, p, and d factors
  const auto fk_string = input.get("fk", std::string{});
  const auto fk_false = qip::ci_compare(fk_string, "false");
  auto fk =
    fk_false ? std::vector<double>{} : input.get("fk", std::vector<double>{});
  bool fk_from_Sigma = false;
  if (fk.empty() && !fk_false && include_Sigma2 && wf.Sigma()) {
    fk = wf.Sigma()->average_fk(2);
    fk_from_Sigma = !fk.empty();
  }

  // Sigma_2 loop diagrams (a1, b1) by the Feynman method, static limit,
  // with the correlation potential's sub-grid and Feynman parameters. The
  // energy dependence of their denominators is kept at second order.
  const auto Feynman_screen_requested =
    include_Sigma2 && input.get("Feynman_screen", false);
  const auto Feynman_screen =
    Feynman_screen_requested && wf.Sigma() &&
    wf.Sigma()->method() == MBPT::SigmaMethod::Feynman;
  const auto all_order_screen = input.get("all_order_screen", true);
  const auto hp_screen = input.get("hp_screen", true);
  // Tag for the sk and ci file hashes
  const auto Feynman_screen_tag =
    Feynman_screen ? std::string("Fscr") + (all_order_screen ? "ao" : "2") +
                       (hp_screen ? "hp" : "") + ";" :
                     "";
  if (Feynman_screen_requested && !Feynman_screen) {
    fmt2::warning();
    std::cout << ": Requested Feynman screening in Sigma_2, but the "
                 "correlation potential does not use the Feynman method."
                 " Not including it\n";
  }

  // Extrapolate Sigma_2 (via average screening) for diagrams outside
  // cis2_basis
  const auto extrapolate_sigma2 = input.get("extrapolate_sigma2", true);

  // Iterative (dSigma/dE) energy correction for Sigma_1
  const auto sigma1_dSdE = input.get("sigma1_dSdE", false);
  const auto BW_each_level = input.get("BW_each_level", false);

  const auto denominators =
    MBPT::parse_Denominators(input.get("denominators", "DFK"s));

  // Shift Sigma_2 to each J/Pi's E0 (BW only; see denominators doc). E0 is
  // found iteratively, with or without the dSigma/dE correction to Sigma_1
  const auto iterative_correction_sigma2 =
    include_Sigma2 && denominators == MBPT::Denominators::BW;

  // Internal reference E0 that the Sigma_2 integrals are tabulated at. The
  // shift to each J/Pi is resummed, so results barely depend on it; it need
  // only be roughly right (zeroth-order ground-state pair energy)
  const auto E0_sigma2 =
    iterative_correction_sigma2 ? 2.0 * DiracSpinor::min_En(ci_sp_basis) : 0.0;

  //----------------------------------------------------------------------------

  if (include_MBPT && !read_only) {
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

      std::cout << "Including <vw|Σ_2|xy> for: vw/xy in "
                << DiracSpinor::state_config(cis2_basis) << ", vw+xy in "
                << DiracSpinor::state_config(cis2_legs_basis) << " (+diag)\n";
      if (extrapolate_sigma2) {
        std::cout << "(and extrapolating beyond using hk=<Sk/Qk>)\n";
      }
      if (exclude_wrong_parity_box) {
        std::cout
          << "Excluding the Σ_2 diagrams that have the 'wrong' parity\n";
      }
      std::cout << "Using: " << MBPT::parse_Denominators(denominators)
                << " denominators\n";
    }
    std::cout << "\n";
  }

  //----------------------------------------------------------------------------

  // The correlation potential holds dSigma/dE matrices (Correlations option
  // `derivative`): build the Sigma_1 correction tables directly from it -
  // consistent with the actual Sigma (any method), and no Sigma_1 qk
  // integrals or Goldstone tabulation needed
  const bool CP_derivative =
    sigma1_dSdE && include_Sigma1 && wf.Sigma() && wf.Sigma()->has_derivative();

  // Lookup table; stores all qk's
  Coulomb::QkTable qk;
  // Often, we don't need to calculate new integrals.
  // It takes time to check if we need to, so faster to skip if we already
  // know all integrals exist.
  // If read_only, we still read the integrals in (the modules need them), but
  // calculate nothing new
  const auto no_new_integralsQ =
    input.get("no_new_integrals", false) || read_only;

  // With no_new_integrals (or read_only), an empty table is not an error
  // (may be intended), but is usually a mistake (e.g., missing file): warn
  const auto warn_if_empty = [no_new_integralsQ](std::size_t count,
                                                 const std::string &name,
                                                 const std::string &filename) {
    if (no_new_integralsQ && count == 0) {
      fmt2::warning();
      fmt::print(": no_new_integrals (or read_only) is set, but no {} "
                 "integrals were read (from: {}).\n"
                 "They will NOT be calculated: zero {} integrals used!\n",
                 name, filename, name);
      std::cout << std::flush;
    }
  };

  {
    std::cout << "\n"
              << (no_new_integralsQ ? "Read" : "Calculate")
              << " two-body Coulomb integrals: Q^k_abcd\n";
    std::cout << std::flush;

    const auto qk_filename =
      input.get("qk_file", wf.identity() + br_string + ".qk.abf");

    // Try to read from disk (may already have calculated Qk)
    qk.read(qk_filename);
    const auto existing = qk.count();
    warn_if_empty(existing, "Coulomb Q^k", qk_filename);

    if (!no_new_integralsQ) {
      // Try to limit number of Coulomb integrals we calculate
      // use whole basis (these are used inside Sigma_2)
      // If not including MBPT, only need to caculate smaller set of integrals

      // First, calculate the integrals between ci basis states (and the
      // hole orbitals, if any):
      {
        const auto ci_and_holes = qip::merge(ci_sp_basis, hole_orbitals);
        std::cout << "For: " << DiracSpinor::state_config(ci_and_holes)
                  << " (for CI)\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(ci_and_holes);
        qk.fill(ci_and_holes, yk, max_k_Coulomb, false);
      }

      // Selection function for which Qk's to calculate.
      // For Sigma, we only need those with 1 or 2 core electrons
      // i.e., no Q_vwxy, Q_vabc, or Q_abcd
      // Note: we *do* need Q_vwxy for the CI part (but with smaller basis)
      const auto select_Q_sigma =
        [eF = wf.FermiLevel()](int, const DiracSpinor &s, const DiracSpinor &t,
                               const DiracSpinor &u, const DiracSpinor &v) {
          // Only calculate Coulomb integrals with 1 or 2 electrons in the core
          auto num = Coulomb::number_below_Fermi(s, t, u, v, eF);
          return num == 1 || num == 2;
        };

      // Then, add those required for Sigma_1 (unless we have matrix!)
      // (With matrix, still required if calculating derivative correction -
      // except when the matrix supplies dSigma/dE itself: CP_derivative)
      const bool need_s1_qk =
        include_Sigma1 && (!wf.Sigma() || (sigma1_dSdE && !CP_derivative));
      if (need_s1_qk) {
        const auto temp_basis = qip::merge(core_s1, excited_s1);
        std::cout << "and: " << DiracSpinor::state_config(temp_basis)
                  << " (for MBPT)\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(temp_basis);
        qk.fill_if(temp_basis, yk, select_Q_sigma, max_k_Coulomb, false);
      }

      // Then, add those required for Sigma_2 (unless we already did Sigma_1)
      if (include_Sigma2 && !need_s1_qk) {
        const auto temp_basis = qip::merge(core_s2, excited_s2);
        std::cout << "and: " << DiracSpinor::state_config(temp_basis)
                  << " (for MBPT)\n"
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
  // nb: not stored on disk, so is re-calculated even if read_only
  std::cout << "Calculate one-body integrals.\n";
  std::cout << std::flush;
  auto h1 = Brueckner ?
              CI::calculate_h1_table(ci_sp_basis, {}, {}, {}, false) :
            wf.Sigma() ?
              CI::calculate_h1_table(ci_sp_basis, *wf.Sigma(), include_Sigma1) :
              CI::calculate_h1_table(ci_sp_basis, core_s1, excited_s1, qk,
                                     include_Sigma1);
  // The hole orbitals: their orbital energies only (no Sigma_1)
  for (const auto &Fa : hole_orbitals) {
    h1.add(Fa, Fa, Fa.en());
  }
  if (!hole_orbitals.empty() && include_Sigma1) {
    std::cout << "nb: no Sigma_1 for the hole orbitals\n";
  }

  // Derivative (dSigma/dE) correction for Sigma_1
  CI::Sigma1Correction s1_corr;
  if (sigma1_dSdE && include_Sigma1) {
    std::cout << "\nIncluding derivative (dSigma/dE) correction for Sigma_1\n";
    std::cout << std::flush;

    if (CP_derivative) {
      // Tables directly from the correlation potential (fast; no file cache
      // needed): consistent with the actual Sigma, whatever the method
      std::cout << "Using Sigma and dSigma/dE from the correlation potential ("
                << wf.Sigma()->method_string() << ")\n";
      s1_corr = CI::calculate_dSdE_correction(ci_sp_basis, *wf.Sigma());
    } else {

      // Second-order (Goldstone) tables, from the qk integrals.
      // Filename encodes everything that changes the tables (cf. sk_file):
      // the ci basis (which pairs + e_sigma convention), the s1 internal
      // lines, and the Coulomb-integral settings
      const auto ds1_settings = DiracSpinor::state_config(ci_sp_basis) + ";" +
                                DiracSpinor::state_config(core_s1) + ";" +
                                std::to_string(n_min_core) + ";" +
                                std::to_string(max_k_Coulomb) + ";";
      const auto ds1_basename = wf.identity() + br_string + "_" +
                                DiracSpinor::state_config(excited_s1) + "_" +
                                qip::hash_string(ds1_settings);
      const auto ds1_filename =
        input.get("ds1_file", ds1_basename + ".ds1.abf");

      const auto read_ds1 = ds1_filename != "false" &&
                            s1_corr.read_write(ds1_filename, IO::FRW::read);
      if (read_ds1) {
        std::cout << "Read Sigma_1 correction tables from " << ds1_filename
                  << "\n";
      } else {
        s1_corr =
          CI::calculate_dSdE_correction(ci_sp_basis, core_s1, excited_s1, qk);
        if (ds1_filename != "false") {
          std::cout << "Writing Sigma_1 correction tables to " << ds1_filename
                    << "\n";
          s1_corr.read_write(ds1_filename, IO::FRW::write);
        }
      }
    }
    std::cout << std::flush;
  }

  // E0 is found iteratively for the dSigma/dE correction, and for BW Sigma_2.
  // BW needs the orbital energies (initial E0; the pair energies, for three
  // electrons), which come with the dSigma/dE tables; without those, they are
  // added on their own
  if (iterative_correction_sigma2 && s1_corr.empty()) {
    for (const auto &v : ci_sp_basis) {
      s1_corr.en[v.nk_index()] = v.en();
    }
  }
  const bool find_E0 = !s1_corr.empty() || iterative_correction_sigma2;

  //----------------------------------------------------------------------------
  // Breit and QED

  if (wf.vHF()->Vrad()) {
    std::cout << "Including QED via HF\n";
  }
  if (wf.vHF()->vBreit()) {
    std::cout << "Including one-body Breit via HF\n";
  }
  std::cout << std::flush;

  // Creat Breit table (only for CI, not MBPT part).
  // Two-body Breit acts between pairs of valence electrons: none for one,
  // unless there are hole CSFs (two particles and a hole)
  Coulomb::WkTable Bk;
  const auto Breit2 =
    input.get("Breit2", true) && (num_electrons > 1 || !hole_orbitals.empty());
  if (wf.vHF()->vBreit() && Breit2) {

    // use a subset of basis for Breit? (plus the hole orbitals, if any)
    const auto Breit_basis_string =
      input.get("Breit_basis", std::to_string(N_max_core + 6) + "spdf");
    const auto Breit_basis = qip::merge(
      CI::basis_subset(ci_sp_basis, Breit_basis_string), hole_orbitals);

    std::cout << (no_new_integralsQ ? "\nRead" : "\nCalculate +")
              << " include two-body Breit integrals for CI: B^k_abcd\n";
    std::cout << "For: " << DiracSpinor::state_config(Breit_basis) << "\n";
    std::cout << std::flush;

    const auto bk_filename =
      input.get("bk_file", wf.identity() + br_string + ".bk.abf");

    Bk = CI::calculate_Bk(bk_filename, wf.vHF()->vBreit(), Breit_basis,
                          max_k_Coulomb, no_new_integralsQ);
    warn_if_empty(Bk.count(), "Breit B^k", bk_filename);
  }

  //----------------------------------------------------------------------------
  // Calculate MBPT corrections to two-body Coulomb integrals

  // The default sk and ci filenames are At_basis_hash, where At identifies the
  // atom, basis is the relevant basis (kept readable), and hash is a short tag
  // for every _other_ setting that changes the stored values. Any change to
  // those settings gives a new filename, so nothing is silently re-used.
  // Settings common to both files:
  std::string common_settings = std::to_string(n_min_core) + ";" +
                                std::to_string(max_k_Coulomb) + ";" +
                                (Brueckner ? "bru" : "") + ";";
  for (const auto f : fk) {
    common_settings += std::to_string(f) + ",";
  }
  common_settings += ";";

  Coulomb::WkTable Sk;
  // dS^k/dE0; empty unless iterative_correction_sigma2
  Coulomb::WkTable dSk;
  // Average S^k/Q^k ratios, for extrapolating Sigma_2 beyond cis2_basis
  std::vector<double> hk;
  if (include_Sigma2) {

    // Here, write basis info into filename, since these are _internal_ lines!
    const auto sk_settings =
      common_settings + MBPT::parse_Denominators(denominators) + ";" +
      (exclude_wrong_parity_box ? "xb" : "") + ";" +
      (iterative_correction_sigma2 ? std::to_string(E0_sigma2) + ";" : "") +
      Feynman_screen_tag;
    const auto sk_basename = wf.identity() + "_" +
                             DiracSpinor::state_config(excited_s2) + "_" +
                             qip::hash_string(sk_settings);
    const auto Sk_filename = input.get("sk_file", sk_basename + ".sk.abf");

    std::cout << (no_new_integralsQ ? "\nRead" : "\nCalculate")
              << " two-body MBPT integrals: Σ^k_abcd\n";

    std::cout << "For: " << DiracSpinor::state_config(cis2_basis) << " (pair), "
              << DiracSpinor::state_config(cis2_legs_basis) << " (all), using "
              << DiracSpinor::state_config(excited_s2) << "\n";

    // output screening factors to check
    if (!fk.empty()) {
      if (fk_from_Sigma) {
        std::cout << "fk from correlation potential (average of lowest "
                     "s, p, d)\n";
      }
      std::cout << "Effective screening into Σ_2:\n fk = [";
      for (std::size_t i = 0; i < fk.size(); ++i) {
        fmt::print("{:.3f}{}", fk.at(i), i + 1 == fk.size() ? "]\n" : ", ");
      }
      if (Feynman_screen) {
        std::cout << " * fk only used for (non-screening) diagrams; Feynman "
                     "for others\n";
      }
    }

    std::cout << std::flush;

    if (iterative_correction_sigma2) {
      fmt::print(
        "Brillouin-Wigner Σ_2: at E0 = {:.2f}, + dΣ^k/dE0 (each J/pi)\n",
        E0_sigma2);
    }

    // Static screened-Coulomb kernels for the loop diagrams. The
    // correlation potential's own Feynman object is used if it has the
    // requested screening and hole-particle options; otherwise one is
    // constructed with its sub-grid and parameters (Q*Pi*Q from the disk
    // cache, if present)
    std::vector<MBPT::RMatrix> screening_kernels;
    if (Feynman_screen && !no_new_integralsQ) {
      fmt::print(
        "Σ_2 direct (screening) diagrams by Feynman method: {} screening{}\n",
        all_order_screen ? "all-orders" : "second-order",
        hp_screen ? " + hole-particle" : "");
      std::cout << std::flush;
      const auto *feynman = wf.Sigma()->feynman();
      if (feynman != nullptr && feynman->screening() == all_order_screen &&
          feynman->hole_particle() == hp_screen) {
        screening_kernels = feynman->dQ_0_k();
      } else {
        const auto new_feynman =
          wf.Sigma()->make_feynman(all_order_screen ? MBPT::Screening::include :
                                                      MBPT::Screening::exclude,
                                   hp_screen ? MBPT::HoleParticle::include :
                                               MBPT::HoleParticle::exclude);
        screening_kernels = new_feynman->dQ_0_k();
      }
    }
    const auto *screening_ptr =
      screening_kernels.empty() ? nullptr : &screening_kernels;

    Sk = MBPT::calculate_Sk(Sk_filename, cis2_legs_basis, cis2_basis, core_s2,
                            excited_s2, qk, max_k_Coulomb,
                            exclude_wrong_parity_box, denominators,
                            no_new_integralsQ, fk, E0_sigma2, screening_ptr);
    warn_if_empty(Sk.count(), "Sigma_2 S^k", Sk_filename);

    if (iterative_correction_sigma2) {
      // dS^k/dE0, by finite difference in E0. The denominators are exactly
      // linear in E0, so the only error is the (second-order) curvature of
      // the sum, which the resummation in CI::corrected_Sk accounts for.
      // nb: the file holds S^k(E0 + delta); it is turned into the derivative
      // here, after reading, so a re-read gives the same thing.
      constexpr auto delta_E0 = 0.01;
      const auto dSk_filename = input.get("dsk_file", sk_basename + ".dsk.abf");
      std::cout << "\nCalculate Σ^k_abcd at E0 + " << delta_E0
                << " au, for dΣ^k/dE0\n"
                << std::flush;
      dSk = MBPT::calculate_Sk(
        dSk_filename, cis2_legs_basis, cis2_basis, core_s2, excited_s2, qk,
        max_k_Coulomb, exclude_wrong_parity_box, denominators,
        no_new_integralsQ, fk, E0_sigma2 + delta_E0, screening_ptr);
      warn_if_empty(dSk.count(), "dSigma_2 (dsk)", dSk_filename);
      // Both tables have the same entries (same fill), so this is safe
      for (std::size_t k = 0; k < dSk->size(); ++k) {
        for (auto &[index, value] : dSk->at(k)) {
          value = (value - Sk.Q(int(k), index)) / delta_E0;
        }
      }
    }

    if (extrapolate_sigma2) {
      // Store the <hk> average ratios.
      // S^k = hk*Q^k formed as needed
      hk = MBPT::average_hk(Sk, qk, cis2_basis, max_k_Coulomb);
      std::cout << "Extrapolating Sigma_2 beyond cis2 basis, using average "
                   "correction ratios:\n ";
      std::cout << " hk = [";
      for (std::size_t k = 0; k < hk.size(); ++k) {
        fmt::print("{:.1e}{}", hk[k], k + 1 == hk.size() ? "]\n" : ", ");
      }
      std::cout << "\n";
    }
  }

  //----------------------------------------------------------------------------
  const auto num_solutions = input.get("num_solutions", 5);
  const auto all_below_cm = input.get<double>("all_below_cm");
  const auto solver = input.get("davidson", num_electrons >= 3) ?
                        Solver::Davidson :
                        Solver::dense;
  const auto sort_output = input.get("sort_output", false);
  const auto print_details = input.get("print_details", true);

  // {2J, parity} for each requested J/Pi
  const auto J_pi_list = requested_J_pi(input, num_electrons);
  levels.resize(J_pi_list.size());

  // CI solutions file: At_basis_hash.ci.abf; see above.
  // Everything that changes the CI Hamiltonian, apart from the ci_basis, goes
  // into the hash. nb: converged E0 is not included: it is derived
  // (self-consistently) from the other settings, so is reproducible.
  std::string ci_settings = common_settings;
  if (!config_space.unrestricted()) {
    ci_settings += "cs;" + config_space.summary();
  }
  if (include_Sigma1) {
    ci_settings += "s1;" + DiracSpinor::state_config(s1_basis) + ";";
  }
  // Sigma_1 method identity (method, screening, ladder): catches Correlation
  // Potential changes not otherwise encoded. Fitted lambdas (rounded) change
  // the Hamiltonian too, but are kept out of method_string
  if (include_Sigma1 && wf.Sigma()) {
    ci_settings += wf.Sigma()->method_string() + ";";
    const auto lambdas = wf.Sigma()->lambda_string();
    if (!lambdas.empty()) {
      ci_settings += lambdas + ";";
    }
  }
  if (!s1_corr.empty()) {
    // Tag the source: CP-provided (all-orders) tables differ from the
    // second-order Goldstone ones
    ci_settings += CP_derivative ? "dSdE-CP;" : "dSdE;";
  }
  if (find_E0 && BW_each_level) {
    ci_settings += "BWeach;";
  }
  if (include_Sigma2) {
    ci_settings +=
      "s2;" + MBPT::parse_Denominators(denominators) + ";" +
      (iterative_correction_sigma2 ? "bw" + std::to_string(E0_sigma2) + ";" :
                                     "") +
      DiracSpinor::state_config(s2_basis) + ";" +
      DiracSpinor::state_config(cis2_basis) + ";" +
      DiracSpinor::state_config(cis2_legs_basis) + ";" +
      (exclude_wrong_parity_box ? "xb" : "") + ";" +
      (extrapolate_sigma2 ? "ex" : "") + ";" + Feynman_screen_tag;
  }
  if (!Bk.emptyQ()) {
    ci_settings +=
      input.get("Breit_basis", std::to_string(N_max_core + 6) + "spdf") + ";";
  }
  if (!hole_orbitals.empty()) {
    ci_settings += "holes;" + DiracSpinor::state_config(hole_orbitals) + ";" +
                   DiracSpinor::state_config(ci_hole_basis) + ";";
  }

  // nb: two electrons is not tagged: two-electron filenames do not depend on
  // num_electrons
  const auto electrons_tag =
    num_electrons == 2 ? "" : std::to_string(num_electrons) + "el_";
  const auto default_ci_fname = wf.identity() + "_" + electrons_tag +
                                basis_string + "_" +
                                qip::hash_string(ci_settings) + ".ci.abf";
  const auto ci_fname_input = input.get("ci_file", default_ci_fname);
  const auto ci_fname =
    (ci_fname_input == "false") ? std::string{} : ci_fname_input;
  if (!ci_fname.empty()) {
    std::cout << "CI solutions file: " << ci_fname << "\n";
  }

  fmt::print("Running CI for {} J/pi's\n\n", J_pi_list.size());
  std::cout << std::flush;

  {
    IO::ChronoTimer t2("CI Eigenvalues");
    for (std::size_t i = 0; i < J_pi_list.size(); ++i) {
      const auto [twoj, pi] = J_pi_list.at(i);

      levels.at(i) =
        run_CI(ci_sp_basis, twoj, pi, num_electrons, num_solutions,
               all_below_cm, h1, qk, Bk, Sk, include_Sigma2, print_details,
               read_only, std::cout, ci_fname, find_E0 ? &s1_corr : nullptr, hk,
               dSk.emptyQ() ? nullptr : &dSk, E0_sigma2, solver, BW_each_level,
               config_space, hole_orbitals, ci_hole_basis);
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

  // E0 is found per (J, parity) inside run_CI, so the value stored here is
  // only a fallback: it is used when a Hci is later re-constructed (from the
  // stored Integrals) for a (J, parity) that was never solved - e.g. the
  // intermediate states of the sum-over-states. Use the ground state.
  if (find_E0 && e0 < 0.0) {
    s1_corr.E0 = e0;
  }

  // This is just for screen output:
  // Sort output in pair {energy, output_string}, so we can optionally sort
  std::vector<std::pair<double, std::string>> E_output;
  for (const auto &Psi_Jpi : levels) {

    for (std::size_t i = 0; i < Psi_Jpi.num_solutions(); ++i) {

      const auto &[config, pc, gJ, L, twoS, L2, S2] = Psi_Jpi.info(i);
      (void)L2, (void)S2;
      const auto iL = (int)std::round(L);
      const auto itwoS = (int)std::round(twoS);

      // Term symbol only where L, S were assigned
      const auto term =
        L < 0.0 ? std::string{} : Term_Symbol(iL, itwoS, Psi_Jpi.parity());
      auto out_string = fmt::format(
        "{:<3} {:+2} {:>2}  {:<9s} {:2.0f}  {:<3s}  {:+12.8f}  {:+12.2f} "
        "{:12.2f}",
        J_string(Psi_Jpi.twoJ()), Psi_Jpi.parity(), i, config, pc * 100.0, term,
        Psi_Jpi.energy(i), Psi_Jpi.energy(i) * PhysConst::Hartree_invcm,
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
    << "J    π  #  conf.     %   Term   Energy(au)   Energy(/cm)   Level(/cm) "
       " gJ\n";
  for (const auto &[E, output] : E_output) {
    std::cout << output << "\n";
  }
  std::cout << std::flush;

  // Keep the integrals: they are required to construct the CI Hamiltonian for
  // any other J/parity later on (e.g., for the mixed states)
  Solutions out;
  out.levels = std::move(levels);
  out.integrals.ci_basis = ci_sp_basis;
  out.integrals.num_electrons = num_electrons;
  out.integrals.config_space = config_space;
  out.integrals.holes = std::move(hole_orbitals);
  out.integrals.holes_basis = std::move(ci_hole_basis);
  out.integrals.brueckner = Brueckner;
  out.integrals.h1 = std::move(h1);
  out.integrals.qk = std::move(qk);
  out.integrals.Bk = std::move(Bk);
  out.integrals.Sk = std::move(Sk);
  out.integrals.dSk = std::move(dSk);
  out.integrals.E0_sigma2 = E0_sigma2;
  out.integrals.hk = std::move(hk);
  out.integrals.s1_corr = std::move(s1_corr);

  return out;
}

//==============================================================================
std::vector<std::pair<int, int>> requested_J_pi(const IO::InputBlock &input,
                                                int num_electrons) {

  // Requested J's, as 2J: 'twoJ' (integers; works for either number of
  // electrons) takes precedence over 'J' (integers or fractions over 2).
  // '+'/'-' versions take precedence over the plain ones, for their parity
  const auto twoJ_list = [&input](const std::string &pi_tag) {
    for (const auto &name : {"twoJ" + pi_tag, std::string{"twoJ"}}) {
      if (input.get<std::vector<int>>(name)) {
        return input.get(name, std::vector<int>{});
      }
    }
    std::vector<int> list;
    for (const auto &name : {"J" + pi_tag, std::string{"J"}}) {
      const auto J_strings = input.get<std::vector<std::string>>(name);
      if (!J_strings)
        continue;
      for (const auto &J_text : *J_strings) {
        const auto twoJ = parse_twoJ(J_text);
        if (twoJ) {
          list.push_back(*twoJ);
        } else {
          fmt2::warning();
          fmt::print(": Could not parse J = '{}' (in {}): skipping\n", J_text,
                     name);
        }
      }
      return list;
    }
    return list;
  };
  const auto twoJ_even_list = twoJ_list("+");
  const auto twoJ_odd_list = twoJ_list("-");

  // 2J must be even for two electrons, odd for three (half-integer J): others
  // cannot be formed
  std::vector<std::pair<int, int>> J_pi_list;
  // Both parity lists may come from the same J list: warn once per J
  std::set<int> impossible_twoJ;
  for (const auto [list, pi] :
       {std::pair{&twoJ_even_list, +1}, std::pair{&twoJ_odd_list, -1}}) {
    for (const auto twoJ : *list) {
      if (twoJ < 0 || (twoJ + num_electrons) % 2 != 0) {
        impossible_twoJ.insert(twoJ);
        continue;
      }
      J_pi_list.push_back({twoJ, pi});
    }
  }
  for (const auto twoJ : impossible_twoJ) {
    fmt2::warning();
    fmt::print(": J = {} is not possible for {} electrons: skipping\n",
               J_string(twoJ), num_electrons);
  }
  return J_pi_list;
}

//==============================================================================
//==============================================================================
//==============================================================================
PsiJPi run_CI(const std::vector<DiracSpinor> &ci_sp_basis, int twoJ, int parity,
              int num_electrons, int num_solutions,
              std::optional<double> all_below_cm,
              const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
              const Coulomb::WkTable &Bk, const Coulomb::WkTable &Sk,
              bool include_Sigma2, bool print_details, bool read_only,
              std::ostream &outstream, const std::string &ci_fname,
              const Sigma1Correction *s1c, const std::vector<double> &hk,
              const Coulomb::WkTable *dSk, double E0_sigma2, Solver solver,
              bool BW_each_level, const ConfigurationSpace &config_space,
              const std::vector<DiracSpinor> &holes,
              const std::vector<DiracSpinor> &holes_basis) {

  auto printPi = [](int pi) { return pi > 0 ? "even" : "odd"; };

  fmt::print(outstream, "CI: J={}, {} parity\n", J_string(twoJ),
             printPi(parity));
  outstream << std::flush;

  if (twoJ < 0 || (twoJ + num_electrons) % 2 != 0) {
    fmt::print(outstream, "No please: J={} is not possible for {} electrons\n",
               J_string(twoJ), num_electrons);
    return PsiJPi{};
  }

  PsiJPi psi{twoJ,         parity, ci_sp_basis, num_electrons,
             config_space, holes,  holes_basis};

  const auto N_CSFs = psi.num_CSFs();
  if (psi.has_holes()) {
    const auto num_hole_CSFs = std::size_t(
      std::count_if(psi.CSFs().cbegin(), psi.CSFs().cend(),
                    [](const auto &csf) { return csf.num_holes() > 0; }));
    outstream << "Total CSFs: " << N_CSFs << " (" << N_CSFs - num_hole_CSFs
              << " + " << num_hole_CSFs << " with a hole)\n"
              << std::flush;
  } else {
    outstream << "Total CSFs: " << N_CSFs << "\n" << std::flush;
  }
  if (N_CSFs == 0) {
    outstream << "No CSFs can be formed for this J/parity from the basis and "
                 "configurations\n\n";
    return psi;
  }
  // Davidson (sparse Hamiltonian) or dense solution
  const auto davidson = solver == Solver::Davidson;

  // The dense Hamiltonian is N^2 doubles, and diagonalising holds a copy of
  // it: warn for very large blocks
  const auto Hci_GB = (8.0 / 1.0e9) * double(N_CSFs) * double(N_CSFs);
  if (!davidson && Hci_GB > 6.0) {
    fmt::print(outstream,
               "nb: CI matrix (dense) requires ~{:.0f} GB; consider using "
               "Davidson\n",
               2.0 * Hci_GB);
  }

  // Try to read existing solutions from file
  bool read_ok =
    !ci_fname.empty() && psi.read_write(ci_fname, IO::FRW::read, outstream);
  const auto num_read_solutions = (int)psi.num_solutions();

  // Solve only if we don't have enough solutions already
  const auto n_required =
    (num_solutions <= 0 || all_below_cm) ? (int)N_CSFs : num_solutions;
  const bool need_solve = !read_ok || num_read_solutions < n_required;

  if (read_only && need_solve) {
    fmt::print(
      outstream,
      "Warning: Requested {} solutions, but only read {} from {}.\n"
      "         Running 'read only' CI, so these will not be computed.\n"
      "         Re-run with `read_only = false` to calculate these levels\n",
      n_required, num_read_solutions, ci_fname);
  }

  if (need_solve && !read_only) {
    // Construct the CI matrix:
    const auto br_ptr = !Bk.emptyQ() ? &Bk : nullptr;
    const auto s2_ptr = include_Sigma2 ? &Sk : nullptr;

    // The dSigma/dE correction and the BW Sigma_2 shift depend on E0, the
    // lowest energy of _this_ (J, parity), which we only know once we have
    // solved: find it self-consistently (see iterate_E0).
    // nb: the final solve re-uses the converged Hci, so no matrix is built
    // that we do not use.
    // Unconverged (Davidson) solutions are used, but not kept for re-use
    bool write_solutions = true;
    if (davidson) {
      using SparseMatrix = LinAlg::SparseSymmetricMatrix<double>;
      SparseMatrix Hci;
      if (s1c != nullptr) {
        Hci = CI::iterate_E0<SparseMatrix>(&psi, *s1c, h1, qk, br_ptr, s2_ptr,
                                           hk, dSk, E0_sigma2, outstream);
      } else {
        // No correction (or iteration disabled): single matrix
        Hci =
          CI::construct_Hci<SparseMatrix>(psi, h1, qk, br_ptr, s2_ptr, s1c, hk);
      }
      // Each stored element: value (double) and column index (32 bit)
      const auto Hci_sparse_GB = double(Hci.num_stored()) * 12.0 / 1.0e9;
      fmt::print(outstream,
                 "CI matrix (sparse): {} non-zero elements ({:.2f} GB)\n",
                 Hci.num_stored(), Hci_sparse_GB);
      if (all_below_cm) {
        fmt::print(outstream,
                   "Finding all solutions below {} cm^-1 (Davidson)\n",
                   *all_below_cm);
      } else if (num_solutions <= 0) {
        fmt::print(outstream, "Finding all solutions (Davidson)\n");
      } else {
        fmt::print(outstream, "Find first {} solutions (Davidson)\n",
                   num_solutions);
      }

      IO::ChronoTimer t2("", false);
      const auto convergence =
        psi.solve(Hci, num_solutions, all_below_cm, outstream);
      fmt::print(outstream,
                 "{} eigenvalues: {} iterations, eps = {:.1e}: T = {}\n",
                 psi.num_solutions(), convergence.iterations,
                 convergence.max_residual, t2.reading_str());
      if (!convergence.converged) {
        fmt::print(outstream,
                   "Warning: Davidson not converged (J={}{}): not saved\n",
                   J_string(twoJ), parity == 1 ? '+' : '-');
        write_solutions = false;
      }

    } else {
      LinAlg::Matrix<double> Hci;
      if (s1c != nullptr) {
        Hci = CI::iterate_E0(&psi, *s1c, h1, qk, br_ptr, s2_ptr, hk, dSk,
                             E0_sigma2, outstream);
      } else {
        // No correction (or iteration disabled): single matrix
        Hci = CI::construct_Hci(psi, h1, qk, br_ptr, s2_ptr, s1c, hk);
      }
      fmt::print(outstream,
                 "CI matrix (dense): {} total elements ({:.2f} GB)\n",
                 Hci.size(), double(Hci.size()) * 8.0 / 1.0e9);

      if (all_below_cm) {
        fmt::print(outstream, "Finding all solutions below {} cm^-1\n",
                   *all_below_cm);
      } else if (num_solutions > 0) {
        fmt::print(outstream, "Find first {} solutions\n", num_solutions);
      } else {
        fmt::print(outstream, "Finding all solutions\n");
      }

      IO::ChronoTimer t2("", false);
      psi.solve(Hci, num_solutions, all_below_cm);
      outstream << psi.num_solutions()
                << " eigenvalues: T = " << t2.reading_str() << "\n";
    }

    // Optionally, each level with its own E0, found in the same way (and
    // solved again with the Hamiltonian built at that E0, as above), so its
    // energy denominators (BW) and dSigma/dE correction are at its own
    // energy. Each level is then an eigenstate of a different Hamiltonian, so
    // the levels are not orthogonal
    if (BW_each_level && s1c != nullptr) {
      for (std::size_t i = 1; i < psi.num_solutions(); ++i) {
        auto s1_level = *s1c;
        s1_level.E0 = psi.energy(i);
        auto psi_level = psi;
        if (davidson) {
          const auto Hci =
            CI::iterate_E0<LinAlg::SparseSymmetricMatrix<double>>(
              &psi_level, s1_level, h1, qk, br_ptr, s2_ptr, hk, dSk, E0_sigma2,
              outstream, i);
          if (!psi_level.solve(Hci, int(i) + 1, {}, outstream).converged) {
            fmt::print(
              outstream,
              "Warning: Davidson not converged (level {}): not saved\n", i);
            write_solutions = false;
          }
        } else {
          const auto Hci =
            CI::iterate_E0(&psi_level, s1_level, h1, qk, br_ptr, s2_ptr, hk,
                           dSk, E0_sigma2, outstream, i);
          psi_level.solve(Hci, int(i) + 1);
        }
        psi.replace_solution(i, psi_level.energy(i), psi_level.coefs(i));
      }
    }

    if (!ci_fname.empty() && write_solutions) {
      psi.read_write(ci_fname, IO::FRW::write, outstream);
    }
  } else {
    outstream << "Using " << psi.num_solutions()
              << " solutions read from file.\n";
  }
  outstream << "\n";
  const auto E0 = psi.num_solutions() > 0 ? psi.energy(0) : 0.0;

  // For calculating g-factors
  DiracOperator::M1 m1{ci_sp_basis.front().grid(), PhysConst::alpha, 0.0};
  // only actually need to do this once..
  const auto m1_tab = Amplitudes::me_table(ci_sp_basis, &m1);

  // Print details of each solution, unless we find all, or read from file:
  const auto print_details_tmp =
    (all_below_cm || num_solutions > 0) && need_solve;
  print_details = print_details && print_details_tmp && !read_only;
  const double minimum_percentage = 5.0; // min % to print

  for (std::size_t i = 0; i < N_CSFs && i < psi.num_solutions(); ++i) {

    const auto pi = parity == 1 ? '+' : '-';
    if (print_details)
      fmt::print(outstream,
                 "{} {} {:<2}  {:+11.8f} au  {:+11.2f} cm^-1  {:11.2f} cm^-1\n",
                 J_string(twoJ), pi, i, psi.energy(i),
                 psi.energy(i) * PhysConst::Hartree_invcm,
                 (psi.energy(i) - E0) * PhysConst::Hartree_invcm);

    if (print_details) {
      for (std::size_t j = 0ul; j < N_CSFs; ++j) {
        const auto cj = 100.0 * std::pow(psi.coef(i, j), 2);
        if (cj > minimum_percentage) {
          fmt::print(outstream, "   {:<9s} {:5.3f}%\n", psi.csf(j).config(true),
                     cj);
        }
      }
    }
    // Leading non-relativistic configuration, and its |c|^2 weight:
    const auto [config, pc] = CI::leading_config(psi.coefs(i), psi.CSFs());

    // g_J <JJz|J|JJz> = <JJz|L + 2*S|JJz>
    // take J=Jz, <JJz|J|JJz> = J
    // then: g_J = <JJ|L + 2*S|JJ> / J
    // And: <JJ|L + 2*S|JJ> = 3js * <A||L+2S||A> (W.E. Theorem)
    const auto m1AA_R =
      CI::ReducedME(psi, i, psi, i, m1_tab, m1.rank(), m1.parity());
    const auto tjs = Angular::threej_2(twoJ, twoJ, 2, twoJ, -twoJ, 0);

    // Calculate g-factors, for line identification. Only defined for J!=0
    const double gJ = twoJ != 0 ? tjs * m1AA_R / (0.5 * twoJ) : 0.0;

    // <L^2> and <S^2> of the CI state (non-rel limit): measure LS-purity
    const auto [L2, S2] = CI::expectation_L2S2(psi.coefs(i), psi.CSFs(), twoJ);

    // Determine Term Symbol, from <L^2> and <S^2>
    const auto [twoS, L] = CI::Term_S_L_from_expectation(L2, S2, twoJ);
    const auto LSeff = [](double x) {
      return 0.5 * (std::sqrt(1.0 + 4.0 * x) - 1.0);
    };

    if (print_details) {
      outstream << "   --------------\n";
      if (twoJ != 0) {
        outstream << "   gJ = " << gJ << "\n";
      }
      // fmt::print(outstream, "   <L^2> = {:.4f}, <S^2> = {:.4f}\n", L2, S2);
      fmt::print(outstream, "   L_eff = {:.4f}, S_eff = {:.4f}\n", LSeff(L2),
                 LSeff(S2));
      fmt::print(outstream, "   {:<9s} {}\n", config,
                 Term_Symbol(twoJ, L, twoS, parity));
      outstream << "\n";
    }

    psi.update_config_info(i, {config, pc, gJ, 1.0 * L, 1.0 * twoS, L2, S2});
  }

  return psi;
}

} // namespace CI