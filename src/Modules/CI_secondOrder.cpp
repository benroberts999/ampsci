#include "Amplitudes/MatrixElements.hpp"
#include "Amplitudes/SecondOrderCI.hpp"
#include "Angular/include.hpp"
#include "CI/SecondOrder.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Modules/Modules.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/format.hpp"
#include <cmath>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
// Is this rank K allowed: triangle rules for the operators, and for the states
static bool allowed_K(int K, int kt, int ks, int twoJb, int twoJa) {
  return K >= std::abs(kt - ks) && K <= kt + ks &&
         Angular::triangle(twoJb, 2 * K, twoJa) != 0;
}

//==============================================================================
// The smallest rank K allowed for the amplitude; negative if there is none
static int smallest_allowed_K(int kt, int ks, int twoJb, int twoJa) {
  for (int K = std::abs(kt - ks); K <= kt + ks; ++K) {
    if (allowed_K(K, kt, ks, twoJb, twoJa))
      return K;
  }
  return -1;
}

//==============================================================================
// Short label for a CI state, e.g., "1-:0  6s6p 3P"
static std::string label(const CI::PsiJPi &psi, std::size_t index) {
  const auto &info = psi.info(index);
  return fmt::format("{:<6} {} {}",
                     CI::to_string({psi.twoJ(), psi.parity(), index}),
                     info.config,
                     CI::Term_Symbol(int(std::round(info.L)),
                                     int(std::round(info.twoS)), psi.parity()));
}

//==============================================================================
namespace Module {

// Declare, register, then define below.
void CI_secondOrder(const IO::InputBlock &input, const Wavefunction &wf);
namespace {
const Register r_CI_secondOrder{
  "CI_secondOrder",
  "Second-order amplitudes between CI states (polarisabilities, PNC), via CI "
  "mixed states",
  &CI_secondOrder};
} // namespace

void CI_secondOrder(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"",
      "Second-order amplitude A^K between two CI states, A -> B, for two "
      "operators t and s, at frequencies omega and omega_s:\n"
      "A^K = sum_n [c1 <B||t||n><n||s||A>/(E_A+omega_s-E_n) "
      "+ c2 <B||s||n><n||t||A>/(E_A+omega-E_n)].\n"
      "Energy conservation fixes omega_s = E_B - E_A - omega, so only omega is "
      "an input: at its default, t carries the whole transition frequency and "
      "s is static. For a dynamic polarisability, set B = A and omega to the "
      "frequency; then omega_s = -omega.\n"
      "The valence sums are evaluated either by sum-over-states over the "
      "solved CI levels (method=SOS; truncated to what the CI{} block "
      "solved), or with CI mixed states (method=MS), which is complete (no "
      "truncation). Run twice to compare the two. Requires a CI{} block "
      "(which is where the CI options are set)."},
     {"A", "Initial CI state, as J{+,-}:index - e.g., '0+' is the lowest J=0 "
           "even-parity solution, and '1-:2' the third J=1 odd-parity one. "
           "May be a comma separated list "
           "(A = 0+:0, 1-:2;) each pair (A_i -> B_i) is calculated in turn, "
           "with a summary of the totals at the end [required]"},
     {"B", "Final CI state(s), as for A [default: same as A]. If B is "
           "shorter than A, the remaining pairs are diagonal: A = a,b,c; "
           "B = d,e; gives (a,d), (b,e), (c,c)"},
     {"t", "The operator that carries the frequency omega [E1]"},
     {"t_options{}", "Options for the t operator"},
     {"s", "The other operator; it carries omega_s = E_b - E_a - omega [E1]"},
     {"s_options{}", "Options for the s operator"},
     {"omega", "Frequency of t. For a transition, leave as default, so that t "
               "carries it all and s is static. For a dynamic polarisability "
               "(B = A), set this to the frequency: s then carries -omega. "
               "[default: E_b - E_a]"},
     {"K", "Rank K of the amplitude. Requires |kt-ks| <= K <= kt+ks, and the "
           "triangle rule for (Jb,K,Ja) [default: smallest allowed]"},
     {"method", "SOS (sum-over-states over the solved CI levels) or MS (CI "
                "mixed states) [MS]"},
     {"rpa", "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram "
             "[true]"},
     {"project_out",
      "List of CI states, as for A, that are removed from the intermediate "
      "states - e.g., '1-:0, 1-:1' removes the two lowest J=1 odd solutions "
      "from the sums [none]"},
     {"StructureRadiation{}",
      "Options for structure radiation and normalisation. If this block is "
      "included, SR and Normalisation is added "
      "to every single-particle matrix element used in the amplitude"}});

  // Check for Structure Radiation
  const auto t_SR_input = input.getBlock("StructureRadiation");
  auto SR_input =
    t_SR_input ? *t_SR_input : IO::InputBlock{"StructureRadiation"};
  if (input.has_option("help")) {
    SR_input.add("help;");
  }
  SR_input.check(
    {{"", "If this block is included, structure radiation will be included"},
     {"Qk_file", "true/false/filename - SR: filename for QkTable file. If "
                 "blank will not use QkTable; if exists, will read it in; if "
                 "doesn't exist, will create it and write to disk. If 'true' "
                 "will use default filename"},
     {"n_minmax", "list; min,max n for core/excited (internal): [1,inf]"},
     {"n_max_legs",
      "Max n of the CI basis states that SR is applied to. SR is only "
      "meaningful between physical states, and the high-n states of the CI "
      "basis are cavity states [default: max_n_core + 3]"},
     {"norm", "Include the normalisation of states, applied to the "
              "single-particle matrix elements [true]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  using namespace std::string_literals;

  //----------------------------------------------------------------------------
  // The CI solutions, and the integrals used to construct them
  const auto &ints = wf.CI_integrals();
  if (wf.CIwfs().empty() || !ints.availableQ()) {
    fmt2::error();
    std::cout << ": Requires CI solutions and the integrals used to construct "
                 "the CI Hamiltonian. Include a CI{} block\n";
    return;
  }

  //----------------------------------------------------------------------------
  // The CI states: lists of A (and B), calculated pairwise (A_i -> B_i).
  // If B is shorter than A, the remaining pairs are diagonal
  const auto list_A = input.get("A", std::vector<std::string>{});
  const auto list_B = input.get("B", std::vector<std::string>{});
  if (list_A.empty()) {
    fmt2::error();
    std::cout << ": No A state(s) given\n";
    return;
  }
  if (list_B.size() > list_A.size()) {
    fmt2::warning();
    std::cout << ": B list is longer than A; the extra B states are ignored\n";
  }

  //----------------------------------------------------------------------------
  // The two operators
  const auto t_name = input.get("t", "E1"s);
  const auto s_name = input.get("s", "E1"s);
  const auto t_options = input.getBlock("t_options");
  const auto s_options = input.getBlock("s_options");
  const auto ht = DiracOperator::generate(
    t_name, t_options ? *t_options : IO::InputBlock(t_name, {}), wf);
  const auto hs = DiracOperator::generate(
    s_name, s_options ? *s_options : IO::InputBlock(s_name, {}), wf);
  const auto kt = ht->rank();
  const auto ks = hs->rank();

  const auto print_op = [](const std::string &op_label,
                           const DiracOperator::TensorOperator *h, double w) {
    fmt::print("{}: {} (rank {}, {} parity), ", op_label, h->name(), h->rank(),
               h->parity() == 1 ? "even" : "odd");
    if (w == 0.0) {
      fmt::print("static\n");
    } else {
      fmt::print("at omega = {:.6f}\n", w);
    }
  };

  //----------------------------------------------------------------------------
  // Method: sum-over-states (over the solved CI levels), or CI mixed states
  const auto method = input.get("method", "MS"s);
  const auto use_sos = qip::ci_compare(method, "SOS");
  if (!use_sos && !qip::ci_compare(method, "MS")) {
    fmt2::error();
    std::cout << ": Unknown method '" << method << "': use SOS or MS\n";
    return;
  }
  std::cout << (use_sos ?
                  "Method: sum-over-states (over the solved CI levels)\n" :
                  "Method: CI mixed states\n");

  //----------------------------------------------------------------------------
  // Levels to remove from the intermediate states
  std::vector<CI::Level> levels_to_remove;
  for (const auto &str : input.get("project_out", std::vector<std::string>{})) {
    const auto level = CI::parse_level(str);
    if (!level) {
      fmt2::error();
      fmt::print(
        ": Could not parse '{}' as a CI level: give as J{{+,-}}:index, "
        "e.g., 1-:0\n",
        str);
      return;
    }
    levels_to_remove.push_back(*level);
  }
  if (!levels_to_remove.empty()) {
    std::cout << "\nRemoving from the intermediate states:";
    for (const auto &level : levels_to_remove) {
      fmt::print(" {}", CI::to_string(level));
    }
    std::cout << "\n";
  }

  //----------------------------------------------------------------------------
  // RPA
  auto rpa_method = input.get("rpa", "true"s);
  if (wf.core().empty()) {
    rpa_method = "false";
  }
  std::cout << "\n";
  auto rpa_t = ExternalField::make_rpa(rpa_method, ht.get(), wf.vHF(), true,
                                       wf.basis(), wf.identity());
  auto rpa_s = ExternalField::make_rpa(rpa_method, hs.get(), wf.vHF(), false,
                                       wf.basis(), wf.identity());

  //----------------------------------------------------------------------------
  // Structure radiation. SR+N is only meaningful between physical states, so
  // it is applied only to the low-n part of the CI basis
  const auto sr_n_max =
    SR_input.get("n_max_legs", DiracSpinor::max_n(wf.core()) + 3);
  const auto sr_norm = SR_input.get("norm", true);

  std::optional<MBPT::StructureRad> sr;
  if (t_SR_input) {
    const auto n_minmax = SR_input.get("n_minmax", std::vector{1});
    const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
    const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
    const auto Qk_file_t = SR_input.get("Qk_file", "false"s);
    const std::string Qk_file =
      Qk_file_t != "false" ?
        (Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t) :
        "";
    std::cout << "\nIncluding structure radiation"
              << (sr_norm ? " and normalisation:\n" : " (no normalisation):\n");
    fmt::print("Applied to CI basis states with n <= {}\n", sr_n_max);
    sr =
      MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max}, Qk_file);
  }

  //----------------------------------------------------------------------------
  // Intermediate states carrying a core hole lie outside the CI space (used
  // for the core and core-valence terms) - common to every pair. They come
  // from the same source as the CI basis, so that the core-valence term
  // removes the excitations into the occupied CI orbitals
  const auto &core_hole_source = ints.brueckner ? wf.spectrum() : wf.basis();
  const auto excited =
    DiracSpinor::split_by_energy(core_hole_source, wf.FermiLevel()).second;

  // The RPA solutions and single-particle matrix-element tables depend only
  // on the frequencies: they are re-used between pairs when omega is
  // unchanged (e.g., a list of static polarisabilities)
  std::optional<std::pair<double, double>> table_freqs{};
  Coulomb::meTable<double> t_me{};
  Coulomb::meTable<double> s_me{};

  //----------------------------------------------------------------------------
  // Loop over the A -> B pairs; store the totals for the summary
  struct Summary {
    std::string A, B, name;
    int K;
    double A_K, value;
  };
  std::vector<Summary> summary;

  for (std::size_t i = 0; i < list_A.size(); ++i) {
    const auto &str_A = list_A.at(i);
    const auto &str_B = i < list_B.size() ? list_B.at(i) : str_A;

    if (list_A.size() > 1) {
      std::cout << "\n";
      IO::print_line('-', 40);
    }

    const auto level_A = CI::parse_level(str_A);
    const auto level_B = CI::parse_level(str_B);
    const auto Psi_a =
      level_A ? wf.CIwf(level_A->twoJ / 2, level_A->parity) : nullptr;
    const auto Psi_b =
      level_B ? wf.CIwf(level_B->twoJ / 2, level_B->parity) : nullptr;
    const auto ia = level_A ? level_A->index : 0;
    const auto ib = level_B ? level_B->index : 0;
    if (Psi_a == nullptr || Psi_b == nullptr || ia >= Psi_a->num_solutions() ||
        ib >= Psi_b->num_solutions()) {
      fmt2::error();
      std::cout << ": Could not find the requested CI state(s): A = " << str_A
                << ", B = " << str_B
                << "\nGive as J{+,-}:index; e.g., A = 0+; B = 1-:2;\n"
                   "Available (J,parity):";
      for (const auto &psi : wf.CIwfs()) {
        fmt::print(" {}{}", psi.twoJ() / 2, psi.parity() == 1 ? '+' : '-');
      }
      std::cout << "\n";
      continue;
    }
    const auto Ea = Psi_a->energy(ia);
    const auto Eb = Psi_b->energy(ib);

    fmt::print("\nSecond-order amplitude, {} -> {}:\n",
               CI::to_string({Psi_a->twoJ(), Psi_a->parity(), ia}),
               CI::to_string({Psi_b->twoJ(), Psi_b->parity(), ib}));
    fmt::print("A: {}  E = {:.8f} au\n", label(*Psi_a, ia), Ea);
    fmt::print("B: {}  E = {:.8f} au\n", label(*Psi_b, ib), Eb);

    const auto omega = input.get("omega", Eb - Ea);

    // s carries whatever frequency t does not, so that energy is conserved:
    // E_b = E_a + omega + omega_s. With omega at its default, s is static, as
    // for a real transition. For a dynamic polarisability, A = B (so
    // E_b = E_a) and omega is set: then omega_s = -omega, giving the usual
    // pair of denominators, E_a -+ omega - E_n
    const auto omega_s = Eb - Ea - omega;
    const auto dynamic = omega != 0.0 && omega_s != 0.0;

    print_op("t", ht.get(), omega);
    print_op("s", hs.get(), omega_s);

    // Overall parity selection rule
    if (Psi_a->parity() * Psi_b->parity() != ht->parity() * hs->parity()) {
      std::cout << ": Amplitude is zero by parity\n";
      continue;
    }

    // Rank K: default is the smallest allowed by the triangle rules
    const auto K_minimum =
      smallest_allowed_K(kt, ks, Psi_b->twoJ(), Psi_a->twoJ());
    if (K_minimum < 0) {
      std::cout << ": No allowed K. Require |kt-ks| <= K <= kt+ks, and the "
                   "triangle rule for (Jb, K, Ja)\n";
      continue;
    }
    const auto K = input.get("K", K_minimum);
    fmt::print("K = {}{}\n", K, K == K_minimum ? " (smallest allowed)" : "");
    if (!allowed_K(K, kt, ks, Psi_b->twoJ(), Psi_a->twoJ())) {
      std::cout << ": K is not allowed: the amplitude is zero\n";
      continue;
    }

    //--------------------------------------------------------------------------
    // RPA, and single-particle matrix element tables, at this pair's
    // frequencies; skipped if they are unchanged from the previous pair
    if (!table_freqs || table_freqs->first != omega ||
        table_freqs->second != omega_s) {
      if (ht->freqDependantQ()) {
        ht->updateFrequency(omega);
      }
      if (hs->freqDependantQ()) {
        hs->updateFrequency(omega_s);
      }
      if (rpa_t) {
        std::cout << "Solve RPA for t at omega = " << omega << "\n";
        rpa_t->solve_core(omega);
      }
      if (rpa_s) {
        std::cout << "Solve RPA for s at omega = " << omega_s << "\n";
        rpa_s->solve_core(omega_s);
      }

      std::cout << "\nFill matrix element tables..." << std::flush;
      if (sr) {
        sr->solve_core(ht.get(), rpa_t.get());
      }
      // The normalisation of states goes onto the single-particle matrix
      // elements, as in the single-valence case
      t_me =
        Amplitudes::me_table(ints.ci_basis, ht.get(), rpa_t.get(),
                             sr ? &*sr : nullptr, omega, sr_n_max, sr_norm);
      if (sr) {
        sr->solve_core(hs.get(), rpa_s.get());
      }
      s_me =
        Amplitudes::me_table(ints.ci_basis, hs.get(), rpa_s.get(),
                             sr ? &*sr : nullptr, omega_s, sr_n_max, sr_norm);
      std::cout << "done\n" << std::flush;
      table_freqs = std::pair{omega, omega_s};
    }

    //--------------------------------------------------------------------------
    // The amplitude: the valence (CI) sum, by SOS or mixed states
    double A_val{0.0};
    if (use_sos) {
      fmt::print("\nContributions to A^{}, by intermediate J and parity:\n", K);
      A_val =
        Amplitudes::sos_ci(K, *Psi_b, ib, *Psi_a, ia, ht.get(), t_me, hs.get(),
                           s_me, omega, omega_s, wf.CIwfs(), levels_to_remove);
    } else {
      fmt::print("\nContributions to A^{}, by intermediate J and parity:\n", K);
      const auto [A_s, A_t] =
        CI::A_K(K, *Psi_b, ib, *Psi_a, ia, ht.get(), t_me, hs.get(), s_me,
                omega, omega_s, ints, levels_to_remove, wf.CIwfs());
      A_val = A_s;
    }

    //--------------------------------------------------------------------------
    // Intermediate states carrying a core hole. These lie outside the CI
    // space, so the mixed states cannot produce them: the polarisation of the
    // core (K=0 and diagonal only), and the core-valence term, which is the
    // Pauli blocking of the core excitations by the valence electrons

    // The core is the same in A and B, so its amplitude needs <B|A>
    const auto diagonal = Psi_a == Psi_b && ia == ib;
    const auto A_core =
      diagonal ?
        CI::A_K_core(K, Psi_a->twoJ(), ht.get(), hs.get(), omega, omega_s,
                     wf.core(), excited, rpa_t.get(), rpa_s.get()) :
        0.0;
    const auto A_cv =
      CI::A_K_cv(K, *Psi_b, ib, *Psi_a, ia, ht.get(), hs.get(), omega, omega_s,
                 wf.core(), ints.ci_basis, rpa_t.get(), rpa_s.get());

    // SR and normalisation enter the valence sum through the single-particle
    // tables; the core and core-valence terms are bare
    const auto A_total = A_val + A_cv + A_core;

    //--------------------------------------------------------------------------
    // The second column of the table is the quantity that was asked for. For
    // the specific cases below (as in the dcp module) that is
    // alpha/beta/E_pnc; otherwise it is the rank-K part of the z-component,
    // A_zz. Only s is tested: t is E1 in all of these cases
    const auto E1_s = hs->name() == "E1";
    // e.g., "Scalar polarisability", "Scalar dynamic polarisability",
    // "Scalar transition polarisability"
    const auto kind =
      " "s + (dynamic ? "dynamic "s : ""s) + (diagonal ? ""s : "transition "s);

    // The z-component of the rank-K amplitude: m_a = m_b = m, and q = 0 for
    // both operators
    const auto two_m = std::min(Psi_a->twoJ(), Psi_b->twoJ());
    auto title = fmt::format("z-component, m = {}", two_m / 2);
    auto name = fmt::format("A_zz (K={})", K);
    auto factor =
      CI::z_component(K, kt, ks, Psi_b->twoJ(), Psi_a->twoJ(), two_m);

    // Scalar polarisability
    if (K == 0 && E1_s) {
      // alpha_0 = (2/3)[J]^-1 sum_n |<a||d||n>|^2/(E_n-E_a), for a = b
      title = fmt::format("Scalar{}polarisability", kind);
      name = "alpha (au)";
      factor = 1.0 / std::sqrt(3.0 * (Psi_b->twoJ() + 1));
    }

    // Tensor polarisability. The normalisation comes from the m-dependence of
    // the Stark shift, so it requires J_a = J_b (as the scalar part does,
    // though there K=0 already implies it), and J >= 1
    if (K == 2 && E1_s && Psi_a->twoJ() == Psi_b->twoJ() &&
        Psi_b->twoJ() >= 2) {
      const auto twoJ = double(Psi_b->twoJ());
      title = fmt::format("Tensor{}polarisability", kind);
      name = "alpha_2 (au)";
      factor = -std::sqrt(2.0 * twoJ * (twoJ - 1.0) /
                          (3.0 * (twoJ + 1.0) * (twoJ + 2.0) * (twoJ + 3.0)));
    }

    // Vector polarisability of a state, alpha_1 = A^1 sqrt(2J/((J+1)(2J+1))):
    // the Stark shift is -(1/2)|E|^2 [alpha_0 + A cos(theta_k) (m/2J) alpha_1
    // + ...], for light of degree of circular polarisation A, as in Le Kien,
    // Schneider, Rauschenbeutel, Eur. Phys. J. D 67, 92 (2013), whose
    // alpha^(1) is -A^1
    if (K == 1 && E1_s && diagonal) {
      const auto J = double(Psi_a->twoJ()) / 2.0;
      title = fmt::format("Vector{}polarisability", kind);
      name = "alpha_1 (au)";
      factor = std::sqrt(2.0 * J / ((J + 1.0) * (2.0 * J + 1.0)));
    }

    // Vector transition polarisability, beta = A^1/(sqrt(2) <b||sigma||a>).
    // The spin matrix element expresses the Wigner-Eckart factor of the
    // rank-1 amplitude, so no radial overlap enters it; see CI::sigma_rme
    if (K == 1 && E1_s && !diagonal) {
      const auto sigma = CI::sigma_rme(*Psi_b, ib, *Psi_a, ia, ints.ci_basis);
      fmt::print("\nVector{}polarisability:\n", kind);
      fmt::print("<B||sigma||A> = {:.6e}\n", sigma);
      if (std::abs(sigma) < 1.0e-12) {
        std::cout << "beta: not defined - the states have no spin-angular "
                     "structure in common\n";
      } else {
        title = fmt::format("Vector{}polarisability", kind);
        name = "beta (au)";
        factor = 1.0 / (std::sqrt(2.0) * sigma);
      }
    }

    // PNC amplitude: the static operator is the PNC interaction
    if (hs->name().substr(0, 3) == "pnc") {
      title = "PNC amplitude";
      name = fmt::format("E_pnc ({})", hs->units());
    }

    //--------------------------------------------------------------------------
    fmt::print("\n{}:\n", title);
    fmt::print("{:<14} {:>16} {:>16}\n", "", fmt::format("A^{}", K), name);
    fmt::print("{:<14} {:16.6e} {:16.6e}\n", "valence (CI)", A_val,
               factor * A_val);
    fmt::print("{:<14} {:16.6e} {:16.6e}\n", "core", A_core, factor * A_core);
    fmt::print("{:<14} {:16.6e} {:16.6e}\n", "core-valence", A_cv,
               factor * A_cv);
    fmt::print("{:<14} {:16.6e} {:16.6e}\n", "total", A_total,
               factor * A_total);

    summary.push_back({CI::to_string({Psi_a->twoJ(), Psi_a->parity(), ia}),
                       CI::to_string({Psi_b->twoJ(), Psi_b->parity(), ib}),
                       name, K, A_total, factor * A_total});
  }

  //----------------------------------------------------------------------------
  // Summary of the totals, for a list of pairs
  if (list_A.size() > 1 && !summary.empty()) {
    std::cout << "\n";
    IO::print_line('-', 40);
    std::cout << "\nSummary of totals:\n";
    fmt::print("{:>6} {:>6} {:>3} {:>16} {:>16}\n", "A", "B", "K", "A^K",
               "value");
    for (const auto &row : summary) {
      fmt::print("{:>6} {:>6} {:>3} {:16.6e} {:16.6e}  {}\n", row.A, row.B,
                 row.K, row.A_K, row.value, row.name);
    }
  }
}

} // namespace Module
