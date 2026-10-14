#include "Amplitudes/SecondOrder.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Angular/include.hpp"
#include "CI/SecondOrder.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/GenerateOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Modules/Modules.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace Module {

// Declare, register, then define below.
void secondOrder(const IO::InputBlock &input, const Wavefunction &wf);
namespace {
const Register r_secondOrder{
  "secondOrder",
  "Second-order amplitudes for a valence state/transition (polarisabilities, "
  "PNC), by sum-over-states or mixed states",
  &secondOrder};
} // namespace

//==============================================================================
void secondOrder(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"",
      "Second-order amplitude A^K for a single-valence transition A -> B, for "
      "two operators t and s, at frequencies omega and omega_s:\n"
      "A^K = sum_n [c1 <B||t||n><n||s||A>/(E_A+omega_s-E_n) "
      "+ c2 <B||s||n><n||t||A>/(E_A+omega-E_n)].\n"
      "Energy conservation fixes omega_s = E_B - E_A - omega, so only omega is "
      "an input: at its default, t carries the whole transition frequency and "
      "s is static. For a dynamic polarisability, set B = A and omega to the "
      "frequency; then omega_s = -omega.\n"
      "The valence sum is evaluated either by sum-over-states (method=SOS), "
      "or with mixed states (method=MS), which is complete (no truncation of "
      "the sum). Run twice to compare the two."},
     {"A", "Initial valence state, e.g., 6s+. May be a list (A = 6s+,7s+;): "
           "each pair (A_i -> B_i) is calculated in turn, with a summary of "
           "the totals at the end [required]"},
     {"B", "Final valence state(s) [default: same as A]. If B is shorter "
           "than A, the remaining pairs are diagonal: A = a,b,c; B = d,e; "
           "gives (a,d), (b,e), (c,c)"},
     {"t", "The operator that carries the frequency omega [E1]"},
     {"t_options{}", "Options for the t operator"},
     {"s", "The other operator; it carries omega_s = E_b - E_a - omega [E1]"},
     {"s_options{}", "Options for the s operator"},
     {"omega", "Frequency of t. For a transition, leave as default, so that t "
               "carries it all and s is static. For a dynamic polarisability "
               "(B = A), set this to the frequency: s then carries -omega. "
               "[default: E_b - E_a]"},
     {"omega_t", "Explicit frequency for t, overriding omega. Energy "
                 "conservation is then up to you; use with care [rare]"},
     {"omega_s", "Explicit frequency for s, overriding omega_s = "
                 "E_b - E_a - omega. Use with care [rare]"},
     {"K", "Rank K of the amplitude. Requires |kt-ks| <= K <= kt+ks, and the "
           "triangle rule for (jb,K,ja) [default: smallest allowed]"},
     {"method", "SOS (sum-over-states) or MS (mixed states) [MS]"},
     {"rpa",
      "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram, cntm "
      "[true]. MS requires TDHF: the diagram method triggers a "
      "warning, and TDHF is used"},
     {"replace_w_valence",
      "SOS: replace the basis/spectrum states with the corresponding valence "
      "states in the sum, where available (e.g., when the valence states were "
      "fitted to experiment) [false]"},
     {"n_main", "The 'main' part of the valence sum: intermediate states up "
                "to this n, printed separately (as the pnc module) "
                "[max_n_core + 4]"},
     {"StructureRadiation{}",
      "Options for structure radiation and normalisation. If this block is "
      "included, SR+N is added to the single-particle matrix elements of the "
      "sum. Sum-over-states method only"}});

  // Check for Structure Radiation
  const auto t_SR_input = input.getBlock("StructureRadiation");
  auto SR_input =
    t_SR_input ? *t_SR_input : IO::InputBlock{"StructureRadiation"};
  if (input.has_option("help")) {
    SR_input.add("help;");
  }
  SR_input.check(
    {{"", "If this block is included, SR + Normalisation corrections will be "
          "included (sum-over-states method only)"},
     {"Qk_file", "true/false/filename - SR: filename for QkTable file. If "
                 "blank will not use QkTable; if exists, will read it in; if "
                 "doesn't exist, will create it and write to disk. If 'true' "
                 "will use default filename"},
     {"n_minmax", "list; min,max n for core/excited (internal): [1,inf]"},
     {"n_max_legs",
      "SR+N is applied to matrix elements whose states both have n <= this "
      "(the valence legs, and the low-n intermediate states). SR+N is only "
      "meaningful between physical states: the high-n basis states are "
      "cavity states [default: max_n_core + 3]"},
     {"norm", "Include the normalisation of states? If false, only the "
              "structure radiation is included [true]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  using namespace std::string_literals;

  //----------------------------------------------------------------------------
  // The valence states: lists of A (and B), calculated pairwise (A_i -> B_i).
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

  const auto print_op = [](const std::string &label,
                           const DiracOperator::TensorOperator *h, double w) {
    fmt::print("{}: {} (rank {}, {} parity), ", label, h->name(), h->rank(),
               h->parity() == 1 ? "even" : "odd");
    if (w == 0.0) {
      fmt::print("static\n");
    } else {
      fmt::print("at omega = {:.6f}\n", w);
    }
  };

  //----------------------------------------------------------------------------
  // Method: sum-over-states, or mixed states
  const auto method = input.get("method", "MS"s);
  const auto use_sos = qip::ci_compare(method, "SOS");
  if (!use_sos && !qip::ci_compare(method, "MS")) {
    fmt2::error();
    std::cout << ": Unknown method '" << method << "': use SOS or MS\n";
    return;
  }
  std::cout << (use_sos ? "Method: sum-over-states\n" :
                          "Method: mixed states\n");

  //----------------------------------------------------------------------------
  // RPA
  auto rpa_method = input.get("rpa", "true"s);
  if (wf.core().empty()) {
    rpa_method = "false";
  }
  const auto rpaQ =
    ExternalField::ParseMethod(rpa_method) != ExternalField::Method::none;

  std::cout << "\n";
  auto rpa_t = ExternalField::make_rpa(rpa_method, ht.get(), wf.vHF(), true,
                                       wf.basis(), wf.identity());
  auto rpa_s = ExternalField::make_rpa(rpa_method, hs.get(), wf.vHF(), false,
                                       wf.basis(), wf.identity());

  // Mixed states require TDHF (it provides the mixed-state solver, whether or
  // not RPA is on). If another RPA method was requested, warn and use TDHF
  auto tdhf_t = dynamic_cast<ExternalField::TDHF *>(rpa_t.get());
  auto tdhf_s = dynamic_cast<ExternalField::TDHF *>(rpa_s.get());
  if (!use_sos && (tdhf_t == nullptr || tdhf_s == nullptr)) {
    if (rpaQ) {
      fmt2::warning();
      std::cout << ": Mixed states require the TDHF method for RPA; using "
                   "TDHF (not "
                << rpa_method << ")\n";
    }
    rpa_t = std::make_unique<ExternalField::TDHF>(ht.get(), wf.vHF());
    rpa_s = std::make_unique<ExternalField::TDHF>(hs.get(), wf.vHF());
    tdhf_t = static_cast<ExternalField::TDHF *>(rpa_t.get());
    tdhf_s = static_cast<ExternalField::TDHF *>(rpa_s.get());
  }

  //----------------------------------------------------------------------------
  // SOS: the spectrum (or basis) that is summed over - common to every pair
  std::vector<DiracSpinor> spectrum{};
  if (use_sos) {
    const auto use_spectrum = !wf.spectrum().empty();
    spectrum = use_spectrum ? wf.spectrum() : wf.basis();
    if (spectrum.empty()) {
      fmt2::error();
      std::cout << ": Sum-over-states requires a basis (or spectrum)\n";
      return;
    }
    std::cout << "SOS: summing over the "
              << (use_spectrum ? "spectrum" : "basis") << " ("
              << DiracSpinor::state_config(spectrum) << ")\n";

    // Optionally, use valence states in place of the corresponding basis
    // states (e.g., when the valence states were fitted to experiment).
    // Not obviously more or less accurate: the downside is orthogonality
    if (input.get("replace_w_valence", false)) {
      std::cout
        << "Replacing spectrum states with corresponding valence states\n";
      for (const auto &v : wf.valence()) {
        auto pv = std::find(spectrum.begin(), spectrum.end(), v);
        if (pv != spectrum.end()) {
          *pv = v;
        }
      }
    }
  }

  // MS: the correlation potential for the mixed states (if available)
  const auto *Sigma = use_sos ? nullptr : wf.Sigma();
  if (Sigma) {
    std::cout << "MS: including the correlation potential in the mixed "
                 "states\n";
  }

  //----------------------------------------------------------------------------
  // Structure radiation + normalisation (sum-over-states only). Constructed
  // once; the matrix-element tables are built per pair
  std::optional<MBPT::StructureRad> sr{};
  auto sr_n_max = 0;
  auto sr_norm = true;
  if (t_SR_input && !use_sos) {
    fmt2::warning();
    std::cout << ": Structure radiation is not (simply) possible with mixed "
                 "states: ignoring. Use method=SOS for SR+N; the correction "
                 "may be added to the MS result by hand, as "
                 "SOS(with SR) - SOS(without)\n";
  }
  if (t_SR_input && use_sos) {
    sr_n_max = SR_input.get("n_max_legs", DiracSpinor::max_n(wf.core()) + 3);
    sr_norm = SR_input.get("norm", true);
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
    fmt::print("Applied to matrix elements with n <= {}\n", sr_n_max);
    sr =
      MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max}, Qk_file);
  }

  const auto n_main = input.get("n_main", DiracSpinor::max_n(wf.core()) + 4);

  //----------------------------------------------------------------------------
  // Loop over the A -> B pairs; store the totals for the summary
  struct Summary {
    std::string A, B, name;
    int K;
    double A_K, value;
  };
  std::vector<Summary> summary;

  for (std::size_t i = 0; i < list_A.size(); ++i) {
    const auto &name_A = list_A.at(i);
    const auto &name_B = i < list_B.size() ? list_B.at(i) : name_A;

    if (list_A.size() > 1) {
      std::cout << "\n";
      IO::print_line('-', 40);
    }

    // Valence states only: getState() would also find core states
    const auto [nA, kA] = AtomData::parse_symbol(name_A);
    const auto [nB, kB] = AtomData::parse_symbol(name_B);
    const auto pFa = wf.isInValence(nA, kA) ? wf.getState(nA, kA) : nullptr;
    const auto pFb = wf.isInValence(nB, kB) ? wf.getState(nB, kB) : nullptr;
    if (pFa == nullptr || pFb == nullptr) {
      fmt2::error();
      std::cout << ": Could not find the requested valence state(s): A = "
                << name_A << ", B = " << name_B << "\nValence:";
      for (const auto &v : wf.valence()) {
        std::cout << " " << v.shortSymbol();
      }
      std::cout << "\n";
      continue;
    }
    const auto &Fa = *pFa;
    const auto &Fb = *pFb;

    fmt::print("\nSecond-order amplitude, {} -> {}:\n", Fa.shortSymbol(),
               Fb.shortSymbol());
    fmt::print("A: {}  E = {:.8f} au\n", Fa.shortSymbol(), Fa.en());
    fmt::print("B: {}  E = {:.8f} au\n", Fb.shortSymbol(), Fb.en());

    // omega_t/omega_s: rare explicit overrides (e.g. comparison to
    // calculations with approximate frequencies); by default, energy
    // conservation holds
    const auto omega =
      input.get("omega_t", input.get("omega", Fb.en() - Fa.en()));
    const auto omega_s = input.get("omega_s", Fb.en() - Fa.en() - omega);
    const auto dynamic = omega != 0.0 && omega_s != 0.0;
    if (std::abs(Fb.en() - Fa.en() - omega - omega_s) > 1.0e-10) {
      std::cout << "Note: omega_t + omega_s does not equal E_b - E_a "
                   "(explicit override)\n";
    }

    print_op("t", ht.get(), omega);
    print_op("s", hs.get(), omega_s);

    // Overall parity selection rule
    if (Fa.parity() * Fb.parity() != ht->parity() * hs->parity()) {
      std::cout << ": Amplitude is zero by parity\n";
      continue;
    }

    // Rank K: default is the smallest allowed by the triangle rules
    const auto K_minimum =
      Amplitudes::smallest_allowed_K(kt, ks, Fb.twoj(), Fa.twoj());
    if (K_minimum < 0) {
      std::cout << ": No allowed K. Require |kt-ks| <= K <= kt+ks, and the "
                   "triangle rule for (jb, K, ja)\n";
      continue;
    }
    const auto K = input.get("K", K_minimum);
    fmt::print("K = {}{}\n", K, K == K_minimum ? " (smallest allowed)" : "");
    if (!Amplitudes::allowed_K(K, kt, ks, Fb.twoj(), Fa.twoj())) {
      std::cout << ": K is not allowed: the amplitude is zero\n";
      continue;
    }

    // Solve RPA at this pair's frequencies
    if (ht->freqDependantQ()) {
      ht->updateFrequency(omega);
    }
    if (hs->freqDependantQ()) {
      hs->updateFrequency(omega_s);
    }
    if (rpaQ && rpa_t) {
      std::cout << "Solve RPA for t at omega = " << omega << "\n";
      rpa_t->solve_core(omega);
    }
    if (rpaQ && rpa_s) {
      std::cout << "Solve RPA for s at omega = " << omega_s << "\n";
      rpa_s->solve_core(omega_s);
    }

    //--------------------------------------------------------------------------
    // The amplitude: valence part, and (for a diagonal amplitude) the
    // polarisation of the closed core
    const auto diagonal = (pFa == pFb);
    double A_v{0.0};
    double A_cv{0.0};
    double A_main{0.0};
    bool have_main{false};
    double A_core{0.0};

    if (use_sos) {
      // Structure radiation + normalisation: added to the single-particle
      // matrix elements, for pairs with n <= n_max_legs (the valence legs and
      // the low-n intermediate states, which are the physical ones)
      Coulomb::meTable<double> t_me{};
      Coulomb::meTable<double> s_me{};
      if (sr) {
        std::cout << "Fill matrix element tables..." << std::flush;
        const std::vector<DiracSpinor> legs{Fa, Fb};
        sr->solve_core(ht.get(), rpa_t.get());
        t_me = Amplitudes::me_table(legs, spectrum, ht.get(), rpa_t.get(), &*sr,
                                    omega, sr_n_max, sr_norm);
        sr->solve_core(hs.get(), rpa_s.get());
        s_me = Amplitudes::me_table(legs, spectrum, hs.get(), rpa_s.get(), &*sr,
                                    omega_s, sr_n_max, sr_norm);
        std::cout << "done\n" << std::flush;
      }

      A_v =
        Amplitudes::sos_valence(K, Fb, Fa, ht.get(), hs.get(), omega, omega_s,
                                spectrum, rpa_t.get(), rpa_s.get(), t_me, s_me);

      // Core-valence (Pauli blocking) part: the below-Fermi terms of the sum.
      // Already included in the valence total; separated for comparison
      const auto [below_fermi, excited] =
        DiracSpinor::split_by_energy(spectrum, wf.FermiLevel());
      A_cv = Amplitudes::sos_valence(K, Fb, Fa, ht.get(), hs.get(), omega,
                                     omega_s, below_fermi, rpa_t.get(),
                                     rpa_s.get(), t_me, s_me);

      // 'Main' part: the low-n (physical) excited states, as the pnc module
      std::vector<DiracSpinor> main_states;
      for (const auto &n : excited) {
        if (n.n() <= n_main) {
          main_states.push_back(n);
        }
      }
      A_main = Amplitudes::sos_valence(K, Fb, Fa, ht.get(), hs.get(), omega,
                                       omega_s, main_states, rpa_t.get(),
                                       rpa_s.get(), t_me, s_me);
      have_main = true;

      if (diagonal && K == 0) {
        A_core =
          Amplitudes::sos_core(K, Fa.twoj(), ht.get(), hs.get(), omega, omega_s,
                               wf.core(), excited, rpa_t.get(), rpa_s.get());
      }
    } else {

      const auto [A_ms_s, A_ms_t] = Amplitudes::ms_valence(
        K, Fb, Fa, ht.get(), hs.get(), omega, omega_s, tdhf_t, tdhf_s, Sigma);
      const auto eps =
        std::abs((A_ms_s - A_ms_t) / std::max(std::abs(A_ms_s), 1.0e-30));
      fmt::print("\n{:>4}  {:>16} {:>16} {:>8}\n", "", "<B_s||t||A_s>",
                 "<B_t||s||A_t>", "eps");
      fmt::print("{:>4}  {:16.6e} {:16.6e} {:8.1e}\n", fmt::format("A^{}", K),
                 A_ms_s, A_ms_t, eps);
      A_v = A_ms_s;

      // Core-valence (Pauli blocking) part: mixed states projected onto the
      // HF core. Already included in the valence total; separated for
      // comparison
      const auto [cv_s, cv_t] =
        Amplitudes::ms_valence(K, Fb, Fa, ht.get(), hs.get(), omega, omega_s,
                               tdhf_t, tdhf_s, Sigma, wf.core());
      A_cv = cv_s;

      // 'Main' part: mixed states projected onto the low-n excited states of
      // the spectrum (or basis), as the pnc module. Requires a basis
      const auto &sob = wf.spectrum().empty() ? wf.basis() : wf.spectrum();
      if (!sob.empty()) {
        std::vector<DiracSpinor> main_states;
        for (const auto &n :
             DiracSpinor::split_by_energy(sob, wf.FermiLevel()).second) {
          if (n.n() <= n_main) {
            main_states.push_back(n);
          }
        }
        A_main =
          Amplitudes::ms_valence(K, Fb, Fa, ht.get(), hs.get(), omega, omega_s,
                                 tdhf_t, tdhf_s, Sigma, main_states)
            .first;
        have_main = true;
      }

      if (diagonal && K == 0) {
        // Core excitations: valence Sigma is not appropriate here (cf
        // sos_core, which uses the HF basis)
        A_core = Amplitudes::ms_core(K, Fa.twoj(), ht.get(), hs.get(), omega,
                                     omega_s, wf.core(), tdhf_t, tdhf_s);
      }
    }

    const auto A_total = A_v + A_core;

    //--------------------------------------------------------------------------
    // The second column of the table is the quantity that was asked for. For
    // the specific cases below (as in the polarisability/pnc/dcp modules)
    // that is alpha/beta/E_pnc; otherwise it is the rank-K part of the
    // z-component, A_zz. Only s is tested: t is E1 in all of these cases
    const auto E1_s = hs->name() == "E1";
    const auto kind =
      " "s + (dynamic ? "dynamic "s : ""s) + (diagonal ? ""s : "transition "s);

    // The z-component of the rank-K amplitude: m_a = m_b = m, and q = 0 for
    // both operators
    const auto two_m = std::min(Fa.twoj(), Fb.twoj());
    auto title = fmt::format("z-component, m = {}/2", two_m);
    auto name = fmt::format("A_zz (K={})", K);
    auto factor = CI::z_component(K, kt, ks, Fb.twoj(), Fa.twoj(), two_m);

    // Scalar polarisability
    if (K == 0 && E1_s) {
      title = fmt::format("Scalar{}polarisability", kind);
      name = "alpha (au)";
      factor = 1.0 / std::sqrt(3.0 * (Fb.twoj() + 1));
    }

    // Tensor polarisability; requires j_a = j_b >= 1
    if (K == 2 && E1_s && Fa.twoj() == Fb.twoj() && Fb.twoj() >= 2) {
      const auto twoJ = double(Fb.twoj());
      title = fmt::format("Tensor{}polarisability", kind);
      name = "alpha_2 (au)";
      factor = -std::sqrt(2.0 * twoJ * (twoJ - 1.0) /
                          (3.0 * (twoJ + 1.0) * (twoJ + 2.0) * (twoJ + 3.0)));
    }

    // Vector polarisability of a state, alpha_1 = A^1 sqrt(2j/((j+1)(2j+1)));
    // see the CI_secondOrder module for the convention
    if (K == 1 && E1_s && diagonal) {
      const auto j = double(Fa.twoj()) / 2.0;
      title = fmt::format("Vector{}polarisability", kind);
      name = "alpha_1 (au)";
      factor = std::sqrt(2.0 * j / ((j + 1.0) * (2.0 * j + 1.0)));
    }

    // Vector transition polarisability, beta = A^1/(sqrt(2) <b||sigma||a>).
    // Single-valence convention: the radial overlap is dropped, so
    // <b||sigma||a> = 2 S_kk (cf CI::sigma_rme and the dcp module)
    if (K == 1 && E1_s && !diagonal) {
      const auto sigma = 2.0 * Angular::S_kk(Fb.kappa(), Fa.kappa());
      fmt::print("\nVector{}polarisability:\n", kind);
      fmt::print("<B||sigma||A> = {:.6e}\n", sigma);
      if (std::abs(sigma) < 1.0e-12) {
        std::cout
          << "beta: not defined - no spin-angular structure in common\n";
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
    fmt::print("{:<18}{:>16} {:>16}\n", "", fmt::format("A^{}", K), name);
    fmt::print("{:<18}{:16.6e} {:16.6e}\n", "valence", A_v, factor * A_v);
    fmt::print("    {:<14}{:16.6e} {:16.6e}   [included in valence]\n",
               "core-valence", A_cv, factor * A_cv);
    if (have_main) {
      fmt::print("    {:<14}{:16.6e} {:16.6e}   [included in valence]\n",
                 fmt::format("main (n<={})", n_main), A_main, factor * A_main);
      const auto A_tail = A_v - A_cv - A_main;
      fmt::print("    {:<14}{:16.6e} {:16.6e}   [included in valence]\n",
                 "tail", A_tail, factor * A_tail);
    }
    if (diagonal && K == 0) {
      fmt::print("{:<18}{:16.6e} {:16.6e}\n", "core", A_core, factor * A_core);
    }
    fmt::print("{:<18}{:16.6e} {:16.6e}\n", "total", A_total, factor * A_total);

    summary.push_back(
      {Fa.shortSymbol(), Fb.shortSymbol(), name, K, A_total, factor * A_total});
  }

  //----------------------------------------------------------------------------
  // Summary of the totals, for a list of pairs
  if (list_A.size() > 1 && !summary.empty()) {
    std::cout << "\n";
    IO::print_line('-', 40);
    std::cout << "\nSummary of totals:\n";
    fmt::print("{:>5} {:>6} {:>3} {:>16} {:>16}\n", "A", "B", "K", "A^K",
               "value");
    for (const auto &row : summary) {
      fmt::print("{:>5} {:>6} {:>3} {:16.6e} {:16.6e}  {}\n", row.A, row.B,
                 row.K, row.A_K, row.value, row.name);
    }
  }
}

} // namespace Module
