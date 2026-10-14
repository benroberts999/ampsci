#include "Amplitudes/Normalisation.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "DiracOperator/GenerateOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Modules/Modules.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace Module {

void normalisation(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_normalisation{
  "normalisation",
  "Normalisation correction to matrix elements, via derivative of Sigma",
  &normalisation};
} // namespace

//==============================================================================
// The matrix elements are calculated by Amplitudes::matrix_elements and the
// Sigma derivative by Amplitudes::dSigma_dE; this module only parses input
// and prints.
void normalisation(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check(
    {{"", "Calculates normalisation correction, via derivative of "
          "Correlation potential. Uses correlation potential from main "
          "wavefunction"},
     {"operator", "e.g., E1, hfs (see ampsci -o for available operators)"},
     {"options{}", "options specific to operator (see ampsci -o 'operator')"},
     {"rpa",
      "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram [true]"},
     {"omega",
      "Text or number. Frequency RPA is solved at. Put 'each' to solve at "
      "correct frequency for each transition. [0.0]"},
     {"delta", "Energy step for the numerical derivative of Sigma [1.0e-4]"},
     {"printBoth", "print <a|h|b> and <b|h|a> [false]"},
     {"use_basis", "If true, will use basis states for valence states [false]"},
     {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
     {"off-diagonal",
      "Calculate off-diagonal matrix elements (if non-zero) [true]"},
     {"what",
      "What to calculate? Options are: Reduced (reduced matric elements), "
      "Stetched (stretched states, with j=m= [j=min(ja,jb) for "
      "off-diagonal]), "
      "or HFConstant for (hyperfine A,B,etc. constants). Default is Reduced, "
      "except for hyperfine operator, for which it is HFConstant"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timer("normalisation");

  const auto delta = input.get("delta", 1.0e-4);

  if (!wf.Sigma()) {
    std::cout << "Correlation potential required to run this module!\n";
    return;
  }

  // Correlation potentials at e_v +/- delta, for the numerical derivative
  std::cout << "Setup Correlation potentials to calculate derivative:\n";
  const auto &Sigma0 = *wf.Sigma();

  auto Sigma1 = *wf.Sigma();
  Sigma1.clear();
  auto Sigma2 = Sigma1;

  std::cout << "\nCalculate correlation potentials:\n";
  for (const auto &v : wf.hf_valence()) {
    Sigma1.formSigma(v.kappa(), v.en() + delta, v.n(), &v);
    Sigma2.formSigma(v.kappa(), v.en() - delta, v.n(), &v);
    const auto lambda = Sigma0.getLambda(v.kappa(), v.n());
    std::cout << "Lambda = " << lambda << "\n\n";
  }

  //-------------------------------------------

  const auto oper = input.get<std::string>("operator", "");
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  const auto h = DiracOperator::generate(oper, h_options, wf);

  // treat hyperfine operator differently: constants instead of RME
  const bool is_hyperfine =
    qip::ci_compare(oper, "hfs") || qip::ci_compare(oper, "MLVP");

  // Determine "what" to calculate:
  const auto what_str =
    input.get<std::string>("what", is_hyperfine ? "HFConstant" : "Reduced");
  const auto matel_type = DiracOperator::parse_MatrixElementType(what_str);

  if (matel_type == DiracOperator::MatrixElementType::Reduced) {
    std::cout << "Reduced matrix elements\n";
  } else if (matel_type == DiracOperator::MatrixElementType::HFConstant) {
    const auto EM = h->rank() % 2 == 0 ? "E" : "M";
    std::cout << "Hyperfine constants " << EM << h->rank() << "\n";
  } else if (matel_type == DiracOperator::MatrixElementType::Stretched) {
    std::cout << "Stretched states with m=J [J=min(j_a,j_b) for off-diag]\n";
  } else {
    fmt2::warning();
    std::cout << " - Unkown matrix element type?\n";
  }
  std::cout << "Units: " << h->units() << "\n";

  const auto use_basis =
    wf.basis().empty() ? false : input.get("use_basis", false);
  if (use_basis) {
    std::cout << "Using basis (instead of valence) for matrix elements\n";
  }

  // RPA:
  auto rpa_method_str = input.get("rpa", std::string("true"));
  if (wf.core().empty())
    rpa_method_str = "false";
  auto rpa = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                     wf.basis(), wf.identity());

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = qip::ci_compare(str_om, "each");
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  // The matrix elements are between HF valence states (the norm correction
  // is for the missing Brueckner normalisation), optionally basis states
  std::vector<DiracSpinor> t_orbs;
  if (use_basis) {
    for (const auto &v : wf.hf_valence()) {
      const auto t = std::find(wf.basis().cbegin(), wf.basis().cend(), v);
      if (t != wf.basis().cend()) {
        t_orbs.push_back(*t);
      }
    }
  }
  const auto &orbs = use_basis ? t_orbs : wf.hf_valence();

  // Matrix elements: the operator follows each transition; the RPA is
  // solved here unless it too follows each transition
  using Amplitudes::Frequency;
  Amplitudes::MEoptions options{Frequency::transition, eachFreqQ ?
                                                         Frequency::transition :
                                                         Frequency::fixed};
  if (!eachFreqQ && rpa) {
    rpa->solve_core(omega, options.rpa_iterations, true);
  }
  options.diagonal = input.get("diagonal", true);
  options.off_diagonal = input.get("off-diagonal", true);
  options.calculate_both = input.get("printBoth", false);
  options.type = matel_type;

  const auto mes =
    Amplitudes::matrix_elements(orbs, h.get(), rpa.get(), nullptr, options);

  // Print, with the normalisation correction for each pair:
  // Norm = (1/2) (t_ab + dV_ab) (dSigma_a/de + dSigma_b/de)
  const auto find_orb = [&orbs](const std::string &sym) {
    const auto it = std::find_if(orbs.cbegin(), orbs.cend(), [&sym](auto &o) {
      return o.shortSymbol() == sym;
    });
    return it == orbs.cend() ? nullptr : &*it;
  };

  std::cout << "\na    b     t_ab         dv_ab        dS_a         dS_b    "
               "     Norm\n";
  for (const auto &m : mes) {
    const auto *Fa = find_orb(m.a);
    const auto *Fb = find_orb(m.b);
    if (Fa == nullptr || Fb == nullptr)
      continue;
    const auto dSa = Amplitudes::dSigma_dE(*Fa, Sigma0, Sigma1, Sigma2, delta);
    const auto dSb = Amplitudes::dSigma_dE(*Fb, Sigma0, Sigma1, Sigma2, delta);
    fmt::print("{:4s} {:4s} {:+.5e} {:+.5e} {:+.5e} {:+.5e} {:+.5e}\n", m.a,
               m.b, m.factor * m.t0, m.factor * m.dv, dSa, dSb,
               0.5 * m.value() * (dSa + dSb));
  }
  std::cout << "\n";
}

} // namespace Module
