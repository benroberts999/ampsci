#include "lifetimes.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

//! E1 partial rate: d = <f||E1||i>, gi = [Ji] = 2Ji+1, w = Ef - Ei
double gamma_E1(double d, double w, double gi) {
  const auto aw = std::abs(PhysConst::alpha * w);
  const auto C = (4.0 / 3.0) * qip::pow<3>(aw);
  return C * qip::pow<2>(d) / gi;
}

//! M1 partial rate: m = <f||M1||i>,  gi = [Ji] = 2Ji+1, w = Ef - Ei
double gamma_M1(double m, double w, double gi) {
  const auto aw = std::abs(PhysConst::alpha * w);
  const auto u2 = PhysConst::muB_CGS * PhysConst::muB_CGS;
  const auto C = (4.0 / 3.0) * qip::pow<3>(aw) * u2;
  return C * qip::pow<2>(m) / gi;
}

//! E2 partial rate: q = <f||E2||i>,  gi = [Ji] = 2Ji+1, w = Ef - Ei
double gamma_E2(double q, double w, double gi) {
  const auto aw = std::abs(PhysConst::alpha * w);
  const auto C = (1.0 / 15.0) * qip::pow<5>(aw);
  return C * qip::pow<2>(q) / gi;
}

namespace Module {

void lifetimes(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nLifetimes (v2):\n";
  std::cout << "Note: Uses _valence_ states - so, must ensure all lower states "
               "have been included in the valence list for accurate results.\n";

  input.check(
      {{"E1", "Include E1 transitions? [true]"},
       {"E2", "Include E2 transitions? [false]"},
       {"M1", "Include M1 transitions? [false]"},
       {"rpa", "Include RPA? [true]"},
       {"SRN", "Include SR+Norm correction? [false]"},
       {"Qk_file", "Qk filename (for SR+N). If blank, will calculate "
                   "from scratch, ortherwise will read/write to Qkfile. "
                   "Note: QkFile assumes spline-legs! [blank]"},
       {"n_minmax", "List: minimum core n, maximum excited n to "
                    "include for SR+N [2,30]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto do_E1 = input.get("E1", true);
  const auto do_E2 = input.get("E2", false);
  const auto do_M1 = input.get("M1", false);
  const auto rpaQ = input.get("rpa", true);
  const auto do_SRN = input.get("SRN", false);
  const auto n_min_max = input.get("n_minmax", std::vector<int>{2, 30});
  const auto Qkfile = input.get("Qk_file", std::string{});
  assert(n_min_max.size() >= 2);

  if (rpaQ && wf.basis().empty() && (do_E2 || do_M1)) {
    std::cout << "\nWARNING: No basis!\nWe use diagram RPA method for E2 and "
                 "M1 - so we require a basis!\n\n";
  }

  // Construct operators:
  DiracOperator::E1 e1(wf.grid());
  DiracOperator::Ek e2(wf.grid(), 2);
  DiracOperator::M1 m1(wf.grid(), wf.alpha());

  // Include RPA (optionally)
  using namespace ExternalField;
  auto dVe1 = (do_E1 && rpaQ) ? std::make_unique<TDHF>(&e1, wf.vHF()) : nullptr;
  // nb: we use DiagramRPA method for E2 and M1
  auto dVe2 = (do_E2 && rpaQ) ? std::make_unique<DiagramRPA>(
                                    &e2, wf.basis(), wf.vHF(), wf.identity()) :
                                nullptr;
  auto dVm1 = (do_M1 && rpaQ) ? std::make_unique<DiagramRPA>(
                                    &m1, wf.basis(), wf.vHF(), wf.identity()) :
                                nullptr;

  // Optionally include SR+N
  if (do_SRN)
    std::cout << "\nIncluding SR+Norm:\n";
  auto SRN = do_SRN ? std::make_unique<MBPT::StructureRad>(
                          wf.basis(), wf.FermiLevel(),
                          std::pair{n_min_max.at(0), n_min_max.at(1)}, Qkfile) :
                      nullptr;

  // Find ground state (so we don't print it)
  const auto ground_state = std::min_element(
      wf.valence().cbegin(), wf.valence().cend(), DiracSpinor::comp_en);
  if (ground_state == wf.valence().cend())
    return;

  // 1. Calculate required matrix elements
  std::vector<MEdata> e1s, e2s, m1s;
  Coulomb::meTable<std::pair<double, double>> SRNe1s, SRNe2s, SRNm1s;
  if (do_E1) {
    e1s = calcMatrixElements(wf.valence(), &e1, dVe1.get(), 0.0, true);
    if (SRN)
      SRNe1s = SRN->srn_table(&e1, wf.valence());
    std::cout << "\nE1 reduced matrix elements:\n";
    std::cout << MEdata::title(rpaQ) << " ";
    if (SRN)
      std::cout << "         SR+N";
    std::cout << "\n";
    for (const auto &each : e1s) {
      std::cout << each;
      auto srn = SRNe1s.get(each.a, each.b);
      if (srn)
        printf(" %11.4e", srn->first);
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  if (do_E2) {
    e2s = calcMatrixElements(wf.valence(), &e2, dVe2.get(), 0.0, true);
    if (SRN)
      SRNe2s = SRN->srn_table(&e2, wf.valence());
    std::cout << "\nE2 reduced matrix elements:\n";
    std::cout << MEdata::title(rpaQ) << " ";
    if (SRN)
      std::cout << "         SR+N";
    std::cout << "\n";
    for (const auto &each : e2s) {
      std::cout << each;
      auto srn = SRNe2s.get(each.a, each.b);
      if (srn)
        printf(" %11.4e", srn->first);
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  if (do_M1) {
    m1s = calcMatrixElements(wf.valence(), &m1, dVm1.get(), 0.0, true);
    if (SRN)
      SRNm1s = SRN->srn_table(&m1, wf.valence());
    std::cout << "\nM1 reduced matrix elements:\n";
    std::cout << MEdata::title(rpaQ) << " ";
    if (SRN)
      std::cout << "         SR+N";
    std::cout << "\n";
    for (const auto &each : m1s) {
      std::cout << each;
      auto srn = SRNm1s.get(each.a, each.b);
      if (srn)
        printf(" %11.4e", srn->first);
      std::cout << "\n";
    }
    std::cout << "\n";
  }

  auto finder = [](const auto &ta, const auto &tb) {
    const auto a = ta.shortSymbol();
    const auto b = tb.shortSymbol();
    return [a, b](auto &el) {
      return (el.a == a && el.b == b) || (el.a == b && el.b == a);
    };
  };

  // Calculate widths (and lifetimes) + print to screen
  std::cout << "\nPartial widths and lifetimes:\n";
  std::vector<std::pair<std::string, double>> results;
  for (const auto &i : wf.valence()) {
    if (i == *ground_state)
      continue;

    std::cout << "\n" << i << ":\n";
    std::cout << "   i    f   |𝜔_if| ";
    if (do_E1)
      std::cout << "       𝛾(E1)";
    if (do_E2)
      std::cout << "       𝛾(E2)";
    if (do_M1)
      std::cout << "       𝛾(M1)";
    std::cout << "\n";

    const auto gi = i.twojp1();
    double G_E1 = 0.0;
    double G_E2 = 0.0;
    double G_M1 = 0.0;
    for (const auto &f : wf.valence()) {

      // only consider transitions to lower energy states:
      if (f.en() >= i.en())
        continue;
      if (e1.isZero(i, f) && e2.isZero(i, f) && m1.isZero(i, f))
        continue;

      printf(" %4s %4s %7.5f ", i.shortSymbol().c_str(),
             f.shortSymbol().c_str(), std::abs(i.en() - f.en()));

      const auto e1_data = std::find_if(e1s.begin(), e1s.end(), finder(i, f));
      const auto e2_data = std::find_if(e2s.begin(), e2s.end(), finder(i, f));
      const auto m1_data = std::find_if(m1s.begin(), m1s.end(), finder(i, f));

      //E1
      if (e1_data != e1s.end()) {
        const auto &[a, b, w, d, dv] = *e1_data;
        const auto t_srn = SRNe1s.get(a, b);
        const auto srn = t_srn ? t_srn->first : 0.0;
        const auto ge1_vf = gamma_E1(d + dv + srn, w, gi);
        printf("  %.4e", ge1_vf);
        G_E1 += ge1_vf;
      } else if (do_E1) {
        printf("  %10i", 0);
      }
      if (e2_data != e2s.end()) {
        const auto &[a, b, w, d, dv] = *e2_data;
        const auto t_srn = SRNe2s.get(a, b);
        const auto srn = t_srn ? t_srn->first : 0.0;
        const auto ge1_vf = gamma_E2(d + dv + srn, w, gi);
        printf("  %.4e", ge1_vf);
        G_E2 += ge1_vf;
      } else if (do_E2) {
        printf("  %10i", 0);
      }
      if (m1_data != m1s.end()) {
        const auto &[a, b, w, d, dv] = *m1_data;
        const auto t_srn = SRNm1s.get(a, b);
        const auto srn = t_srn ? t_srn->first : 0.0;
        const auto ge1_vf = gamma_M1(d + dv + srn, w, gi);
        printf("  %.4e", ge1_vf);
        G_M1 += ge1_vf;
      } else if (do_M1) {
        printf("  %10i", 0);
      }
      std::cout << "\n";
    }
    std::cout << "Total              ";
    if (do_E1)
      printf("  %.4e", G_E1);
    if (do_E2)
      printf("  %.4e", G_E2);
    if (do_M1)
      printf("  %.4e", G_M1);
    std::cout << "\n";

    // std::cout << "\n";
    printf("𝚪(%3s) = %.6e au\n", i.shortSymbol().c_str(), G_E1 + G_E2 + G_M1);
    const auto tau_s = PhysConst::time_s / (G_E1 + G_E2 + G_M1);
    printf("𝜏(%3s) = %.6e s\n", i.shortSymbol().c_str(), tau_s);
    results.push_back({i.shortSymbol(), tau_s});
  }

  // summary of results:
  // Convert from s to ns or ms
  std::cout << "\nSummary:\n";
  for (const auto &[state, tau_s] : results) {

    const auto tau_s2 = tau_s; // Structured binding cannot be captured?
    const auto [tau_u, u] = [tau_s2]() {
      if (tau_s2 > 1.0)
        return std::pair{tau_s2, "s"};
      if (tau_s2 > 1.0e-6)
        return std::pair{tau_s2 * 1.0e3, "ms"};
      return std::pair{tau_s2 * 1.0e9, "ns"};
    }();

    printf("𝜏(%3s) = %.5e s = %.3f %s\n", state.c_str(), tau_s, tau_u, u);
  }
}
} // namespace Module
