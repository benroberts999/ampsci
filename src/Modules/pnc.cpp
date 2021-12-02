#include "pnc.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>

namespace Module {

using namespace Pnc;

//******************************************************************************
void calculatePNC(const IO::InputBlock &input, const Wavefunction &wf) {
  const std::string ThisModule = "Module::PNC";

  input.checkBlock_old({"t", "c", "transition", "nmain", "rpa", "omega",
                        "E1_rpa_it", "pnc_rpa_it"});

  // input: nuc parameters for rho:
  const auto c_dflt =
      Nuclear::c_hdr_formula_rrms_t(Nuclear::find_rrms(wf.Znuc(), wf.Anuc()));
  const auto t = input.get("t", Nuclear::default_t);
  const auto c = input.get("c", c_dflt);

  // input: rpa? and which n to consider for 'main'
  const auto main_n = input.get("nmain", wf.maxCore_n() + 4);
  const auto rpaQ = input.get("rpa", true);

  // input: transition
  const auto Fb_vec = input.get<std::vector<std::string>>("transition", {});
  const auto Fb_ok = (Fb_vec.size() >= 2);
  const auto [na, ka] =
      Fb_ok ? AtomData::parse_symbol(Fb_vec[0]) : std::pair{0, 0};
  const auto [nb, kb] =
      Fb_ok ? AtomData::parse_symbol(Fb_vec[1]) : std::pair{0, 0};

  const auto pA = wf.getState(na, ka);
  const auto pB = wf.getState(nb, kb);
  if (!pA || !pB) {
    std::cerr << "\nFAIL in Module:PNC\n"
              << "Couldn't find requested state: ("
              << input.get<std::string>("transition", "")
              << "). Is it in valence list?\n";
    return;
  }

  const auto &Fa = *pA;
  const auto &Fb = *pB;

  std::cout << "\n********************************************** \n";
  std::cout << "E_pnc (B->A): " << wf.atom() << ":   A = " << Fa.symbol()
            << " ,  B = " << Fb.symbol() << "\n"
            << "z-component, z=min(ja,jb). units: i(-Qw/N)10^-11.\n";
  std::cout << "w/ c=" << c << ", t=" << t << "\n";

  // Generate operators:
  const auto N_nuc = wf.Anuc() - wf.Znuc();
  DiracOperator::PNCnsi hpnc(c, t, *wf.rgrid, -N_nuc, "i(Qw/-N)*e-11");
  DiracOperator::E1 he1(*wf.rgrid);

  // Find core/valence energy: allows distingush core/valence states
  const auto ec_max =
      std::max_element(cbegin(wf.core), cend(wf.core), DiracSpinor::comp_en)
          ->en();
  const auto ev_min = std::min_element(cbegin(wf.valence), cend(wf.valence),
                                       DiracSpinor::comp_en)
                          ->en();
  const auto en_core = 0.5 * (ev_min + ec_max);

  // TDHF (nb: need object even if not doing RPA)
  auto dVE1 = ExternalField::TDHF(&he1, wf.getHF());
  auto dVpnc = ExternalField::TDHF(&hpnc, wf.getHF());
  if (rpaQ) {
    const auto omega_dflt = std::abs(Fa.en() - Fb.en());
    const auto omega = input.get("omega", omega_dflt);
    auto E1_it = input.get("E1_rpa_it", 99);
    auto pnc_it = input.get("pnc_rpa_it", 99);
    dVE1.solve_core(omega, E1_it);
    dVpnc.solve_core(0.0, pnc_it);
  }

  // SOS, use HF
  std::cout << "\nUsing HF states:\n";
  auto hf_basis = wf.core;
  hf_basis.insert(end(hf_basis), cbegin(wf.valence), cend(wf.valence));
  pnc_sos(Fa, Fb, &hpnc, &dVpnc, &he1, &dVE1, hf_basis, main_n, en_core, true);

  // SOS, use spectrum
  std::cout << "\nUsing spectrum:\n";
  const auto [sos1, sos2] = pnc_sos(Fa, Fb, &hpnc, &dVpnc, &he1, &dVE1,
                                    wf.spectrum, main_n, en_core, true);
  // Can swap order of the E1/PNCthe operators: ("trivially" the same for sos)
  // pnc_sos(Fa, Fb, &he1, &dVE1, &hpnc, &dVpnc, wf.spectrum, main_n, en_core,
  //         true);

  // Solving equations method (E1 amplitude, PNC perturbed)
  std::cout << "\n";
  const auto [se_d1, se_d2] =
      pnc_tdhf(Fa, Fb, &hpnc, &dVpnc, &he1, &dVE1, wf.getSigma(), wf.spectrum,
               main_n, en_core, true);

  // Solving equations method (PNC amplitude, E1 perturbed)
  std::cout << "\n";
  const auto [se_h1, se_h2] =
      pnc_tdhf(Fa, Fb, &he1, &dVE1, &hpnc, &dVpnc, wf.getSigma(), wf.spectrum,
               main_n, en_core, true);

  // Calculate relative difference (numerical accuracy)
  const auto eps_sos1 = std::abs((sos1 - se_d1) / (sos1 + se_d1));
  const auto eps_sos2 = std::abs((sos2 - se_d2) / (sos2 + se_d2));
  const auto eps_sos =
      std::abs((sos1 + sos2 - se_d1 - se_d2) / (sos1 + sos2 + se_d1 + se_d2));

  // note: for SEs, the terms appear swapped.. (just order changes, sign same)
  const auto eps_se1 = std::abs((se_h1 - se_d2) / (se_h1 + se_d2));
  const auto eps_se2 = std::abs((se_h1 - se_d2) / (se_h1 + se_d2));
  const auto eps_se = std::abs((se_h1 + se_h2 - se_d1 - se_d2) /
                               (se_h1 + se_h2 + se_d1 + se_d2));
  printf("\neps(sos/SEs): %.0e, %.0e  %.1e\n", eps_sos1, eps_sos2, eps_sos);
  printf("eps(SEs)    : %.0e, %.0e  %.1e\n", eps_se1, eps_se2, eps_se);
}

//******************************************************************************
//******************************************************************************
namespace Pnc {

//******************************************************************************
std::pair<double, double> pnc_sos(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                  const DiracOperator::TensorOperator *hpnc,
                                  const ExternalField::TDHF *dVpnc,
                                  const DiracOperator::TensorOperator *he1,
                                  const ExternalField::TDHF *dVE1,
                                  const std::vector<DiracSpinor> &spectrum,
                                  int main_n, double en_core, bool print) {

  if (print) {
    std::cout << "Sum-over-states method.\n";
    std::cout << "<" << Fa.shortSymbol() << "|" << he1->name() << "|n><n|"
              << hpnc->name() << "|" << Fb.shortSymbol() << ">/dE + <"
              << Fa.shortSymbol() << "|" << hpnc->name() << "|n><n|"
              << he1->name() << "|" << Fb.shortSymbol() << ">/dE\n";
  }

  const bool conj = Fa.en() < Fb.en() ? true : false;

  // Print each core+tail term (for testing)
  const auto print_all = false;

  // For angular factors (RME -> z-comp, z=min(ja,jb))
  const auto tja = Fa.twoj();
  const auto tjb = Fb.twoj();
  // z-comp defined as min(a,b)
  const auto twom = std::min(tja, tjb);

  // This omega is so we can swap the E1 and PNC operators!
  const auto w_SE = hpnc->imaginaryQ() ? 0.0 : (Fa.en() - Fb.en());

  double pnc_1 = 0.0, core_1 = 0.0, main_1 = 0.0;
  double pnc_2 = 0.0, core_2 = 0.0, main_2 = 0.0;
  for (auto &np : spectrum) {
    if (np == Fb || np == Fa)
      continue;
    if (hpnc->isZero(np.k, Fa.k) && hpnc->isZero(np.k, Fb.k))
      continue;
    if (he1->isZero(np.k, Fa.k) && he1->isZero(np.k, Fb.k))
      continue;
    const auto coreQ = np.en() < en_core;
    const auto mainQ = !coreQ && np.n <= main_n;

    // nb: need 'conj' here, since w = |w|, and want work for a->b and b->a ?
    const auto dAp = he1->reducedME(Fa, np) + dVE1->dV(Fa, np, conj);
    const auto hpB = hpnc->reducedME(np, Fb) + dVpnc->dV(np, Fb, conj);
    const auto hAp = hpnc->reducedME(Fa, np) + dVpnc->dV(Fa, np, conj);
    const auto dpB = he1->reducedME(np, Fb) + dVE1->dV(np, Fb, conj);

    // Angular factors:
    // nb: these are always the same for pnc... but in general may not be?
    // so, can use for other operators
    const auto tjn = np.twoj();
    const auto c10 = he1->rme3js(tja, tjn, twom) * hpnc->rme3js(tjn, tjb, twom);
    const auto c01 = hpnc->rme3js(tja, tjn, twom) * he1->rme3js(tjn, tjb, twom);

    const double pnc1 = c10 * dAp * hpB / (Fb.en() - np.en() + w_SE);
    const double pnc2 = c01 * hAp * dpB / (Fa.en() - np.en() - w_SE);
    if (np.n <= main_n && print_all)
      printf("%7s, pnc= %12.5e + %12.5e = %12.5e\n", np.symbol().c_str(), pnc1,
             pnc2, pnc1 + pnc2);

    pnc_1 += pnc1;
    pnc_2 += pnc2;
    if (coreQ) {
      core_1 += pnc1;
      core_2 += pnc2;
    } else if (mainQ) {
      main_1 += pnc1;
      main_2 += pnc2;
    }
  }

  const auto printer = [](const char *str, double p1, double p2) {
    if (p1 != 0.0 || p1 != 0.0)
      printf("%5s: %11.6f + %11.6f = %11.6f\n", str, p1, p2, p1 + p2);
  };
  if (print) {
    printer("Core ", core_1, core_2);
    printer("Main ", main_1, main_2);
    printer("Tail ", pnc_1 - main_1 - core_1, pnc_2 - main_2 - core_2);
    printer("Total", pnc_1, pnc_2);
  }

  return {pnc_1, pnc_2};
}

//******************************************************************************
DiracSpinor orthog_to_core(DiracSpinor dF,
                           const std::vector<DiracSpinor> &in_orbs,
                           double en_core) {
  for (const auto &Fc : in_orbs) {
    const auto coreQ = Fc.en() < en_core;
    if (dF.k == Fc.k && coreQ)
      dF -= (dF * Fc) * Fc;
  }
  return dF;
}

DiracSpinor orthog_to_coremain(DiracSpinor dF,
                               const std::vector<DiracSpinor> &in_orbs,
                               double en_core, int n_main) {
  for (const auto &Fc : in_orbs) {
    const auto coreQ = Fc.en() < en_core;
    const auto mainQ = !coreQ && Fc.n <= n_main;
    if (dF.k == Fc.k && (coreQ || mainQ))
      dF -= (dF * Fc) * Fc;
  }
  return dF;
}

//******************************************************************************
std::pair<double, double> pnc_tdhf(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                   const DiracOperator::TensorOperator *hpnc,
                                   const ExternalField::TDHF *dVpnc,
                                   const DiracOperator::TensorOperator *he1,
                                   const ExternalField::TDHF *dVE1,
                                   const MBPT::CorrelationPotential *Sigma,
                                   const std::vector<DiracSpinor> &spectrum,
                                   int main_n, double en_core, bool print) {
  // Note: calling with e1 and pnc swapped is valid! (and a good consistancy
  // check!)

  if (print) {
    std::cout << "TDHF (Solving-equations) method. hw: " << hpnc->name()
              << "\n";
  }

  // For angular factors (RME -> z-comp, z=min(ja,jb))
  const auto tja = Fa.twoj();
  const auto tjb = Fb.twoj();
  // z-comp defined as min(a,b)
  const auto twom = std::min(tja, tjb);

  const bool conj = Fa.en() < Fb.en() ? true : false;

  // Allow swapping the 'pnc' and 'E1' operators:
  // note: MUST be this, not RPA w
  const auto w_SE = hpnc->imaginaryQ() ? 0.0 : (Fa.en() - Fb.en());

  //
  std::cout << "<" << Fa.shortSymbol() << "|" << he1->name() << "|d"
            << Fb.shortSymbol() << "> + "
            << "<d" << Fa.shortSymbol() << "|" << he1->name() << "|"
            << Fb.shortSymbol() << "> ; w/ h = " << hpnc->name() << "\n";
  const auto XB = dVpnc->solve_dPsis(Fb, w_SE, ExternalField::dPsiType::X,
                                     Sigma, ExternalField::StateType::ket);
  const auto YA = dVpnc->solve_dPsis(Fa, w_SE, ExternalField::dPsiType::Y,
                                     Sigma, ExternalField::StateType::bra);

  // get total PNC amplitude:
  double pnc1 = 0.0, pnc2 = 0.0;
  for (auto &xb : XB) {
    const auto c10 =
        he1->rme3js(tja, xb.twoj(), twom) * hpnc->rme3js(xb.twoj(), tjb, twom);
    pnc1 += c10 * (he1->reducedME(Fa, xb) + dVE1->dV(Fa, xb, conj));
  }
  for (auto &ya : YA) {
    const auto c01 =
        hpnc->rme3js(tja, ya.twoj(), twom) * he1->rme3js(ya.twoj(), tjb, twom);
    pnc2 += c01 * (he1->reducedME(ya, Fb) + dVE1->dV(ya, Fb, conj));
  }

  // Find main +tail, by forcing "mixed states" to be orthog to core/main:
  double pnc1_t = 0.0, pnc2_t = 0.0;
  double pnc1_m = 0.0, pnc2_m = 0.0;
  for (const auto &xb : XB) {
    // check if zero, skip
    const auto c10 =
        he1->rme3js(tja, xb.twoj(), twom) * hpnc->rme3js(xb.twoj(), tjb, twom);
    // Force orthogonal to core (leaving main+tail):
    const auto x_mt = orthog_to_core(xb, spectrum, en_core);
    const auto mt = c10 * (he1->reducedME(Fa, x_mt) + dVE1->dV(Fa, x_mt, conj));
    // Force orthogonal to main (leaving just tail):
    const auto x_t = orthog_to_coremain(xb, spectrum, en_core, main_n);
    const auto tail = c10 * (he1->reducedME(Fa, x_t) + dVE1->dV(Fa, x_t, conj));
    pnc1_t += tail;
    pnc1_m += (mt - tail);
  }
  for (const auto &ya : YA) {
    const auto c01 =
        hpnc->rme3js(tja, ya.twoj(), twom) * he1->rme3js(ya.twoj(), tjb, twom);
    // Force orthogonal to core (leaving main+tail):
    const auto y_mt = orthog_to_core(ya, spectrum, en_core);
    const auto mt = c01 * (he1->reducedME(y_mt, Fb) + dVE1->dV(y_mt, Fb, conj));
    // Force orthogonal to main (leaving just tail):
    const auto y_t = orthog_to_coremain(ya, spectrum, en_core, main_n);
    const auto tail = c01 * (he1->reducedME(y_t, Fb) + dVE1->dV(y_t, Fb, conj));
    pnc2_t += tail;
    pnc2_m += (mt - tail);
  }

  // print results:
  const auto printer = [](const char *str, double p1, double p2) {
    printf("%5s: %11.6f + %11.6f = %11.6f\n", str, p1, p2, p1 + p2);
  };
  const auto pnc1_c = pnc1 - pnc1_m - pnc1_t;
  const auto pnc2_c = pnc2 - pnc2_m - pnc2_t;
  if (print) {
    printer("Core ", pnc1_c, pnc2_c);
    printer("Main ", pnc1_m, pnc2_m);
    printer("Tail ", pnc1_t, pnc2_t);
    printer("Total", pnc1, pnc2);
  }

  return {pnc1, pnc2};
}
} // namespace Pnc

} // namespace Module
