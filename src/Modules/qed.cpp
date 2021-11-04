#include "Modules/qed.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/matrixElements.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include <numeric>
#include <string>
#include <vector>

namespace Module {

//****************************************************************************
// Used for finding A and b for effective vertex QED operator
void QED(const IO::InputBlock &input, const Wavefunction &wf) {

  std::cout << "\nQED Module:\n";

  // Check input options for spelling mistakes etc.:
  input.check(
      {{"A_vertex", "A vtx factor; blank means dflt"},
       {"b_vertex", "B vtx factor; =1 by default"},
       {"rrms", "double; effective rrms used in radiative potential"},
       {"out_file", "Results appended to this file (if given)"},
       {"matrixElements", "sub-block; same as Module::matrixElements"},
       {"coreQED", "bool; Include QED into core? Or just valence"},
       {"scale_l",
        "list; Scale factors for each l; e.g., 1,0 means s, but no p,d. "
        "default = 1, meaning include all"},
       {"vertex", "bool; Calculate vertex corrections?"}});

  std::cout << "\nDO NOT include QED into wavefunction in AMPSCI\n";

  // filename for output file
  const auto fname = input.get<std::string>("out_file", "");

  std::ofstream of_en;

  auto first_line_en =
      IO::FRW::file_exists(fname + ".energy")
          ? ""
          : "#          E0            dE(VP)        dE(SE)        dE(tot)\n";

  if (fname != "")
    of_en.open(fname + ".energy", std::ios_base::app);
  of_en << first_line_en;
  // append instead of overwrite

  // allow rrms to be different from charge distribution of wavefunction
  const auto r_rmsfm0 = input.get("rrms", wf.get_rrms());
  const auto scale_rN = r_rmsfm0 / wf.get_rrms();

  // New Wavefunctions, which include the QED
  Wavefunction wf_VP(wf.rgrid->params(), wf.get_nuclearParameters());
  Wavefunction wf_SE(wf.rgrid->params(), wf.get_nuclearParameters());

  // Diference scale for each l (mainly for core/relaxation)
  const auto x_spd = input.get("scale_l", std::vector{1.0});

  const auto coreQED = input.get("coreQED", true);
  if (coreQED) {
    wf_VP.radiativePotential({1.0, 0.0, 0.0, 0.0, 0.0}, 15.0, scale_rN, x_spd);
    wf_SE.radiativePotential({0.0, 1.0, 1.0, 1.0, 0.0}, 15.0, scale_rN, x_spd);
  }

  // Do Hartree-Fock on new wavefunctions
  const auto wfhf = wf.getHF();
  if (!wf.core.empty()) {
    std::cout << "HF, including VP, SE: " << HF::parseMethod(wfhf->method())
              << "\n";
  }
  wf_VP.solve_core(HF::parseMethod(wfhf->method()), wfhf->x_Breit(),
                   wf.coreConfiguration());
  wf_SE.solve_core(HF::parseMethod(wfhf->method()), wfhf->x_Breit(),
                   wf.coreConfiguration());

  if (!coreQED) {
    std::cout << "no core QED; add QED _after_ core HartreeFock\n";
    wf_VP.radiativePotential({1.0, 0.0, 0.0, 0.0, 0.0}, 15.0, scale_rN, x_spd);
    wf_SE.radiativePotential({0.0, 1.0, 1.0, 1.0, 0.0}, 15.0, scale_rN, x_spd);
  }

  // Copy correlation potential to new wfs
  if (wf.getSigma() != nullptr) {
    std::cout << "\nIncluding correlations.\n"
                 "Copying correlation potential across to new wavefunctions\n"
                 "Note: SAME correlation potential; i.e., QED not included in "
                 "correlation potential.\n";
    wf_VP.copySigma(wf.getSigma());
    wf_SE.copySigma(wf.getSigma());
  }

  // Solve valence states for new wfs
  if (wf.core.empty()) {
    // H-like:
    wf_VP.localValence(DiracSpinor::state_config(wf.valence));
    wf_SE.localValence(DiracSpinor::state_config(wf.valence));
  } else if (wfhf->method() == HF::Method::KohnSham) {
    // Kohn Sham input valence list has different format
    std::string list{};
    for (const auto &Fv : wf.valence) {
      if (Fv.k < 0)
        list += Fv.shortSymbol().substr(0, Fv.shortSymbol().length() - 1);
    }
    std::cout << list << "\n";
    wf_VP.localValence(list, true);
    wf_SE.localValence(list, true);
  } else {
    // Hartree-Fock:
    wf_VP.solve_valence(DiracSpinor::state_config(wf.valence));
    wf_SE.solve_valence(DiracSpinor::state_config(wf.valence));
  }
  // Solve Brueckner orbitals
  if (wf.getSigma() != nullptr) {
    wf_VP.hartreeFockBrueckner();
    wf_SE.hartreeFockBrueckner();
  }

  // print the new valence energies to screen:
  std::cout << "\nEnergies, including VP:\n";
  wf_VP.printValence();
  std::cout << "Energies, including SE:\n";
  wf_SE.printValence();

  // Copy over basis (if exists) - only used for RPA (basis/diagram)
  if (!wf.basis.empty()) {
    std::cout << "\nCopying basis across; Note: basis does NOT include "
                 "QED\n(Basis may be used for RPA, otherwise not used)\n";
    wf_VP.basis = wf.basis;
    wf_SE.basis = wf.basis;
  }

  // Operator (for first-order shift) - l-dependent
  std::cout
      << "\nFirst-order QED contribution to valence energies; atomic units\n";
  std::cout << "       dE(VP)        dE(SE)        dE(tot)\n";
  for (const auto &Fv : wf.valence) {
    const auto h_SE =
        DiracOperator::Hrad(wf_SE.qed->Vel(Fv.l()), wf_SE.qed->Hmag(Fv.l()));
    const auto h_VP =
        DiracOperator::Hrad(wf_VP.qed->Vel(Fv.l()), wf_VP.qed->Hmag(Fv.l()));

    const auto dese = h_SE.radialIntegral(Fv, Fv);
    const auto devp = h_VP.radialIntegral(Fv, Fv);
    printf("%4s  %12.5e  %12.5e  %12.5e\n", Fv.shortSymbol().c_str(), devp,
           dese, devp + dese);
  }

  // Output: energy results
  std::cout << "\nQED contribution to valence energies; atomic units (HF)\n";
  if (coreQED) {
    std::cout << "Includes core relaxation\n";
  }
  std::cout << "       E0            dE(VP)        dE(SE)        dE(tot)\n";
  for (const auto &Fv : wf.valence) {
    // Find corresponding state without QED: (can't assume in same order)
    const auto Fv_SE = wf_SE.getState(Fv.n, Fv.k);
    const auto Fv_VP = wf_VP.getState(Fv.n, Fv.k);
    if (!Fv_SE || !Fv_SE)
      continue;

    const auto de_vp = Fv_VP->en() - Fv.en();
    const auto de_se = Fv_SE->en() - Fv.en();
    const auto o = qip::fstring("%4s  %12.5e  %12.5e  %12.5e  %12.5e\n",
                                Fv.shortSymbol().c_str(), Fv.en(), de_vp, de_se,
                                de_vp + de_se);
    std::cout << o;
    of_en << wf.Znuc() << " " << o;
  }

  // QED to matrix elements (perturbed orbital part):
  const auto me_input = input.getBlock("matrixElements");
  if (me_input) {
    me_input->check(
        {{"operator", "e.g., E1, hfs"},
         {"options", "options specific to operator; blank by dflt"},
         {"rpa", "true(=TDHF), false, TDHF, basis, diagram"},
         {"omega", "freq. for RPA"},
         {"radialIntegral", "false by dflt (means red. ME)"},
         {"printBoth", "print <a|h|b> and <b|h|a> (dflt false)"},
         {"onlyDiagonal", "only <a|h|a> (dflt false)"}});

    const auto oper = me_input->get<std::string>("operator", "");
    // Get optional 'options' for operator
    auto h_options = IO::InputBlock(oper, {});
    const auto tmp_opt = me_input->getBlock("options");
    if (tmp_opt) {
      h_options = *tmp_opt;
    }
    const auto h = generateOperator(oper, h_options, wf, true);
    const bool diagonal_only = me_input->get("onlyDiagonal", false);

    const bool rpaQ = me_input->get("rpa", false);
    std::unique_ptr<ExternalField::CorePolarisation> rpa0{nullptr},
        rpa_vp{nullptr}, rpa_se{nullptr};
    if (rpaQ) {
      std::cout << "RPA: TDHF method, at zero frequency\n";
      rpa0 = std::make_unique<ExternalField::TDHF>(h.get(), wf.getHF());
      rpa_vp = std::make_unique<ExternalField::TDHF>(h.get(), wf_VP.getHF());
      rpa_se = std::make_unique<ExternalField::TDHF>(h.get(), wf_SE.getHF());
    }

    auto first_line_me_po =
        IO::FRW::file_exists(fname + ".me_po")
            ? ""
            : "#                ME(0)         d(VP)         "
              "d(SE)         d(tot)\n";

    std::ofstream of_me;
    if (fname != "")
      of_me.open(fname + ".me_po", std::ios_base::app);
    of_me << first_line_me_po;

    // std::cout << "\nQED correction to Matrix elements (PO)\n";
    // std::cout << "\nNo QED:";
    const auto me0 = ExternalField::calcMatrixElements(
        wf.valence, h.get(), rpa0.get(), 0.0, false, diagonal_only);

    // std::cout << "\nVacuum polarisation (PO):";
    const auto mevp = ExternalField::calcMatrixElements(
        wf_VP.valence, h.get(), rpa_vp.get(), 0.0, false, diagonal_only);

    // std::cout << "\nSelf-energy (PO):";
    const auto mese = ExternalField::calcMatrixElements(
        wf_SE.valence, h.get(), rpa_se.get(), 0.0, false, diagonal_only);

    std::cout << "\nQED contribution matrix elements (Perturbed Orbital)\n";
    std::cout
        << "             ME(0)         d(VP)         d(SE)         d(tot)\n";
    for (const auto &[a, b, x0, dv1, dv] : me0) {
      // find corresponding ME calculations in list:
      // (can't assume order will be the same, though usually will be)
      auto l = [a = a, b = b](auto e) { return (a == e.a && b == e.b); };
      // nb: a=a here due to clang error capturing bindings...
      const auto vp = std::find_if(begin(mevp), end(mevp), l);
      const auto se = std::find_if(begin(mese), end(mese), l);
      if (vp != mevp.end() && se != mese.end()) {
        const auto d_vp = (vp->hab + vp->dv) - (x0 + dv);
        const auto d_se = (se->hab + se->dv) - (x0 + dv);
        const auto o = qip::fstring(
            "%4s %4s:  %12.5e  %12.5e  %12.5e  %12.5e\n", a.c_str(), b.c_str(),
            x0 + dv, d_vp, d_se, d_vp + d_se);
        std::cout << o;
        of_me << wf.Znuc() << " " << o;
      }
    }

    // For vertexQED:
    if (input.get("vertex", false)) {

      std::cout << "\nQED contribution (Vertex) to radial integrals\n";
      if (input.get("rpa", false)) {
        std::cout << "Note: NO RPA in vertex (yet)!!\n";
      }

      std::cout << "Unit conversions:\n";
      const auto a2 = std::pow(PhysConst::alpha, 2);
      const auto z3 = std::pow(wf.Znuc(), 3);
      const auto gImmPoI = h->getc();
      // gImmPoI = mu * PhysConst::muN_CGS_MHz / IN (mu here in units of muN)
      const auto factor_xRad =
          (4.0 / 3) * a2 * z3 * gImmPoI / PhysConst::muB_CGS;
      const auto factor_eF = (PhysConst::alpha / M_PI) * factor_xRad; //?????

      std::cout << wf.atom() << " " << wf.nuclearParams() << "\n";
      // "muN * PhysConst::muN_CGS_MHz / IN"
      // "(4/3)Z^3a^2 (g/mp)"
      std::cout << "factor_xRad = A_Fermi = (4/3)Z^3a^2 (g/mp) = "
                << factor_xRad << " MHz\n";
      std::cout << "(alpha/pi)A_Fermi = " << factor_eF << " MHz\n";

      auto first_line_vx = IO::FRW::file_exists(fname + ".me_vx")
                               ? ""
                               : "#                h(0)         d(MLVP)       "
                                 "d(SEvx)       sum\n";

      std::ofstream of_vx;
      if (fname != "")
        of_vx.open(fname + ".me_vx", std::ios_base::app);
      of_vx << first_line_vx;

      auto vtx_in = *me_input; // has 'operator/options'
      // Add options from outer block to 'vertex block'
      for (auto &str : {"rrms", "A_vertex", "b_vertex"}) {
        const auto opt = input.getOption(str);
        if (opt)
          vtx_in.add(*opt);
      }

      const auto r_rmsfm = vtx_in.get("rrms", wf.get_rrms());
      const auto r_nucau = std::sqrt(5.0 / 3.0) * r_rmsfm / PhysConst::aB_fm;
      const auto h_MLVP =
          std::make_unique<DiracOperator::MLVP>(h.get(), *wf.rgrid, r_nucau);

      // Vertex QED term:
      const auto A_z_default = DiracOperator::VertexQED::a(wf.Znuc());
      const auto A_vertex = vtx_in.get("A_vertex", A_z_default);
      const auto b_vertex = vtx_in.get("b_vertex", 1.0);
      const auto hVertexQED = std::make_unique<DiracOperator::VertexQED>(
          h.get(), *wf.rgrid, A_vertex, b_vertex);
      if (A_vertex != 0.0) {
        std::cout << "Including effective vertex (SE) QED with: A=" << A_vertex
                  << ", b=" << b_vertex << "\n";
      }

      const auto vx_vp = ExternalField::calcMatrixElements(
          wf_VP.valence, h_MLVP.get(), nullptr, 0.0, false, diagonal_only);
      const auto vx_se = ExternalField::calcMatrixElements(
          wf_VP.valence, hVertexQED.get(), nullptr, 0.0, false, diagonal_only);

      std::cout
          << "\n             h(0)         d(MLVP)       d(SEvx)       sum\n";
      for (const auto &[a, b, x0, dv1, dv] : me0) {
        // find corresponding ME calculations in list:
        // (can't assume order will be the same, though usually will be)
        auto l = [a = a, b = b](auto e) { return (a == e.a && b == e.b); };
        // nb: a=a here due to clang error capturing bindings...
        const auto vp = std::find_if(begin(vx_vp), end(vx_vp), l);
        const auto se = std::find_if(begin(vx_se), end(vx_se), l);

        const auto Fa =
            std::find_if(begin(wf.valence), end(wf.valence),
                         [a = a](auto &x) { return x.shortSymbol() == a; });

        auto Aconst = DiracOperator::HyperfineA::convertRMEtoA(*Fa, *Fa);

        if (vp != vx_vp.end() && se != vx_se.end()) {
          const auto d_vp = vp->hab * Aconst;
          const auto d_se = se->hab * Aconst;
          const auto o =
              qip::fstring("%4s %4s:  %12.5e %12.5e  %12.5e  %12.5e\n",
                           a.c_str(), b.c_str(), x0, d_vp, d_se, d_vp + d_se);
          std::cout << o;
          of_vx << wf.Znuc() << " " << o;
        }
      }
    }

  } // me_input
} // namespace Module

//****************************************************************************

void vertexQED(const IO::InputBlock &input, const Wavefunction &wf) {
  calc_vertexQED(input, wf);
}

std::vector<std::string> calc_vertexQED(const IO::InputBlock &input,
                                        const Wavefunction &wf,
                                        const Wavefunction *wf_VP,
                                        const Wavefunction *wf_SE) {

  std::vector<std::string> out; // ??

  // wf_VP should include Vap Pol (only); wf_SE include SE only
  // If not given, assume same as wf (which should include none!)
  if (wf_VP == nullptr)
    wf_VP = &wf;
  if (wf_SE == nullptr)
    wf_SE = &wf;

  input.check({{"operator", "operator (e.g., E1 or hfs)"},
                    {"options", "operator options (same as matrixElements)"},
                    {"rrms", "nuclear rms, for QED part"},
                    {"onlyDiagonal", "only print <a|h|a>"},
                    {"radialIntegral", "false by default (means red. mat. el)"},
                    {"A_vertex", "A vtx factor; blank=default"},
                    {"b_vertex", "A vtx factor; =1 by default"},
                    {"rpa", "include RPA? NOT USED FOR NOW"},
                    {"omega", "freq. for RPA; NOT USED FOR NOW"}});

  const auto oper = input.get<std::string>("operator", "");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }
  const auto h = generateOperator(oper, h_options, wf, true);

  const bool radial_int = input.get("radialIntegral", false);

  // spacial case: HFS A (MHz)
  const bool AhfsQ = (oper == "hfs" && !radial_int);

  const auto which_str =
      radial_int ? " (radial integral)." : AhfsQ ? " A (MHz)." : " (reduced).";

  std::cout << "\n"
            << input.name() << which_str << " Operator: " << h->name() << "\n";
  std::cout << "Units: " << h->units() << "\n";

  // const bool print_both = input.get("printBoth", false);
  const bool diagonal_only = input.get("onlyDiagonal", false);

  const auto r_rmsfm = input.get("rrms", wf.get_rrms());
  const auto r_nucfm = std::sqrt(5.0 / 3.0) * r_rmsfm;
  const auto r_nucau = r_nucfm / PhysConst::aB_fm;
  const auto h_MLVP =
      std::make_unique<DiracOperator::MLVP>(h.get(), *wf.rgrid, r_nucau);

  // Vertex QED term:
  const auto A_z_default = DiracOperator::VertexQED::a(wf.Znuc());
  const auto A_vertex = input.get("A_vertex", A_z_default);
  const auto b_vertex = input.get("b_vertex", 1.0);
  const auto hVertexQED = std::make_unique<DiracOperator::VertexQED>(
      h.get(), *wf.rgrid, A_vertex, b_vertex);
  if (A_vertex != 0.0) {
    std::cout << "Including effective vertex (SE) QED with: A=" << A_vertex
              << ", b=" << b_vertex << "\n";
  }

  if (input.get("rpa", false)) {
    std::cout << "\nNote: NO RPA in vertex (yet)!!\n";
  }

  std::cout << "\nUnit conversions:\n";
  const auto a2 = std::pow(PhysConst::alpha, 2);
  const auto z3 = std::pow(wf.Znuc(), 3);
  const auto gImmPoI = h->getc();
  // gImmPoI = mu * PhysConst::muN_CGS_MHz / IN (mu here in units of muN)
  const auto factor_xRad = (4.0 / 3) * a2 * z3 * gImmPoI / PhysConst::muB_CGS;
  const auto factor_eF = (PhysConst::alpha / M_PI) * factor_xRad; //?????

  std::cout << wf.atom() << " " << wf.nuclearParams() << "\n";
  // "muN * PhysConst::muN_CGS_MHz / IN"
  // "(4/3)Z^3a^2 (g/mp)"
  std::cout << "factor_xRad = A_Fermi = (4/3)Z^3a^2 (g/mp) = " << factor_xRad
            << " MHz\n";
  std::cout << "(alpha/pi)A_Fermi = " << factor_eF << " MHz\n";

  std::cout << "\n";
  std::cout << "State    :    <h_0>         SE_vertex     MLVP\n";
  for (const auto &Fb : wf.valence) {
    for (const auto &Fa : wf.valence) {

      const auto a = AhfsQ ? DiracOperator::HyperfineA::convertRMEtoA(Fa, Fb)
                           : radial_int ? 1.0 / h->angularF(Fa.k, Fb.k) : 1.0;

      if (h->isZero(Fa.k, Fb.k))
        continue;
      if (diagonal_only && Fb != Fa)
        continue;
      if (Fb > Fa)
        continue;

      auto o = qip::fstring("%4s %4s: ", Fb.shortSymbol().c_str(),
                            Fa.shortSymbol().c_str());
      const auto h0 = h->reducedME(Fa, Fb) * a;

      // ok
      const auto Fa_SE = wf_SE->getState(Fa.n, Fa.k);
      const auto Fa_VP = wf_VP->getState(Fa.n, Fa.k);
      const auto Fb_SE = wf_SE->getState(Fb.n, Fb.k);
      const auto Fb_VP = wf_VP->getState(Fb.n, Fb.k);
      if (!Fa_SE || !Fb_SE || !Fa_VP || !Fb_VP)
        continue;

      const auto se_vtx = hVertexQED->reducedME(*Fa_SE, *Fb_SE) * a;
      const auto mlvp = h_MLVP->reducedME(*Fa_VP, *Fb_VP) * a;
      o += qip::fstring("%13.6e %13.6e %13.6e\n", h0, se_vtx, mlvp);
      // Add RPA? Might be important?
      std::cout << o;
      out.push_back(o);
    }
  }
  return out;
}

} // namespace Module
