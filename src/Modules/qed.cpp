#include "Modules/qed.hpp"
#include "DiracOperator/Operators.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/matrixElements.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/RadiativePotential.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Format.hpp"
#include <numeric>
#include <string>
#include <vector>

namespace Module {
//****************************************************************************
void vertexQED(const IO::InputBlock &input, const Wavefunction &wf) {

  input.checkBlock(
      {"operator", "options", "rrms", "onlyDiagonal", "A_vertex", "b_vertex"});

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

  const auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto r_rmsfm = input.get("rrms", isotope.r_rms);
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

  std::cout << "\n";
  std::cout << "State     :    <h_0>         SE_vertex     MLVP\n";
  for (const auto &Fb : wf.valence) {
    for (const auto &Fa : wf.valence) { // Special case: HFS A:
      const auto a = AhfsQ ? DiracOperator::HyperfineA::convertRMEtoA(Fa, Fb)
                           : radial_int ? 1.0 / h->angularF(Fa.k, Fb.k) : 1.0;

      if (h->isZero(Fa.k, Fb.k))
        continue;
      if (diagonal_only && Fb != Fa)
        continue;
      if (Fb > Fa)
        continue;

      printf("%4s %4s : ", Fb.shortSymbol().c_str(), Fa.shortSymbol().c_str());
      const auto h0 = h->reducedME(Fa, Fb) * a;
      const auto se_vtx = hVertexQED->reducedME(Fa, Fb) * a;
      const auto mlvp = h_MLVP->reducedME(Fa, Fb) * a;
      printf("%13.6e %13.6e %13.6e\n", h0, se_vtx, mlvp);
      // Add RPA? Might be important?
    }
  }
}

//****************************************************************************
// Used for finding A and b for effective vertex QED operator
void hyperfine_vertex_test(const IO::InputBlock &input,
                           const Wavefunction &wf) {

  // Check input options for spelling mistakes etc.:
  input.check({},
              {"options", "A_vertex", "b_vertex", "rrms", "out_file", "units"});

  // units: maybe A0 (gives dA/A0), xrad, or eF

  std::cout << "\n"
            << "vertex QED corrections to hyperfine: Test + fit\n"
            << "Solve new wavefunction, without QED:\n";

  // Note: 'wf' should include QED (for perturbed orbitals)
  // This is the wavefunction calculated in the main routine
  // Note: As of now, this ONLY works for H-like systems
  // (can be easily updated in future for general case)

  // Double-check wf is H-like (i.e., has no 'core' states)
  if (!wf.core.empty()) {
    std::cout << "Note: only works for H-like systems!\n";
    return;
  }

  // const auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto r_rmsfm0 = input.get("rrms", wf.get_rrms());
  const auto r_nucfm = std::sqrt(5.0 / 3.0) * r_rmsfm0;
  const auto r_nucau = r_nucfm / PhysConst::aB_fm;

  const auto scale_rN = r_rmsfm0 / wf.get_rrms();

  // New wavefunction, but without QED, using same parameters as original
  Wavefunction wf0(wf.rgrid->params(), wf.get_nuclearParameters());
  Wavefunction wf_VP(wf.rgrid->params(), wf.get_nuclearParameters());
  Wavefunction wf_SE(wf.rgrid->params(), wf.get_nuclearParameters());

  wf_VP.radiativePotential(0.0, 1.0, 0.0, 0.0, 0.0, 5.0, scale_rN, {1.0});
  wf_SE.radiativePotential(0.0, 0.0, 1.0, 1.0, 1.0, 5.0, scale_rN, {1.0});

  // Solve for same valence states (but, without QED)
  // note: This would be editted to allow for HartreeFock (non-H-like systems)
  wf0.localValence(DiracSpinor::state_config(wf.valence));
  wf_VP.localValence(DiracSpinor::state_config(wf.valence));
  wf_SE.localValence(DiracSpinor::state_config(wf.valence));
  // print the new valence energies to screen:
  wf0.printValence();

  // Generate the hyperfine structure operator (use "generate_hfs" function)
  // 'h' is the hyperfine operator (without QED);
  // nb: for now, just hyperfine. Can change easily to work for any operator
  // const auto hfs_options = *input.getBlock("options");
  auto options = input.getBlock("options");
  auto hfs_options = IO::InputBlock("hfs", {});
  if (options) {
    hfs_options.add(options->options());
    auto hfs_rrms = hfs_options.get("rrms");
    if (!hfs_rrms) {
      auto rstr = qip::fstring("%.6f", r_rmsfm0);
      hfs_options.add("rrms = " + rstr + ";");
    }
  }

  const auto h = generate_hfsA(hfs_options, wf, true);

  // Form the vertex QED operator, called "hVertexQED":
  const auto A_z_default = DiracOperator::VertexQED::a(wf.Znuc());
  const auto A_x = input.get("A_vertex", A_z_default);
  const auto b_x = input.get("b_vertex", 1.0);
  const DiracOperator::VertexQED hVertexQED(h.get(), *wf.rgrid, A_x, b_x);
  const auto h_MLVP = DiracOperator::MLVP(h.get(), *wf.rgrid, r_nucau);

  std::cout << "\nIncluding effective vertex QED with: A=" << A_x
            << ", b=" << b_x << "\n";

  std::string s_units = input.get<std::string>("units", "A0");
  enum class Units { A0, xrad, eF, D };
  Units units = Units::A0;
  if (s_units == "xrad")
    units = Units::xrad;
  else if (s_units == "eF")
    units = Units::eF;
  else if (s_units == "D")
    units = Units::D;

  std::ofstream outfile;
  const auto fname = input.get<std::string>("out_file", "");
  if (fname != "")
    outfile.open(fname, std::ios_base::app); // append instead of overwrite

  const auto a2 = std::pow(PhysConst::alpha, 2);
  const auto z3 = std::pow(wf.Znuc(), 3);
  const auto gImmPoI = h->getc();
  const auto factor_xRad = (4.0 / 3) * a2 * z3 * gImmPoI / PhysConst::muB_CGS;
  const auto factor_eF = (PhysConst::alpha / M_PI) * factor_xRad; //?????

  std::cout << wf0.atom() << " " << wf0.nuclearParams() << "\n";
  if (units == Units::xrad) {
    std::cout << "Units: x_rad\n";
    std::cout << "factor_xRad = " << factor_xRad << " MHz\n";
  } else if (units == Units::eF) {
    std::cout << "Units: (alpha/pi)E_Fermi\n";
    std::cout << "(alpha/pi)E_Fermi = " << factor_eF << " MHz\n";
  } else if (units == Units::D) {
    std::cout << "Units: D = (alpha/pi)E_Fermi/n^2\n";
    std::cout << "(alpha/pi)E_Fermi = " << factor_eF << " MHz\n";
  } else {
    std::cout << "Units: dA/A0\n";
  }
  std::cout << "State   Z   A0(MHz)       PO(VP)        MLVP          PO(SE)   "
               "     Vertex(SE)\n";

  // Loop through each valence state, and calculate various QED corrections:
  for (const auto &Fv : wf0.valence) {

    // Factor to convert "reduced matrix element" to "hyperfine constant A"
    //(note: cancels in ratio)
    const auto a = DiracOperator::HyperfineA::convertRMEtoA(Fv, Fv);

    // Find corresponding state without QED: (can't assume in same order)
    // const auto &Fv0 = *wf0.getState(Fv.n, Fv.k);
    const auto &Fv_SE = *wf_SE.getState(Fv.n, Fv.k);
    const auto &Fv_VP = *wf_VP.getState(Fv.n, Fv.k);

    // Zeroth-order A (no QED)
    const auto A0 = h->reducedME(Fv, Fv) * a;
    const auto A_vp = h->reducedME(Fv_VP, Fv_VP) * a;
    const auto A_se = h->reducedME(Fv_SE, Fv_SE) * a;

    // Including perturbed orbital QED:
    // (subtract A0 to get just PO contribution)
    // nb: this is the part we need high-precission for, since A and A0 are
    // similar in magnitude
    auto n3 = std::pow(Fv.n, 3);
    auto unit_factor = (units == Units::xrad)
                           ? factor_xRad
                           : (units == Units::eF)
                                 ? factor_eF
                                 : (units == Units::D) ? factor_eF / n3 : A0;

    const auto A_po_SE = (A_se - A0) / unit_factor;
    const auto A_po_VP = (A_vp - A0) / unit_factor;
    // Just the vertex part:
    const auto A_vertex = hVertexQED.reducedME(Fv, Fv) * a / unit_factor;
    const auto A_mlvp = h_MLVP.reducedME(Fv, Fv) * a / unit_factor;

    auto str =
        qip::fstring(" %4s, %2i, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e\n",
                     Fv.shortSymbol().c_str(), wf0.Znuc(), A0, A_po_VP, A_mlvp,
                     A_po_SE, A_vertex);

    // Write to screen, AND (optionally) append to file:
    std::cout << str;
    outfile << str;
  }
}
} // namespace Module
