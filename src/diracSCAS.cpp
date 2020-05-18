#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //for 'ExtraPotential'
#include "IO/UserInput.hpp"
#include "Maths/Interpolator.hpp"          //for 'ExtraPotential'
#include "Maths/NumCalc_quadIntegrate.hpp" //for 'ExtraPotential'
#include "Modules/Module_runModules.hpp"
#include "Physics/PhysConst_constants.hpp" //for fit_energies
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <string>
//
#include "MBPT/CorrelationPotential.hpp"

int main(int argc, char *argv[]) {
  IO::ChronoTimer timer("\ndiracSCAS");
  const std::string input_file = (argc > 1) ? argv[1] : "diracSCAS.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // MBPT::GMatrix A(2, false);
  // A.ff[0][0] = 0;
  // A.ff[0][1] = 1;
  // A.ff[1][0] = 10;
  // A.ff[1][1] = 11;
  // A.ff.print();
  // std::cout << "\n";
  // auto A2 = A;
  // { A2.ff = A.ff.transpose(); }
  // A2.ff.print();
  // std::cout << "\n";
  // MBPT::ComplexGMatrix B(2, false);
  // B.ff = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, A.ff);
  // B.ff.print();
  // std::cout << "\n";
  // MBPT::ComplexGMatrix B2(2, false);
  // B2.ff = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, A2.ff);
  // B2.ff.print();
  // std::cout << "\n";
  // B.mult_elements_by(B2);
  // B.ff.print();
  // std::cout << "\n";
  //
  // return 1;

  // Rean in input options file
  const IO::UserInput input(input_file);

  // Get + setup atom parameters
  auto input_ok = input.check("Atom", {"Z", "A", "varAlpha2"});
  const auto atom_Z = input.get<std::string>("Atom", "Z");
  auto atom_A = input.get("Atom", "A", -1);
  const auto var_alpha = [&]() {
    const auto varAlpha2 = input.get("Atom", "varAlpha2", 1.0);
    return (varAlpha2 > 0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Get + setup Grid parameters
  input_ok = input_ok && input.check("Grid", {"r0", "rmax", "num_points",
                                              "type", "b", "fixed_du"});
  const auto r0 = input.get("Grid", "r0", 1.0e-6);
  const auto rmax = input.get("Grid", "rmax", 120.0);
  const auto du_tmp =
      input.get("Grid", "fixed_du", -1.0); // >0 means calc num_points
  const auto num_points =
      (du_tmp > 0) ? 0 : input.get("Grid", "num_points", 1600ul);
  const auto b = input.get("Grid", "b", 0.33 * rmax);
  const auto grid_type =
      (b <= r0 || b >= rmax)
          ? "logarithmic"
          : input.get<std::string>("Grid", "type", "loglinear");

  // Get + setup nuclear parameters
  input_ok =
      input_ok && input.check("Nucleus", {"A", "rrms", "skin_t", "type"});
  atom_A = input.get("Nucleus", "A", atom_A); // over-writes "atom" A
  const auto nuc_type = input.get<std::string>("Nucleus", "type", "Fermi");
  const auto rrms = input.get("Nucleus", "rrms", -1.0); // <0 means get default
  const auto skint = input.get("Nucleus", "skin_t", -1.0);

  // Create wavefunction object
  Wavefunction wf({num_points, r0, rmax, b, grid_type, du_tmp},
                  {atom_Z, atom_A, nuc_type, rrms, skint}, var_alpha);

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n";

  // Parse input for HF method
  input_ok =
      input_ok &&
      input.check("HartreeFock", {"core", "valence", "convergence", "method",
                                  "orthonormaliseValence", "sortOutput"});
  if (!input_ok)
    return 1;
  const auto str_core = input.get<std::string>("HartreeFock", "core", "[]");
  const auto eps_HF = input.get("HartreeFock", "convergence", 1.0e-12);
  const auto HF_method =
      input.get<std::string>("HartreeFock", "method", "HartreeFock");

  if (HF_method == "Hartree")
    std::cout << "Using Hartree Method (no Exchange)\n";

  // Inlcude QED radiatve potential
  const auto qed_ok =
      input.check("RadPot", {"RadPot", "Simple", "Ueh", "SE_h", "SE_l", "SE_m",
                             "rcut", "scale_rN", "scale_l", "core_qed"});
  const auto include_qed = input.get("RadPot", "RadPot", false);
  const auto x_Simple = input.get("RadPot", "Simple", 0.0);
  const auto x_Ueh = input.get("RadPot", "Ueh", 0.0);
  const auto x_SEe_h = input.get("RadPot", "SE_h", 0.0);
  const auto x_SEe_l = input.get("RadPot", "SE_l", 0.0);
  const auto x_SEm = input.get("RadPot", "SE_m", 0.0);
  const auto rcut = input.get("RadPot", "rcut", 1.0);
  const auto scale_rN = input.get("RadPot", "scale_rN", 1.0);
  const auto x_spd = input.get_list("RadPot", "scale_l", std::vector{1.0});
  const bool core_qed = input.get("RadPot", "core_qed", true);

  if (include_qed && qed_ok && core_qed) {
    wf.radiativePotential(x_Simple, x_Ueh, x_SEe_h, x_SEe_l, x_SEm, rcut,
                          scale_rN, x_spd);
    std::cout << "Including QED into Hartree-Fock core (and valence)\n\n";
  }

  // Inlcude extra potential (read in from text file):
  // Note: interpolated onto grid, but NOT extrapolated (zero outside region!)
  const auto extra_ok =
      input.check("ExtraPotential", {"filename", "factor", "beforeHF"});
  const auto ep_fname =
      input.get<std::string>("ExtraPotential", "filename", "");
  const auto ep_factor = input.get("ExtraPotential", "factor", 0.0);
  const auto ep_beforeHF = input.get("ExtraPotential", "beforeHF", false);
  const auto extra_pot =
      ep_fname != "" && std::abs(ep_factor) > 0.0 && extra_ok;
  std::vector<double> Vextra;
  if (extra_pot) {
    const auto &[x, y] = IO::FRW::readFile_xy_PoV("testIn.txt");
    Vextra = Interpolator::interpolate(x, y, wf.rgrid.r);
    NumCalc::scaleVec(Vextra, ep_factor);
  }

  // Add "extra potential", before HF (core + valence)
  if (extra_pot && ep_beforeHF) {
    wf.vnuc = NumCalc::add_vectors(wf.vnuc, Vextra);
  }

  { // Solve Hartree equations for the core:
    IO::ChronoTimer t(" core");
    wf.hartreeFockCore(HF_method, str_core, eps_HF);
  }

  if (include_qed && qed_ok && !core_qed) {
    wf.radiativePotential(x_Simple, x_Ueh, x_SEe_h, x_SEe_l, x_SEm, rcut,
                          scale_rN, x_spd);
    std::cout << "Including QED into Valence only\n\n";
  }

  // Add "extra potential", after HF (only valence)
  if (extra_pot && !ep_beforeHF) {
    wf.vdir = NumCalc::add_vectors(wf.vdir, Vextra);
  }

  // Adds effective polarision potential to direct potential
  // (After HF core, before HF valence)
  const auto Vpol_ok = input.check("dVpol", {"a_eff", "r_cut"});
  const auto a_eff = input.get("dVpol", "a_eff", 0.0);
  if (std::abs(a_eff) > 0.0 && Vpol_ok) {
    const auto r_cut = input.get("dVpol", "r_cut", 1.0);
    const auto a4 = r_cut * r_cut * r_cut * r_cut;
    auto dV = [=](auto x) { return -0.5 * a_eff / (x * x * x * x + a4); };
    for (auto i = 0u; i < wf.rgrid.num_points; ++i) {
      wf.vdir[i] += dV(wf.rgrid.r[i]);
    }
  }

  // Solve for the valence states:
  const auto valence_list =
      (wf.Ncore() < wf.Znuc())
          ? input.get<std::string>("HartreeFock", "valence", "")
          : "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad happens if below are called
    IO::ChronoTimer t("  val");
    wf.hartreeFockValence(valence_list);
    if (input.get("HartreeFock", "orthonormaliseValence", false))
      wf.orthonormaliseOrbitals(wf.valence, 2);
  }

  // Output Hartree Fock energies:
  std::cout << "\nHartree Fock: " << wf.atom() << "\n";
  const auto sorted = input.get("HartreeFock", "sortOutput", true);
  wf.printCore(sorted);
  wf.printValence(sorted);

  // Construct B-spline basis:
  const auto basis_ok =
      input.check("Basis", {"number", "order", "r0", "r0_eps", "rmax", "states",
                            "print", "positron"});
  if (basis_ok)
    wf.formBasis({input.get("Basis")});
  if (input.get("Basis", "print", false) && !wf.basis.empty()) {
    std::cout << "Basis:\n";
    wf.printBasis(wf.basis);
  }

  // Correlations: read in options
  const auto Sigma_ok = input.check(
      "Correlations", {"Brueckner", "energyShifts", "n_min_core", "fitTo_cm",
                       "lambda_k", "io_file", "stride"});
  const bool do_energyShifts = input.get("Correlations", "energyShifts", false);
  const bool do_brueckner = input.get("Correlations", "Brueckner", false);
  const auto sigma_stride = input.get("Correlations", "stride", 4);
  const auto sigma_file = input.get<std::string>("Correlations", "io_file", "");

  const auto n_min_core = input.get("Correlations", "n_min_core", 1);
  auto fit_energies =
      input.get_list("Correlations", "fitTo_cm", std::vector<double>{});
  // energies given in cm^-1, convert to au:
  NumCalc::scaleVec(fit_energies, 1.0 / PhysConst::Hartree_invcm);
  const auto lambda_k =
      input.get_list("Correlations", "lambda_k", std::vector<double>{});

  // Form correlation potential:
  if (do_energyShifts || do_brueckner) {
    IO::ChronoTimer t("Sigma");
    wf.formSigma(n_min_core, do_brueckner, sigma_stride, lambda_k, sigma_file);
  }

  auto Sk = wf.getSigma();

  if (false) {
    Sk->fill_qhat();
    std::cout << "\nTest (real-valued) Green's function:\n";
    std::cout << "(should read: 1.0  0.0)\n";
    double w = -0.3;
    std::cout << "w = " << w << "\n";
    for (auto kappa : {-1, 1, -2, 2, -3}) {

      const auto g0 = Sk->Green_hf(kappa, w).get_real();
      for (const auto &orbs : {&wf.core, &wf.valence}) {
        for (const auto &a : *orbs) {
          if (a.k != kappa)
            continue;
          const auto pfv = orbs == &wf.core ? wf.firstValenceState(kappa)
                                            : wf.firstCoreState(kappa);
          if (!pfv)
            continue;
          const auto &fv = *pfv;
          auto dv = Sk->Sigma_G_Fv(g0, a);
          auto nn = (a * dv) * (w - a.en);
          auto n2 = (fv * dv);
          auto s1 = std::string("<" + a.shortSymbol() + "|g");
          auto s2 = std::string("|" + a.shortSymbol() + ">");
          auto s3 = std::string("|" + fv.shortSymbol() + ">");
          printf("%s%s: %7.5f , %s: %6.1e\n", s1.c_str(), s2.c_str(), nn,
                 s3.c_str(), n2);
        }
      }
    }

    std::cout << "\nTest (complex-valued) Green's function:\n";
    std::cout << "(should read: 1.0+1.0i ,  0.0+0.0i)\n";
    double wr = -0.3;
    std::cout << "Re(w) = " << wr << "\n";
    for (auto wi : {0.01, 0.1, 1.0, 10.0}) {
      std::cout << "Im(w) = " << wi << "\n";

      for (auto kappa : {-1, 1, -2, 2, -3}) {
        const auto g0 = Sk->Green_hf(kappa, wr);
        const auto g = Sk->ComplexG(g0, wi);
        // const auto g = Sk->Green_core(kappa, wr, wi);
        // std::cout << "Green_core:\n";

        const auto REg = g.get_real();
        const auto IMg = g.get_imaginary();

        for (const auto &orbs : {&wf.core, &wf.valence}) {
          for (const auto &a : *orbs) {
            if (a.k != kappa)
              continue;
            const auto pfv = orbs == &wf.core ? wf.firstValenceState(kappa)
                                              : wf.firstCoreState(kappa);
            if (!pfv)
              continue;
            const auto &fv = *pfv;
            auto dvR = Sk->Sigma_G_Fv(REg, a);
            auto dvI = Sk->Sigma_G_Fv(IMg, a);
            // auto fR = (wr - a.en) + wi / (wr - a.en);
            // auto fI = (wi == 0) ? 1.0 : -((wr - a.en) * (wr - a.en) / wi +
            // wi);
            auto a2b2 = (wr - a.en) * (wr - a.en) + wi * wi;
            auto fR = a2b2 / (wr - a.en);
            auto fI = (wi == 0) ? 1.0 : -a2b2 / wi;
            double R1 = a * dvR; // * fR;
            double R0 = fv * dvR;
            double I1 = a * dvI; // * fI;
            double I0 = fv * dvI;
            auto s1 = std::string("<" + a.shortSymbol() + "|g");
            auto s2 = std::string("|" + a.shortSymbol() + ">");
            auto s3 = std::string("|" + fv.shortSymbol() + ">");
            printf("%s%s: %7.5f +%8.5f , %s:%8.1e +%8.1e", s1.c_str(),
                   s2.c_str(), R1 * fR, I1 * fI, s3.c_str(), R0, I0);
            if (std::abs(R1 * fR + I1 * fI - 2.0) > 0.2)
              std::cout << " *";
            if (std::abs(R1 * fR + I1 * fI - 2.0) > 0.35)
              std::cout << "*";
            if (std::abs(R1 * fR + I1 * fI - 2.0) > 0.75)
              std::cout << "*";
            std::cout << "\n";
            // printf("%10.3e  %10.3ei  ---  ", R1, I1);
            // printf("%10.3e  %10.3ei\n\n", 1.0 / fR, 1.0 / fI);
          }
        }
      }
    }
  }

  {
    IO::ChronoTimer t("FeynmanDirect");
    Sk->fill_qhat();
    // Sk->FeynmanDirect(-1, -0.127368);
    Sk->FeynmanEx_1(-1, -0.127368);
    // Sk->FeynmanDirect(1, -0.085616);
    // Sk->FeynmanDirect(-2, -0.083785);
    // Sk->FeynmanDirect(2, -0.064419651);
    // Sk->FeynmanDirect(-3, -0.064529789);
  }

  // Just energy shifts
  if (!wf.valence.empty() && do_energyShifts && Sigma_ok) {
    IO::ChronoTimer t("de");
    wf.SOEnergyShift();
  }
  // Brueckner orbitals (optionally, fit Sigma to exp energies)
  if (!wf.valence.empty() && do_brueckner && Sigma_ok) {
    IO::ChronoTimer t("Br");
    if (!fit_energies.empty())
      wf.fitSigma_hfBrueckner(valence_list, fit_energies);
    else
      wf.hartreeFockBrueckner();
  }
  // Print out info for new "Brueckner" valence orbitals:
  if (!wf.valence.empty() && do_brueckner && Sigma_ok) {
    std::cout << "\nBrueckner orbitals:\n";
    wf.printValence(sorted);
  }

  // Construct B-spline Spectrum:
  const auto spectrum_ok =
      input.check("Spectrum", {"number", "order", "r0", "r0_eps", "rmax",
                               "states", "print", "positron"});
  if (spectrum_ok)
    wf.formSpectrum({input.get("Spectrum")});
  if (input.get("Spectrum", "print", false) && !wf.spectrum.empty()) {
    std::cout << "Spectrum:\n";
    wf.printBasis(wf.spectrum);
  }

  // run each of the modules with the calculated wavefunctions
  Module::runModules(input, wf);

  return 0;
}

//******************************************************************************
