#include "Modules/polarisability.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <fstream>
#include <iomanip>

/*
ToDO:

 * alpha and beta (transition)

*/

//==============================================================================
namespace Module {

//==============================================================================
void polarisability(const IO::InputBlock &input, const Wavefunction &wf) {

  std::cout << "\n----------------------------------------------------------\n";
  std::cout << "Calculate atomic polarisabilities at single frequency\n";

  input.check(
      {{"rpa", "Include RPA? [true]"},
       {"omega", "frequency (for single w) [0.0]"},
       {"StrucRad", "SR: include SR+Norm correction [false]"},
       {"n_min_core", "SR: Minimum n to include in SR+N [1]"},
       {"max_delta_n",
        "SR: Maximum delta n to include in SR+N sum-over-states [3]"},
       {"Qk_file",
        "SR: filename for QkTable file. If blank will not use QkTable; if "
        "exists, "
        "will read it in; if doesn't exist, will create it and write to disk. "
        "Save time (10x) at cost of memory."}});

  const auto omega = input.get("omega", 0.0);
  const auto rpaQ = input.get("rpa", true);

  const auto he1 = DiracOperator::E1(wf.grid());
  auto dVE1 = ExternalField::TDHF(&he1, wf.vHF());

  // We should use _spectrum_ for the sos - but if it is empty, just use basis
  const auto &spectrum = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

  // Solve TDHF for core, is doing RPA.
  // nb: even if not doing RPA, need TDHF object for tdhf method
  if (rpaQ) {
    dVE1.solve_core(omega);
  }

  // calculate core contribution (single omega):
  const auto ac_sos = alphaD::core_sos(wf.core(), spectrum, he1, &dVE1, omega);
  const auto ac_ms = alphaD::core_tdhf(wf.core(), he1, dVE1, omega, wf.Sigma());

  const auto eps = std::abs(2.0 * (ac_sos - ac_ms) / (ac_sos + ac_ms));
  std::cout << "\nCore polarisability (at w=" << omega << "):\n";
  std::cout << "         SOS           MS             eps\n";
  printf("Core :  %12.5e  %12.5e   [%.0e]\n", ac_sos, ac_ms, eps);

  // Valence contributions and total polarisabilities (single omega)
  if (!wf.valence().empty()) {
    std::cout << "\nValence states (at w=" << omega
              << ") [includes core contribution]:\n";
    std::cout << "         SOS           MS             eps\n";
    for (auto &Fv : wf.valence()) {
      const auto av_sos = alphaD::valence_sos(Fv, spectrum, he1, &dVE1, omega);
      const auto av_ms = alphaD::valence_tdhf(Fv, he1, dVE1, omega, wf.Sigma());
      const auto epsv = std::abs(2.0 * (ac_sos + av_sos - ac_ms - av_ms) /
                                 (ac_sos + av_sos + ac_ms + av_ms));
      printf("%4s :  %12.5e  %12.5e   [%.0e]\n", Fv.shortSymbol().c_str(),
             ac_sos + av_sos, ac_ms + av_ms, epsv);
    }
  }

  // Optionally calculate SR+N contribution
  if (input.get("StrucRad", false)) {
    const auto n_min_core = input.get("n_min_core", 1);
    const auto max_delta_n = input.get("max_delta_n", 3);
    const auto Qk_file = input.get("Qk_file", std::string{""});
    std::vector<std::pair<std::string, double>> sr_summary;
    for (const auto &Fv : wf.valence()) {
      const auto srn_v = alphaD::valence_SRN(
          Fv, spectrum, he1, &dVE1, omega, n_min_core, max_delta_n, wf.basis(),
          wf.en_coreval_gap(), Qk_file);
      sr_summary.push_back({Fv.shortSymbol(), srn_v});
    }
    // print summary:
    std::cout << "\nSummary of SR+N contributions:\n";
    for (const auto &[symbol, srn] : sr_summary) {
      printf("%4s :  %12.5e\n", symbol.c_str(), srn);
    }
  }
}

//==============================================================================
void dynamicPolarisability(const IO::InputBlock &input,
                           const Wavefunction &wf) {

  std::cout << "\n----------------------------------------------------------\n";
  std::cout << "Calculate atomic dynamic polarisabilities\n";

  input.check(
      {{"rpa", "Include RPA? [true]"},
       {"num_steps", "number of steps for dynamic polarisability [10]"},
       {"omega_minmax",
        "list (frequencies): omega_min, omega_max (in au) [0.01, 0.1]"},
       {"lambda_minmax", "list (wavelengths, will override omega_minmax): "
                         "lambda_min, lambda_max (in nm) [600, 1800]"},
       {"method",
        "Method used for dynamic pol. Either 'SOS' (sum-over-states) or 'MS' "
        "(mixed-states=TDHF). MS can be unstable for dynamic pol. [SOS]"},
       {"filename", "output filename for dynamic polarisability (if blank, "
                    "will not wite to file) []"}});

  const auto rpaQ = input.get("rpa", true);

  const auto he1 = DiracOperator::E1(wf.grid());
  auto dVE1 = ExternalField::TDHF(&he1, wf.vHF());

  // We should use _spectrum_ for the sos - but if it is empty, just use basis
  const auto &spectrum = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

  //-------------------------------------------------
  // dynamic polarisability:
  const auto num_w_steps = input.get("num_steps", 10);

  std::vector<double> w_list;
  if (input.get("lambda_minmax") != std::nullopt) {
    const auto l_minmax = input.get("lambda_minmax", std::vector{600, 1800});

    const auto w_min =
        (2.0 * M_PI * PhysConst::c / l_minmax.at(1)) * PhysConst::aB_nm;
    const auto w_max =
        (2.0 * M_PI * PhysConst::c / l_minmax.at(0)) * PhysConst::aB_nm;
    w_list = qip::uniform_range(w_min, w_max, num_w_steps);
  } else {
    const auto w_minmax = input.get("omega_minmax", std::vector{0.01, 0.1});
    w_list = qip::uniform_range(w_list.at(0), w_list.at(1), num_w_steps);
  }

  if (rpaQ) {
    // solve RPA using all iterations for initial frequency
    // then, solve the rest with just a few iterations
    dVE1.solve_core(w_list.front());
  }

  // Parse method to use for dynamic pol:
  using namespace std::string_literals;
  const auto method = input.get("method", "SOS"s);
  if (method == "SOS") {
    std::cout << "Using sum-over-states method for a(w)\n";
  } else if (method == "MS") {
    std::cout << "Using Mixed-States (TDHF) method for a(w)\n";
    std::cout << "Warning: can be unstable\n";
  } else {
    std::cout << "\nWARNING: unkown method: " << method
              << ". Available options are 'MS' or 'SOS'\n";
    std::cout << "Defaulting to SOS\n\n";
  }

  // Setup output optional file
  const auto of_name = input.get("filename", ""s);
  std::ofstream ofile;
  if (!of_name.empty()) {
    ofile.open(of_name);
    ofile << std::fixed << std::setprecision(9);
    std::cout << "Writing dynamic polarisability to file: " << of_name << "\n";
  }

  // Calculate dynamic polarisability and write to screen+file
  std::string title = " w(au)     lamda(nm) core      ";
  for (auto &Fv : wf.valence()) {
    title += (Fv.shortSymbol() + "       ");
  }
  if (rpaQ) {
    title += "eps(dV)";
  }
  std::cout << title << "\n";
  ofile << title << "\n";
  for (auto ww : w_list) {
    const auto lambda = (2.0 * M_PI * PhysConst::c / ww) * PhysConst::aB_nm;

    if (rpaQ) {
      if (dVE1.get_eps() > 1.0e-2) {
        // if tdhf didn't converge well last time, start from scratch
        // (otherwise, start from where we left off, since much faster)
        dVE1.clear();
        dVE1.solve_core(ww, 128, false);
      } else {
        dVE1.solve_core(ww, 5, false);
      }
    }
    // const auto ac =
    //     method == "MS" ?
    //         core_tdhf(wf.core(), he1, dVE1, ww, wf.Sigma()) :
    //         core_sos(wf.core(), spectrum, he1, &dVE1, ww);
    // MS method is fine for the core, and _much_ faster, and core contributed
    // negligably..so fine.
    const auto ac = alphaD::core_tdhf(wf.core(), he1, dVE1, ww, wf.Sigma());
    printf("%9.2e %9.2e %9.2e ", ww, lambda, ac);
    ofile << ww << " " << lambda << " " << ac << " ";
    std::vector<double> avs(wf.valence().size());
#pragma omp parallel for
    for (auto iv = 0ul; iv < wf.valence().size(); ++iv) {
      const auto &Fv = wf.valence().at(iv);
      const auto av =
          (ww < -Fv.en()) ?
              (method == "MS" ?
                   ac + alphaD::valence_tdhf(Fv, he1, dVE1, ww, wf.Sigma()) :
                   ac + alphaD::valence_sos(Fv, spectrum, he1, &dVE1, ww)) :
              0.0;
      avs.at(iv) = av;
    }
    for (auto &av : avs) {
      printf("%9.2e ", av);
      ofile << av << " ";
    }
    if (rpaQ) {
      printf("[%.0e]", dVE1.get_eps());
      ofile << dVE1.get_eps();
    }
    std::cout << "\n";
    ofile << "\n";
  }
}

//==============================================================================
//==============================================================================
namespace alphaD {

// J. Mitroy, M. S. Safronova, and C. W. Clark, Theory and Applications of
// Atomic and Ionic Polarizabilities, J. Phys. B 43, 44 (2010).

//==============================================================================
// Calculates polarisability of atomic core, using some-over-states method.
// CorePolarisation (dVE1) is optional - assumed to be solved already
double core_sos(const std::vector<DiracSpinor> &core,
                const std::vector<DiracSpinor> &spectrum,
                const DiracOperator::E1 &he1,
                const ExternalField::CorePolarisation *dVE1, double omega) {

  const auto f = (-2.0 / 3.0); // more general?

  // XXX For the "core" part - should we use HF core?
  // OR, we should use core part from spectrum?
  // i.e...should be orthog to valence!?

  // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states
  auto alpha_core = 0.0;
#pragma omp parallel for reduction(+ : alpha_core)
  for (std::size_t ic = 0; ic < core.size(); ++ic) {
    const auto &Fc = core.at(ic);
    const auto x = Fc.occ_frac(); // not sure if this works for open-shells
    for (const auto &Fn : spectrum) {
      if (he1.isZero(Fc.kappa(), Fn.kappa()))
        continue;
      const auto d_nc = he1.reducedME(Fn, Fc);
      const auto dv_nc = dVE1 ? dVE1->dV(Fn, Fc) : 0.0;
      // For core, only have dV for one ME
      const auto de = Fc.en() - Fn.en();
      const auto denom =
          omega == 0.0 ? 1.0 / de : de / (de * de - omega * omega);
      alpha_core += x * std::abs(d_nc * (d_nc + dv_nc)) * denom;
    }
  }

  return f * alpha_core;
}

//------------------------------------------------------------------------------
double valence_sos(const DiracSpinor &Fv,
                   const std::vector<DiracSpinor> &spectrum,
                   const DiracOperator::E1 &he1,
                   const ExternalField::CorePolarisation *dVE1, double omega) {

  const auto f = (-2.0 / 3.0) / (Fv.twoj() + 1);

  auto alpha_v = 0.0;
#pragma omp parallel for reduction(+ : alpha_v)
  for (std::size_t in = 0; in < spectrum.size(); ++in) {
    const auto &Fn = spectrum.at(in);
    if (he1.isZero(Fv.kappa(), Fn.kappa()))
      continue;
    const auto d0_nv = he1.reducedME(Fn, Fv);
    const auto dv_nv = dVE1 ? dVE1->dV(Fn, Fv) : 0.0;
    const auto d_nv = d0_nv + dv_nv;
    const auto de = Fv.en() - Fn.en();
    const auto denom = omega == 0.0 ? 1.0 / de : de / (de * de - omega * omega);
    alpha_v += std::abs(d_nv * d_nv) * denom;
  }

  return f * alpha_v;
}

//==============================================================================
// Calculates polarisability of atomic core, using TDHF (mixed states) method.
// TDHF (dVE1) is required - assumed to be solved already. If it's not solved,
// equivilant to no RPA.?
double core_tdhf(const std::vector<DiracSpinor> &core,
                 const DiracOperator::E1 &he1, const ExternalField::TDHF &dVE1,
                 double omega, const MBPT::CorrelationPotential *const Sigma) {

  // V. A. Dzuba, J. C. Berengut, J. S. M. Ginges, and V. V. Flambaum, Screening
  // of an Oscillating External Electric Field in Atoms, Phys. Rev. A 98, 043411
  // (2018).

  // NOTE: There is a question as to _which_ sigma should be used here.
  // i.e., at which energy Sigma should be solved at
  // This will use the first sigma of correct kappa (which is probably fine)
  // ...*but* Sigma should _probably_ be evaluated at valence energy

  const auto f = (-1.0 / 3.0); // More general?

  auto alpha_core = 0.0;
#pragma omp parallel for reduction(+ : alpha_core)
  for (std::size_t ic = 0; ic < core.size(); ++ic) {
    const auto &Fc = core.at(ic);
    const auto x = Fc.occ_frac(); // not sure if this works for open-shells
    // this will include dV
    const auto Xc =
        dVE1.solve_dPsis(Fc, omega, ExternalField::dPsiType::X, Sigma);
    const auto Yc =
        omega == 0.0 ?
            Xc :
            dVE1.solve_dPsis(Fc, omega, ExternalField::dPsiType::Y, Sigma);
    for (const auto &Xbeta : Xc) {
      // no dV here (for closed-shell core)
      alpha_core += x * he1.reducedME(Xbeta, Fc);
      // alpha_core += he1.reducedME(Fc, Xbeta);
    }
    for (const auto &Ybeta : Yc) {
      alpha_core += x * he1.reducedME(Ybeta, Fc);
    }
  }
  return f * alpha_core;
}

//------------------------------------------------------------------------------
double valence_tdhf(const DiracSpinor &Fv, const DiracOperator::E1 &he1,
                    const ExternalField::TDHF &dVE1, double omega,
                    const MBPT::CorrelationPotential *const Sigma) {

  const auto f = (-1.0 / 3.0) / (Fv.twoj() + 1);
  const auto Xv =
      dVE1.solve_dPsis(Fv, omega, ExternalField::dPsiType::X, Sigma);
  const auto Yv =
      omega == 0.0 ?
          Xv :
          dVE1.solve_dPsis(Fv, omega, ExternalField::dPsiType::Y, Sigma);

  double alpha_v = 0.0;
  for (const auto &Xbeta : Xv) {
    alpha_v += he1.reducedME(Xbeta, Fv) + dVE1.dV(Xbeta, Fv);
  }
  for (const auto &Ybeta : Yv) {
    alpha_v += he1.reducedME(Ybeta, Fv) + dVE1.dV(Ybeta, Fv);
  }
  return f * alpha_v;
}

//==============================================================================
double valence_SRN(const DiracSpinor &Fv,
                   const std::vector<DiracSpinor> &spectrum,
                   const DiracOperator::E1 &he1,
                   const ExternalField::CorePolarisation *dVE1, double omega,
                   // SR+N part:
                   int delta_n_max_sum, int n_min_core,
                   const std::vector<DiracSpinor> &hf_basis,
                   const double en_core, const std::string &Qk_fname) {

  // NOTE: Basis should be HF basis (used for MBPT), NOT spectrum

  std::cout << "\nStructure Radiation + Normalisation contribution: "
            << Fv.symbol() << "\n";

  auto alpha_v0 = 0.0;
  auto alpha_v1 = 0.0;
  auto alpha_v2 = 0.0;
  const auto f = (-2.0 / 3.0) / (Fv.twoj() + 1);

  auto sr = MBPT::StructureRad(hf_basis, en_core, {n_min_core, 99}, Qk_fname);
  std::cout << "Including core states from n>=" << n_min_core
            << " in diagrams\n";
  std::cout << "Calculating SR+N for terms with n=n(val)+/-" << delta_n_max_sum
            << " in the sum-over-states\n";

  // Should be possible to sum over external states _first_, then to SR?
  // Maybe not
  // Issue is depends on energy of legs

  std::cout << "  n    |d|       dEn      |  a0(n)     da(HF)    da(RPA)  |  "
               "%(HF)     %(RPA)\n";
  for (const auto &Fn : spectrum) {
    if (Fn.en() < en_core)
      continue;
    // Only do for terms with small delta_n
    if (std::abs(Fn.n() - Fv.n()) > delta_n_max_sum)
      continue;
    if (he1.isZero(Fv.kappa(), Fn.kappa()))
      continue;
    const auto d0 = he1.reducedME(Fn, Fv) + (dVE1 ? dVE1->dV(Fn, Fv) : 0.0);

    const auto [srn, srndV] = sr.srn(&he1, Fn, Fv, 0.0, dVE1);
    const auto d1 = d0 + srn;
    const auto d2 = d0 + srndV;

    const auto de = Fv.en() - Fn.en();
    const auto da_v0 = f * std::abs(d0 * d0) * de / (de * de - omega * omega);
    const auto da_v1 = f * std::abs(d1 * d1) * de / (de * de - omega * omega);
    const auto da_v2 = f * std::abs(d2 * d2) * de / (de * de - omega * omega);
    alpha_v0 += da_v0;
    alpha_v1 += da_v1;
    alpha_v2 += da_v2;

    printf(" %4s %9.2e %9.2e | %9.2e %9.2e %9.2e | %8.1e%% %8.1e%%\n",
           Fn.shortSymbol().c_str(), d0, de, da_v0, da_v1 - da_v0,
           da_v2 - da_v0, (da_v1 - da_v0) / da_v0 * 100,
           (da_v2 - da_v0) / da_v0 * 100);
    std::cout << std::flush;
  }

  const auto srn = (alpha_v1 - alpha_v0);
  const auto srn_x = (alpha_v2 - alpha_v0);

  std::cout << "a(main): " << alpha_v0 << "\n";
  std::cout << "StrucRad+Norm: " << Fv.symbol() << "\n";
  std::cout << "No RPA:  " << srn << " (" << srn / alpha_v0 * 100.0 << "%)\n";
  std::cout << "w/ RPA:  " << srn_x << " (" << srn_x / alpha_v0 * 100 << "%)\n";
  std::cout << std::flush;

  return srn;
}
} // namespace alphaD

} // namespace Module
