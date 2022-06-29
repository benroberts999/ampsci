#include "Modules/polarisability.hpp"
#include "Coulomb/meTable.hpp"
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
       {"tensor", "Also calculate tensor alpha_2(w) (as well as a0) [false]"},
       {"StrucRad", "SR: include SR+Norm correction [false]"},
       {"n_min_core", "SR: Minimum n to include in SR+N [1]"},
       {"max_n_SR",
        "SR: Maximum n to include in the sum-over-states for SR+N [9]"},
       {"Qk_file",
        "SR: filename for QkTable file. If blank will not use QkTable; if "
        "exists, "
        "will read it in; if doesn't exist, will create it and write to disk. "
        "Save time (10x) at cost of memory."}});

  const auto omega = input.get("omega", 0.0);
  const auto rpaQ = input.get("rpa", true);
  const auto do_tensor = input.get("tensor", false);

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

  if (do_tensor && !wf.valence().empty()) {
    std::cout << "\nTensor polarasability alpha_2 (at w=" << omega << "):\n";
    std::cout << "         SOS\n";
    for (auto &Fv : wf.valence()) {
      const auto a2_sos = alphaD::tensor2_sos(Fv, spectrum, he1, &dVE1, omega);
      printf("%4s :  %12.5e\n", Fv.shortSymbol().c_str(), a2_sos);
    }
  }

  // Optionally calculate SR+N contribution
  if (input.get("StrucRad", false)) {
    const auto n_min_core = input.get("n_min_core", 1);
    const auto max_n_SR = input.get("max_n_SR", 9);
    const auto Qk_file = input.get("Qk_file", std::string{""});
    std::vector<std::tuple<std::string, double, double>> sr_summary;
    for (const auto &Fv : wf.valence()) {
      const auto [srn_v, srn_2] = alphaD::valence_SRN(
          Fv, spectrum, he1, &dVE1, omega, max_n_SR, n_min_core, wf.basis(),
          wf.en_coreval_gap(), Qk_file);
      sr_summary.push_back({Fv.shortSymbol(), srn_v, srn_2});
    }
    // print summary:
    std::cout << "\nSummary of SR+N contributions:\n";
    for (const auto &[symbol, srn, srn_2] : sr_summary) {
      printf("%4s :  %12.5e", symbol.c_str(), srn);
      if (do_tensor) {
        printf("  %12.5e", srn_2);
      }
      std::cout << "\n";
    }
  }
}

//==============================================================================
void dynamicPolarisability(const IO::InputBlock &input,
                           const Wavefunction &wf) {

  std::cout << "\n----------------------------------------------------------\n";
  std::cout << "Calculate atomic dynamic polarisabilities\n";

  input.check(
      {{"tensor", "Do tensor polarisability a2(w) (as well as a0) [false]"},
       {"rpa", "Include RPA? [true]"},
       {"core_omega",
        "Frequency-dependent core? If true, core part evaluated at each "
        "frequency. If false, core evaluated once at w=0 [true]"},
       {"rpa_omega", "Frequency-dependent RPA? If true, RPA solved at each "
                     "frequency. If false, RPA solved once at w=0 [true]"},
       {"num_steps", "number of steps for dynamic polarisability [10]"},
       {"omega_minmax",
        "list (frequencies): omega_min, omega_max (in au) [0.01, 0.1]"},
       {"lambda_minmax", "list (wavelengths, will override omega_minmax): "
                         "lambda_min, lambda_max (in nm) [600, 1800]"},
       {"method", "Method used for dynamic pol. for a0(w). Either 'SOS' "
                  "(sum-over-states) or 'MS' (mixed-states=TDHF). MS can be "
                  "unstable for dynamic pol. [SOS]"},
       {"filename", "output filename for dynamic polarisability (a0_ and/or "
                    "a2_ will be appended to start of filename) [identity.txt "
                    "(e.g., CsI.txt)]"},
       {"StrucRad", "SR: include SR+Norm correction [false]"},
       {"n_min_core", "SR: Minimum n to include in SR+N [1]"},
       {"max_n_SR",
        "SR: Maximum n to include in the sum-over-states for SR+N [9]"},
       {"Qk_file",
        "SR: filename for QkTable file. If blank will not use QkTable; if "
        "exists, "
        "will read it in; if doesn't exist, will create it and write to disk. "
        "Save time (10x) at cost of memory."}});

  const auto do_tensor = input.get("tensor", false);
  const auto rpaQ = input.get("rpa", true);
  const auto rpa_omegaQ = input.get("rpa_omega", true);
  const auto core_omegaQ = input.get("core_omega", true);

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
    if (w_minmax.size() != 2) {
      std::cout << "\nFailure 162 in dynamicPolarisability.\nIssue with "
                   "omega_minmax option: must have exactly 2 (or 0) entires. "
                   "Missing or extra comma?\n";
      return;
    }
    w_list = qip::uniform_range(w_minmax.at(0), w_minmax.at(1), num_w_steps);
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
  if (rpaQ && rpa_omegaQ) {
    std::cout << "Solving RPA at each frequency\n";
  } else if (rpaQ && !rpa_omegaQ) {
    std::cout << "Solving RPA once at zero frequency\n";
  } else {
    std::cout << "Not including RPA\n";
  }
  if (core_omegaQ) {
    std::cout << "Including frequency-dependent core part.\n";
  } else {
    std::cout << "Core part evaluated once at zero frequency.\n";
  }

  if (rpaQ) {
    // solve RPA using all iterations for initial frequency
    // then, solve the rest with just a few iterations
    const auto w_initial = rpa_omegaQ ? w_list.front() : 0.0;
    dVE1.solve_core(w_initial);
  }
  // if not rpa_omegaQ, then dV included in meTable
  auto *const dVptr = rpa_omegaQ ? &dVE1 : nullptr;

  // static (w=0) core part.
  const auto ac0 = core_omegaQ ?
                       0.0 :
                       alphaD::core_tdhf(wf.core(), he1, dVE1, 0.0, wf.Sigma());

  const auto StrucRadQ = input.get("StrucRad", false);

  std::optional<MBPT::StructureRad> sr{std::nullopt};
  const auto max_n_SR = input.get("max_n_SR", 9);
  if (StrucRadQ) {
    const auto n_min_core = input.get("n_min_core", 1);
    const auto Qk_file = input.get("Qk_file", std::string{""});
    sr = MBPT::StructureRad(wf.basis(), wf.en_coreval_gap(), {n_min_core, 99},
                            Qk_file);
    std::cout << "Including core states from n>=" << n_min_core
              << " in diagrams\n";
    std::cout << "Calculating SR+N for terms up to n=" << max_n_SR
              << " in the sum-over-states\n";
  }

  // build tables of matrix elements;
  // Calculate these once, saves much time for calculation
  // *note: relies on assumption that alphaD::..._sos always calls (n,v) or
  // (n,c) and never (v,n) or (c,n)
  Coulomb::meTable metab{};
  if (method != "MS") {
    std::cout << "Building table of matrix elements.." << std::flush;
    for (const auto &Fn : spectrum) {
      // core part:
      for (const auto &Fc : wf.core()) {
        if (he1.isZero(Fn, Fc))
          continue;
        auto me = he1.reducedME(Fn, Fc);
        if (rpaQ && !rpa_omegaQ) {
          me += dVE1.dV(Fn, Fc);
        }
        metab.add(Fn, Fc, me);
        // *
        // metab_c.add(Fc, Fn, he1.symm_sign(Fn, Fc) * me);
      }
      // valence part:
      for (const auto &Fv : wf.valence()) {
        if (he1.isZero(Fn, Fv))
          continue;
        auto me = he1.reducedME(Fn, Fv);
        if (rpaQ && !rpa_omegaQ) {
          me += dVE1.dV(Fn, Fv);
        }
        // Adds SR+Norm! (ignores freq. dependence)
        if (sr && Fn.n() <= max_n_SR) {
          me += sr->srn(&he1, Fn, Fv).first;
        }
        metab.add(Fn, Fv, me);
        // *
        // metab_v.add(Fv, Fn, he1.symm_sign(Fn, Fv) * me);
      }
    }
    std::cout << " Done.\n" << std::flush;
  }

  // Setup output optional file
  const auto of_name = input.get("filename", wf.identity() + ".txt");
  std::ofstream ofile, o2file;
  if (!of_name.empty()) {
    ofile.open("a0_" + of_name);
    ofile << std::scientific << std::setprecision(9);
    std::cout << "Writing dynamic polarisability a0(w) to file: a0_" << of_name
              << "\n";
    if (do_tensor) {
      o2file.open("a2_" + of_name);
      o2file << std::scientific << std::setprecision(9);
      std::cout << "Writing dynamic polarisability a2(w) to file: a2_"
                << of_name << "\n";
    }
  }

  // Calculate dynamic polarisability and write to screen+file
  std::string title = "w(au)      lamda(nm) core     ";
  for (auto &Fv : wf.valence()) {
    title += (" "s + Fv.shortSymbol() + "      "s);
  }
  if (rpaQ && rpa_omegaQ) {
    title += "eps(dV)";
  }
  std::cout << title << "\n";
  ofile << title << "\n";
  o2file << title << "\n";
  int count = 0;
  for (auto ww : w_list) {
    const auto lambda = (2.0 * M_PI * PhysConst::c / ww) * PhysConst::aB_nm;

    // if <20, print all; otherwise, first + last + every 10th
    count++;
    const auto print =
        w_list.size() < 20 ?
            true :
            (ww == w_list.front() || ww == w_list.back() || (count % 10 == 0));

    if (rpaQ && rpa_omegaQ) {
      if (dVE1.get_eps() > 1.0e-2) {
        // if tdhf didn't converge well last time, start from scratch
        // (otherwise, start from where we left off, since much faster)
        dVE1.clear();
        dVE1.solve_core(ww, 128, false);
      } else {
        dVE1.solve_core(ww, 5, false);
      }
    }
    // MS method is fine for the core, and _much_ faster, and core contributes
    // negligably..so fine.
    const auto ac =
        !core_omegaQ ?
            ac0 :
            (method == "MS" ?
                 alphaD::core_tdhf(wf.core(), he1, dVE1, ww, wf.Sigma()) :
                 alphaD::core_sos(wf.core(), spectrum, he1, dVptr, ww, &metab));

    if (print)
      printf("%9.2e %9.2e %9.2e ", ww, lambda, ac);
    ofile << ww << " " << lambda << " " << ac << " ";
    // no core contrib to a2, but write zero so columns align
    o2file << ww << " " << lambda << " " << 0.0 << " ";
    std::vector<double> avs(wf.valence().size());
    std::vector<double> a2s(wf.valence().size());
#pragma omp parallel for if (method == "MS")
    for (auto iv = 0ul; iv < wf.valence().size(); ++iv) {
      const auto &Fv = wf.valence().at(iv);
      const auto av =
          method == "MS" ?
              ac + alphaD::valence_tdhf(Fv, he1, dVE1, ww, wf.Sigma()) :
              ac + alphaD::valence_sos(Fv, spectrum, he1, dVptr, ww, &metab);
      avs.at(iv) = av;
      if (do_tensor) {
        const auto a2 =
            alphaD::tensor2_sos(Fv, spectrum, he1, dVptr, ww, &metab);
        a2s.at(iv) = a2;
      }
    }
    for (auto &av : avs) {
      if (print)
        printf("%9.2e ", av);
      ofile << av << " ";
    }
    if (rpaQ && rpa_omegaQ) {
      if (print)
        printf("[%.0e]", dVE1.get_eps());
      ofile << dVE1.get_eps();
    }
    if (print)
      std::cout << "\n";
    ofile << "\n";

    if (do_tensor) {
      if (print) {
        std::cout << "          "
                  << "          "
                  << " a2(w)    ";
      }
      for (auto &a2 : a2s) {
        if (print)
          printf("%9.2e ", a2);
        o2file << a2 << " ";
      }
      if (rpaQ && rpa_omegaQ)
        o2file << dVE1.get_eps();
      if (print)
        std::cout << "\n";
      o2file << "\n";
    }
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
                const ExternalField::CorePolarisation *dVE1, double omega,
                const Coulomb::meTable<double> *dtab) {

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
      const auto dtab_nc = dtab ? dtab->get(Fn, Fc) : nullptr;
      const auto d_nc = dtab_nc ? *dtab_nc : he1.reducedME(Fn, Fc);
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
                   const ExternalField::CorePolarisation *dVE1, double omega,
                   const Coulomb::meTable<double> *dtab) {

  const auto f = (-2.0 / 3.0) / (Fv.twoj() + 1);

  auto alpha_v = 0.0;
#pragma omp parallel for reduction(+ : alpha_v)
  for (std::size_t in = 0; in < spectrum.size(); ++in) {
    const auto &Fn = spectrum.at(in);
    if (he1.isZero(Fv.kappa(), Fn.kappa()))
      continue;
    const auto dtab_nv = dtab ? dtab->get(Fn, Fv) : nullptr;
    const auto d0_nv = dtab_nv ? *dtab_nv : he1.reducedME(Fn, Fv);
    const auto dv_nv = dVE1 ? dVE1->dV(Fn, Fv) : 0.0;
    const auto d_nv = d0_nv + dv_nv;
    const auto de = Fv.en() - Fn.en();
    const auto denom = omega == 0.0 ? 1.0 / de : de / (de * de - omega * omega);
    alpha_v += std::abs(d_nv * d_nv) * denom;
  }

  return f * alpha_v;
}

//------------------------------------------------------------------------------
double tensor2_sos(const DiracSpinor &Fv,
                   const std::vector<DiracSpinor> &spectrum,
                   const DiracOperator::E1 &he1,
                   const ExternalField::CorePolarisation *dVE1, double omega,
                   const Coulomb::meTable<double> *dtab) {

  // B. Arora, M. S. Safronova, and C. W. Clark, Phys. Rev. A 76, 052509 (2007).
  // J. Mitroy, M. S. Safronova, and C. W. Clark, Theory and Applications of
  // Atomic and Ionic Polarizabilities, J. Phys. B 43, 44 (2010).

  const auto ctop = 2.5 * (Fv.twoj() * (Fv.twoj() - 1));
  const auto cbot = 3.0 * ((Fv.twoj() + 2) * (Fv.twoj() + 1) * (Fv.twoj() + 3));
  const auto C = +4.0 * std::sqrt(ctop / cbot); // nb: diff de sign
  if (ctop <= 0.0)
    return 0.0;

  auto alpha_2 = 0.0;
#pragma omp parallel for reduction(+ : alpha_2)
  for (std::size_t in = 0; in < spectrum.size(); ++in) {
    const auto &Fn = spectrum.at(in);
    if (he1.isZero(Fv.kappa(), Fn.kappa()))
      continue;

    const auto sj = Angular::sixj_2(Fv.twoj(), 2, Fn.twoj(), 2, Fv.twoj(), 4);
    if (sj == 0.0)
      continue;
    const auto s = Angular::neg1pow_2(Fv.twoj() + Fn.twoj() + 2);

    // const auto d0_nv = he1.reducedME(Fn, Fv);
    const auto dtab_nv = dtab ? dtab->get(Fn, Fv) : nullptr;
    const auto d0_nv = dtab_nv ? *dtab_nv : he1.reducedME(Fn, Fv);
    const auto dv_nv = dVE1 ? dVE1->dV(Fn, Fv) : 0.0;
    const auto d_nv = d0_nv + dv_nv;
    const auto de = Fv.en() - Fn.en();
    const auto denom = omega == 0.0 ? 1.0 / de : de / (de * de - omega * omega);
    alpha_2 += s * sj * std::abs(d_nv * d_nv) * denom;
  }

  return C * alpha_2;
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
std::pair<double, double>
valence_SRN(const DiracSpinor &Fv, const std::vector<DiracSpinor> &spectrum,
            const DiracOperator::E1 &he1,
            const ExternalField::CorePolarisation *dVE1, double omega,
            // SR+N part:
            int max_n_SOS, int n_min_core,
            const std::vector<DiracSpinor> &hf_basis, const double en_core,
            const std::string &Qk_fname) {

  // NOTE: Basis should be HF basis (used for MBPT), NOT spectrum

  std::cout << "\nStructure Radiation + Normalisation contribution: "
            << Fv.symbol() << "\n";

  auto alpha_v0 = 0.0;
  auto alpha_v1 = 0.0;
  const auto f = (-2.0 / 3.0) / (Fv.twoj() + 1);

  // for tensor
  auto alpha2_v0 = 0.0;
  auto alpha2_v1 = 0.0;
  const auto ctop = 2.5 * (Fv.twoj() * (Fv.twoj() - 1));
  const auto cbot = 3.0 * ((Fv.twoj() + 2) * (Fv.twoj() + 1) * (Fv.twoj() + 3));
  const auto C = ctop >= 0.0 ? +4.0 * std::sqrt(ctop / cbot) : 0.0;

  auto sr = MBPT::StructureRad(hf_basis, en_core, {n_min_core, 99}, Qk_fname);
  std::cout << "Including core states from n>=" << n_min_core
            << " in diagrams\n";
  std::cout << "Calculating SR+N for terms up to n=" << max_n_SOS
            << " in the sum-over-states\n";

  // Should be possible to sum over external states _first_, then to SR?
  // Maybe not
  // Issue is depends on energy of legs

  bool do_tensor = true;

  std::cout << "  n    |d|       dEn      |  a0(n)     dSR(n)    [%]      ";
  if (do_tensor)
    std::cout << "|  a2(n)     dSR2(n)";
  std::cout << "\n";
  for (const auto &Fn : spectrum) {
    if (Fn.en() < en_core)
      continue;
    // Only do for terms with small delta_n
    if (Fn.n() > max_n_SOS)
      continue;
    if (he1.isZero(Fv.kappa(), Fn.kappa()))
      continue;
    const auto d0 = he1.reducedME(Fn, Fv) + (dVE1 ? dVE1->dV(Fn, Fv) : 0.0);

    // don't include RPA into SR+N (simpler)
    const auto [srn, srndV] = sr.srn(&he1, Fn, Fv, omega, nullptr);
    const auto d1 = d0 + srn;

    // alpha_0(w)
    const auto de = Fv.en() - Fn.en();
    const auto da_v0 = f * std::abs(d0 * d0) * de / (de * de - omega * omega);
    const auto da_v1 = f * std::abs(d1 * d1) * de / (de * de - omega * omega);
    alpha_v0 += da_v0;
    alpha_v1 += da_v1;

    printf(" %4s %9.2e %9.2e | %9.2e %9.2e %8.1e%%", Fn.shortSymbol().c_str(),
           d0, de, da_v0, da_v1 - da_v0, (da_v1 - da_v0) / da_v0 * 100);

    // alpha_2(w) [tensor]
    if (do_tensor && C != 0.0) {
      const auto sj = Angular::sixj_2(Fv.twoj(), 2, Fn.twoj(), 2, Fv.twoj(), 4);
      const auto s = Angular::neg1pow_2(Fv.twoj() + Fn.twoj() + 2);
      const auto da2_v0 =
          s * sj * C * std::abs(d0 * d0) * de / (de * de - omega * omega);
      const auto da2_v1 =
          s * sj * C * std::abs(d1 * d1) * de / (de * de - omega * omega);
      alpha2_v0 += da2_v0;
      alpha2_v1 += da2_v1;
      printf(" | %9.2e %9.2e", da2_v0, da2_v1 - da2_v0);
    }

    std::cout << '\n' << std::flush;
  }

  const auto srn = (alpha_v1 - alpha_v0);
  const auto srn2 = (alpha2_v1 - alpha2_v0);

  std::cout << "StrucRad+Norm: " << Fv.symbol() << "\n";
  std::cout << "a0(main): " << alpha_v0 << "\n";
  std::cout << "dSR_0   :  " << srn << " (" << srn / alpha_v0 * 100.0 << "%)\n";
  if (do_tensor) {
    std::cout << "a2(main): " << alpha2_v0 << "\n";
    std::cout << "dSR_2   :  " << srn2 << " (" << srn2 / alpha2_v0 * 100.0
              << "%)\n";
  }
  std::cout << std::flush;

  return {srn, srn2};
}
} // namespace alphaD

} // namespace Module
