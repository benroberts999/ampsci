#include "Modules/polarisability.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <cassert>
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
  IO::ChronoTimer t("polarisability");

  std::cout << "Atomic polarisabilities, 𝛼, at single frequency\n";

  input.check(
      {{"rpa", "Include RPA? [true]"},
       {"omega", "frequency (for single w) [0.0]"},
       {"tensor", "Also calculate tensor alpha_2(w) (as well as a0) [false]"},
       {"drop_continuum", "Discard states from the spectrum with e>0 - these "
                          "can cause spurious resonances [false]"},
       {"drop_states",
        "List. Discard these states from the spectrum for "
        "sum-over-states for valence part of alpha, and "
        "from TDHF by orthogonality (must be in core/valence) []"},
       {"SRN", "SR: include SR+Norm correction [false]"},
       {"n_min_core", "SR: Minimum n to include in SR+N [1]"},
       {"max_n_SR",
        "SR: Maximum n to include in the sum-over-states for SR+N [9]"},
       {"Qk_file",
        "SR: filename for QkTable file. If blank will not use QkTable; if "
        "exists, "
        "will read it in; if doesn't exist, will create it and write to disk. "
        "Save time (10x) at cost of memory."}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto omega = input.get("omega", 0.0);
  const auto rpaQ = input.get("rpa", true);
  const auto do_tensor = input.get("tensor", false);

  const auto he1 = DiracOperator::E1(wf.grid());
  auto dVE1 = ExternalField::TDHF(&he1, wf.vHF());

  // We should use _spectrum_ for the sos - but if it is empty, just use basis
  auto spectrum = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

  // Solve TDHF for core, is doing RPA.
  // nb: even if not doing RPA, need TDHF object for tdhf method
  if (rpaQ) {
    dVE1.solve_core(omega);
  }
  const auto eFemi = wf.FermiLevel();
  const auto spectrum_core =
      qip::select_if(spectrum, [=](const auto &f) { return f.en() < eFemi; });
  const auto spectrum_excited =
      qip::select_if(spectrum, [=](const auto &f) { return f.en() > eFemi; });

  // calculate core contribution (single omega):
  const auto ac_sos_0 =
      alphaD::core_sos(wf.core(), spectrum, he1, &dVE1, omega);
  const auto ac_sos =
      alphaD::core_sos(spectrum_core, spectrum, he1, &dVE1, omega);
  const auto ac_ms_0 =
      alphaD::core_tdhf(wf.core(), he1, dVE1, omega, wf.Sigma());
  const auto ac_ms =
      alphaD::core_tdhf(spectrum_core, he1, dVE1, omega, wf.Sigma());

  const auto eps_0 =
      std::abs(2.0 * (ac_sos_0 - ac_ms_0) / (ac_sos_0 + ac_ms_0));
  const auto eps = std::abs(2.0 * (ac_ms - ac_sos) / (ac_ms + ac_sos));
  std::cout << "\nCore polarisability (at w=" << omega << "):\n";
  std::cout << "         SOS           MS             eps\n";
  printf("HF   :  %12.5e  %12.5e   [%.0e]\n", ac_sos_0, ac_ms_0, eps_0);
  printf("Spect:  %12.5e  %12.5e   [%.0e]\n", ac_sos, ac_ms, eps);

  // Drop states from spectrum _AFTER_ we calculate core contribution
  // This is used to find "tail", excluding resonant terms
  const auto drop_continuum = input.get("drop_continuum", false);
  const auto drop_states = input.get("drop_states", std::vector<std::string>{});
  if (drop_continuum) {
    std::cout << "Dropping continuum states (e>0) from the sum-over-states "
                 "spectrum.\n";
    auto is_continuum = [](const auto &a) { return a.en() > 0.0; };
    auto it = std::remove_if(spectrum.begin(), spectrum.end(), is_continuum);
    spectrum.erase(it, spectrum.end());
  }
  std::vector<DiracSpinor> force_orthog;
  if (!drop_states.empty()) {
    std::cout
        << "Dropping following states from spectrum for sum-over-states:\n ";
    for (const auto &state : drop_states) {
      std::cout << state << ", ";
      // first, drop from spectrum
      const auto [nn, kk] = AtomData::parse_symbol(state);
      const auto n = nn;
      const auto k = kk; // structured binding cannot be captured
      const auto is_nk = [n, k](const auto &a) {
        return a.n() == n && a.kappa() == k;
      };
      auto it = std::remove_if(spectrum.begin(), spectrum.end(), is_nk);
      spectrum.erase(it, spectrum.end());
      // then, find core/valence state to
      auto pFnk = wf.getState(n, k);
      if (pFnk)
        force_orthog.push_back(*pFnk);
    }
    std::cout << "\n";
  }

  //Separate out cv contribution
  // Note: below way *alomst* the same, but not exact.
  // Since it includes the "extra" dV in the core-valence part!
  if (!wf.valence().empty()) {
    const int n_main = wf.valence().front().n() + 2;
    std::cout << "\nSOS, with core-valence seperate (at w=" << omega << ")\n";
    std::cout << "Main includes up to n=" << n_main << "\n";
    std::cout << "         cv            main          tail          tot\n";

    for (auto &Fv : wf.valence()) {
      const auto main = qip::select_if(spectrum, [=](const auto &f) {
        return f.en() > eFemi && f.n() <= n_main;
      });
      const auto tail = qip::select_if(spectrum, [=](const auto &f) {
        return f.en() > eFemi && f.n() > n_main;
      });

      const auto a_vc =
          -alphaD::core_sos(spectrum_core, {Fv}, he1, &dVE1, omega) /
          Fv.twojp1();
      const auto av_sos2_main =
          alphaD::valence_sos(Fv, main, he1, &dVE1, omega);
      const auto av_sos2_tail =
          alphaD::valence_sos(Fv, tail, he1, &dVE1, omega);
      printf("%4s :  %12.5e  %12.5e  %12.5e  %12.5e  \n",
             Fv.shortSymbol().c_str(), a_vc, av_sos2_main, av_sos2_tail,
             ac_sos + a_vc + av_sos2_main + av_sos2_tail);
    }
  }

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
  if (input.get("SRN", false)) {
    const auto n_min_core = input.get("n_min_core", 1);
    const auto max_n_SR = input.get("max_n_SR", 9);
    const auto Qk_file = input.get("Qk_file", std::string{""});
    std::vector<std::tuple<std::string, double, double>> sr_summary;
    for (const auto &Fv : wf.valence()) {
      const auto [srn_v, srn_2] = alphaD::valence_SRN(
          Fv, spectrum, he1, &dVE1, omega, do_tensor, max_n_SR, n_min_core,
          wf.basis(), wf.FermiLevel(), Qk_file);
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
  std::cout << "\n";
}

//==============================================================================
void dynamicPolarisability(const IO::InputBlock &input,
                           const Wavefunction &wf) {
  IO::ChronoTimer t1("dynamicPolarisability");

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
       {"replace_w_valence",
        "Replace corresponding spectrum states with valence states - "
        "circumvents spectrum issue! [false]"},
       {"drop_continuum", "Discard states from the spectrum with e>0 - these "
                          "can cause spurious resonances [false]"},
       {"drop_states",
        "List. Discard these states from the spectrum for sum-over-states []"},
       {"filename", "output filename for dynamic polarisability (a0_ and/or "
                    "a2_ will be appended to start of filename) [identity.txt "
                    "(e.g., CsI.txt)]"},
       {"SRN", "SR: include SR+Norm correction [false]"},
       {"n_min_core", "SR: Minimum n to include in SR+N [1]"},
       {"max_n_SR",
        "SR: Maximum n to include in the sum-over-states for SR+N [9]"},
       {"Qk_file",
        "SR: filename for QkTable file. If blank will not use QkTable; if "
        "exists, "
        "will read it in; if doesn't exist, will create it and write to disk. "
        "Save time (10x) at cost of memory."}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto do_tensor = input.get("tensor", false);
  const auto rpaQ = input.get("rpa", true);
  const auto rpa_omegaQ = input.get("rpa_omega", true);
  const auto core_omegaQ = input.get("core_omega", true);

  const auto he1 = DiracOperator::E1(wf.grid());
  auto dVE1 = ExternalField::TDHF(&he1, wf.vHF());

  // We should use _spectrum_ for the sos - but if it is empty, just use basis
  auto spectrum = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

  const auto replace_w_valence = input.get("replace_w_valence", false);
  const auto drop_continuum = input.get("drop_continuum", false);
  const auto drop_states = input.get("drop_states", std::vector<std::string>{});
  if (replace_w_valence) {
    std::cout
        << "Replacing spectrum states with corresponding valence states\n";
    for (const auto &Fv : wf.valence()) {
      auto it = std::find(spectrum.begin(), spectrum.end(), Fv);
      *it = Fv;
    }
  }
  if (drop_continuum) {
    std::cout << "Dropping continuum states (e>0) from the sum-over-states "
                 "spectrum.\n";
    auto is_continuum = [](const auto &a) { return a.en() > 0.0; };
    auto it = std::remove_if(spectrum.begin(), spectrum.end(), is_continuum);
    spectrum.erase(it, spectrum.end());
  }
  if (!drop_states.empty()) {
    std::cout
        << "Dropping following states from spectrum for sum-over-states:\n ";
    for (const auto &state : drop_states) {
      std::cout << state << ", ";
      const auto [nn, kk] = AtomData::parse_symbol(state);
      const auto n = nn;
      const auto k = kk; // structured binding cannot be captured
      const auto is_nk = [n, k](const auto &a) {
        return a.n() == n && a.kappa() == k;
      };
      auto it = std::remove_if(spectrum.begin(), spectrum.end(), is_nk);
      spectrum.erase(it, spectrum.end());
    }
    std::cout << "\n";
  }

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

  const auto StrucRadQ = input.get("SRN", false);

  std::optional<MBPT::StructureRad> sr{std::nullopt};
  const auto max_n_SR = input.get("max_n_SR", 9);
  if (StrucRadQ) {
    const auto n_min_core = input.get("n_min_core", 1);
    const auto Qk_file = input.get("Qk_file", std::string{""});
    sr = MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min_core, 99},
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
    IO::ChronoTimer t("Build meTable");
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
        metab.add(Fn, Fc, me); // *
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
        metab.add(Fn, Fv, me); // *
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
            (ww == w_list.front() || ww == w_list.back() || (count % 20 == 0));

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
void transitionPolarisability(const IO::InputBlock &input,
                              const Wavefunction &wf) {
  IO::ChronoTimer t("transitionPolarisability");

  input.check(
      {{"transition", "List. states (e.g., 6s,6s) []"},
       {"rpa", "Include RPA? [true]"},
       {"omega", "frequency (for single w) [default: transition freq.]"},
       {"SRN", "SR: include SR+Norm correction [false]"},
       {"n_min_core", "SR: Minimum n to include in SR+N [1]"},
       {"max_n_SR",
        "SR: Maximum n to include in the sum-over-states for SR+N [9]"},
       {"Qk_file",
        "SR: filename for QkTable file. If blank will not use QkTable; if "
        "exists, will read it in; if doesn't exist, will create it and write "
        "to disk. Save time (10x) at cost of memory."}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto states = input.get("transition", std::vector<std::string>{});
  if (states.size() != 2) {
    std::cout << "Error 491 in transitionPolarisability(): transition option "
                 "must have exactly two states comma-separated\n";
    return;
  }
  const auto pv = wf.getState(states.at(0));
  const auto pw = wf.getState(states.at(1));
  if (!pv) {
    std::cout << "Error: Couldn't find state: " << states.at(0) << "?\n";
    return;
  }
  if (!pw) {
    std::cout << "Error: Couldn't find state: " << states.at(1) << "?\n";
    return;
  }
  const auto &Fv = *pv;
  const auto &Fw = *pw;

  const auto omega_default = std::abs(Fv.en() - Fw.en());
  const auto omega = input.get("omega", omega_default);
  const auto rpaQ = input.get("rpa", true);
  const auto srnQ = input.get("SRN", false);

  const auto alpha_text =
      " 𝛼(" + Fv.shortSymbol() + ", " + Fw.shortSymbol() + ")";
  const auto beta_text =
      " 𝛽(" + Fv.shortSymbol() + ", " + Fw.shortSymbol() + ")";
  const auto ratio_text =
      " 𝛼/𝛽(" + Fv.shortSymbol() + ", " + Fw.shortSymbol() + ")";

  std::cout << "Scalar + vector transition polarisabilities\n";
  if (rpaQ) {
    std::cout << "Including RPA at 𝜔 = " << omega << "\n";
  }
  if (srnQ) {
    std::cout << "Including Structure Radiation + Normalisation\n";
  }
  std::cout << "\n";

  const auto he1 = DiracOperator::E1(wf.grid());
  auto dVE1 = ExternalField::TDHF(&he1, wf.vHF());

  // We should use _spectrum_ for the sos - but if it is empty, just use basis
  const auto &spectrum = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

  // Solve TDHF for core, is doing RPA.
  // nb: even if not doing RPA, need TDHF object for tdhf method
  if (rpaQ) {
    dVE1.solve_core(omega);
    std::cout << "\n";
  }

  {
    std::cout << "Seperate core/main/tail (using SOS)\n";
    const int n_main = wf.valence().front().n() + 2;
    const auto eFemi = wf.FermiLevel();
    const auto core = qip::select_if(
        spectrum, [=](const auto &f) { return f.en() <= eFemi; });
    const auto main = qip::select_if(spectrum, [=](const auto &f) {
      return f.en() > eFemi && f.n() <= n_main;
    });
    const auto tail = qip::select_if(spectrum, [=](const auto &f) {
      return f.en() > eFemi && f.n() > n_main;
    });

    const auto a_c = alphaD::transition_sos(Fv, Fw, core, he1, &dVE1);
    const auto a_m = alphaD::transition_sos(Fv, Fw, main, he1, &dVE1);
    const auto a_t = alphaD::transition_sos(Fv, Fw, tail, he1, &dVE1);

    const auto b_c = alphaD::beta_sos(Fv, Fw, core, he1, &dVE1);
    const auto b_m = alphaD::beta_sos(Fv, Fw, main, he1, &dVE1);
    const auto b_t = alphaD::beta_sos(Fv, Fw, tail, he1, &dVE1);

    std::cout << "Main includes up to n=" << n_main << "\n";
    std::cout << "   core          main          tail          tot\n";
    printf("𝛼 %12.5e  %12.5e  %12.5e  %12.5e\n", a_c, a_m, a_t,
           a_c + a_m + a_t);
    printf("𝛽 %12.5e  %12.5e  %12.5e  %12.5e\n", b_c, b_m, b_t,
           b_c + b_m + b_t);
    std::cout << "\n";
  }

  // Valence contributions and total polarisabilities (single omega)
  std::cout
      << "                 SOS           MS(vw)        MS(wv)        eps\n";
  const auto avw_sos = alphaD::transition_sos(Fv, Fw, spectrum, he1, &dVE1);
  const auto avw_ms = alphaD::transition_tdhf(Fv, Fw, he1, dVE1, wf.Sigma());
  const auto awv_ms = alphaD::transition_tdhf(Fw, Fv, he1, dVE1, wf.Sigma());
  const auto eps1 = std::abs(2.0 * (avw_sos - avw_ms) / (avw_sos + avw_ms));
  const auto eps2 = std::abs(2.0 * (awv_ms - avw_ms) / (awv_ms + avw_ms));
  const auto eps = std::max(eps1, eps2);
  printf("  %11s  %12.5e  %12.5e  %12.5e  [%.0e]\n", alpha_text.c_str(),
         avw_sos, avw_ms, awv_ms, eps);

  const auto Bvw_sos = alphaD::beta_sos(Fv, Fw, spectrum, he1, &dVE1);
  const auto Bvw_ms = alphaD::beta_tdhf(Fv, Fw, he1, dVE1, wf.Sigma());
  const auto Bwv_ms = -alphaD::beta_tdhf(Fw, Fv, he1, dVE1, wf.Sigma());
  const auto Beps1 = std::abs(2.0 * (Bvw_sos - Bvw_ms) / (Bvw_sos + Bvw_ms));
  const auto Beps2 = std::abs(2.0 * (Bwv_ms - Bvw_ms) / (Bwv_ms + Bvw_ms));
  const auto Beps = std::max(Beps1, Beps2);
  printf("  %11s  %12.5e  %12.5e  %12.5e  [%.0e]\n", beta_text.c_str(), Bvw_sos,
         Bvw_ms, Bwv_ms, Beps);

  const auto R1 = avw_sos / Bvw_sos;
  const auto R2 = avw_ms / Bvw_ms;
  const auto R3 = awv_ms / Bwv_ms;
  const auto Reps1 = std::abs(2.0 * (R1 - R2) / (R1 + R2));
  const auto Reps2 = std::abs(2.0 * (R2 - R3) / (R2 + R3));
  const auto Reps = std::max(Reps1, Reps2);
  printf("%13s  %12.5e  %12.5e  %12.5e  [%.0e]\n", ratio_text.c_str(), R1, R2,
         R3, Reps);

  std::cout << "\n";

  // Optionally calculate SR+N contribution
  if (srnQ) {
    const auto n_min_core = input.get("n_min_core", 1);
    const auto max_n_SR = input.get("max_n_SR", 9);
    const auto Qk_file = input.get("Qk_file", std::string{""});

    const auto [srn_v, beta_x] = alphaD::transition_SRN(
        Fv, Fw, spectrum, he1, &dVE1, max_n_SR, n_min_core, wf.basis(),
        wf.FermiLevel(), Qk_file);

    const auto pc_A = 100.0 * srn_v / ((avw_sos + avw_ms + awv_ms) / 3.0);
    const auto pc_B = 100.0 * beta_x / ((Bvw_sos + Bvw_ms + Bwv_ms) / 3.0);
    std::cout
        << "\nInclude SR:      SOS           MS(vw)        MS(wv)        %\n";

    printf("  %11s  %12.5e  %12.5e  %12.5e  %8.1e%%\n", alpha_text.c_str(),
           avw_sos + srn_v, avw_ms + srn_v, awv_ms + srn_v, pc_A);
    printf("  %11s  %12.5e  %12.5e  %12.5e  %8.1e%%\n", beta_text.c_str(),
           Bvw_sos + beta_x, Bvw_ms + beta_x, Bwv_ms + beta_x, pc_B);

    const auto R1sr = (avw_sos + srn_v) / (Bvw_sos + beta_x);
    const auto R2sr = (avw_ms + srn_v) / (Bvw_ms + beta_x);
    const auto R3sr = (awv_ms + srn_v) / (Bwv_ms + beta_x);
    const auto pc_R = 100.0 * (R1sr - R1) / ((R1 + R2 + R3) / 3.0);
    printf("%13s  %12.5e  %12.5e  %12.5e  %8.1e%%\n", ratio_text.c_str(), R1sr,
           R2sr, R3sr, pc_R);

    std::cout << "\n";
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

  const auto f = (-2.0 / 3.0);

  // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states
  auto alpha_core = 0.0;
#pragma omp parallel for reduction(+ : alpha_core)
  for (std::size_t ic = 0; ic < core.size(); ++ic) {
    const auto &Fc = core.at(ic);
    // const auto x = Fc.occ_frac(); // not sure if this works for open-shells
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
      alpha_core += std::abs(d_nc * (d_nc + dv_nc)) * denom;
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
double transition_sos(const DiracSpinor &Fv, const DiracSpinor &Fw,
                      const std::vector<DiracSpinor> &spectrum,
                      const DiracOperator::E1 &he1,
                      const ExternalField::CorePolarisation *dVE1) {
  const auto two_m = 1; // assumes m=1/2

  double alpha_ss = 0.0;
  for (const auto &n : spectrum) {
    if (he1.isZero(Fv, n) || he1.isZero(Fw, n))
      continue;
    const auto d_vn = he1.reducedME(Fv, n) + (dVE1 ? dVE1->dV(Fv, n) : 0.0);
    const auto d_nw = he1.reducedME(n, Fw) + (dVE1 ? dVE1->dV(n, Fw) : 0.0);
    const auto f_de = 1.0 / (Fv.en() - n.en()) + 1.0 / (Fw.en() - n.en());
    const auto f = he1.rme3js(Fv.twoj(), n.twoj(), two_m) *
                   he1.rme3js(n.twoj(), Fw.twoj(), two_m);
    alpha_ss += f * d_vn * d_nw * f_de;
  }
  return alpha_ss;
}

//------------------------------------------------------------------------------
double beta_sos(const DiracSpinor &Fv, const DiracSpinor &Fw,
                const std::vector<DiracSpinor> &spectrum,
                const DiracOperator::E1 &he1,
                const ExternalField::CorePolarisation *dVE1) {
  assert(Fv.kappa() == Fw.kappa() && Fv.kappa() == -1); //only s-states

  double beta_vw = 0.0;
  for (const auto &n : spectrum) {
    if (he1.isZero(Fv, n) || he1.isZero(Fw, n))
      continue;
    const auto d_vn = he1.reducedME(Fv, n) + (dVE1 ? dVE1->dV(Fv, n) : 0.0);
    const auto d_nw = he1.reducedME(n, Fw) + (dVE1 ? dVE1->dV(n, Fw) : 0.0);
    const auto f_de = 1.0 / (Fv.en() - n.en()) - 1.0 / (Fw.en() - n.en());
    const auto f = 1.0 / (3.0 * n.twojp1());
    beta_vw += f * d_vn * d_nw * f_de;
  }
  return beta_vw;
}

//------------------------------------------------------------------------------
double transition_tdhf(const DiracSpinor &Fv, const DiracSpinor &Fw,
                       const DiracOperator::E1 &he1,
                       const ExternalField::TDHF &dVE1,
                       const MBPT::NewSigma *const Sigma) {
  const auto two_m = 1; // assumes m=1/2

  const auto X_v = dVE1.solve_dPsis(Fv, 0.0, ExternalField::dPsiType::X, Sigma);
  const auto X_w = dVE1.solve_dPsis(Fw, 0.0, ExternalField::dPsiType::X, Sigma);

  double alpha_ss = 0.0;
  for (const auto &x_w : X_w) {
    const auto f = he1.rme3js(Fv.twoj(), x_w.twoj(), two_m) *
                   he1.rme3js(x_w.twoj(), Fw.twoj(), two_m);
    alpha_ss += f * (he1.reducedME(Fv, x_w) + dVE1.dV(Fv, x_w));
  }
  for (const auto &x_v : X_v) {
    const auto f = he1.rme3js(Fw.twoj(), x_v.twoj(), two_m) *
                   he1.rme3js(x_v.twoj(), Fv.twoj(), two_m);
    alpha_ss += f * (he1.reducedME(Fw, x_v) + dVE1.dV(Fw, x_v));
  }
  return alpha_ss;
}

//------------------------------------------------------------------------------
double beta_tdhf(const DiracSpinor &Fv, const DiracSpinor &Fw,
                 const DiracOperator::E1 &he1, const ExternalField::TDHF &dVE1,
                 const MBPT::NewSigma *const Sigma) {
  assert(Fv.kappa() == Fw.kappa() && Fv.kappa() == -1); //only s-states

  const auto dv = dVE1.solve_dPsis(Fv, 0.0, ExternalField::dPsiType::X, Sigma);
  const auto dw = dVE1.solve_dPsis(Fw, 0.0, ExternalField::dPsiType::X, Sigma);

  double beta_vw = 0.0;
  for (const auto &dw_x : dw) {
    const auto f = 1.0 / (3.0 * dw_x.twojp1());
    beta_vw -= f * (he1.reducedME(Fv, dw_x) + dVE1.dV(Fv, dw_x));
  }
  for (const auto &dv_x : dv) {
    const auto f = 1.0 / (3.0 * dv_x.twojp1());
    beta_vw += f * (he1.reducedME(Fw, dv_x) + dVE1.dV(Fw, dv_x));
  }
  return beta_vw;
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
                 double omega, const MBPT::NewSigma *const Sigma) {

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
    // const auto x = Fc.occ_frac(); // not sure if this works for open-shells
    // nb: "spectrum" doesn't have occ_frac!
    // this will include dV
    const auto Xc =
        dVE1.solve_dPsis(Fc, omega, ExternalField::dPsiType::X, Sigma);
    const auto Yc =
        omega == 0.0 ?
            Xc :
            dVE1.solve_dPsis(Fc, omega, ExternalField::dPsiType::Y, Sigma);
    for (const auto &Xbeta : Xc) {
      // no dV here (for closed-shell core)
      alpha_core += he1.reducedME(Xbeta, Fc);
    }
    for (const auto &Ybeta : Yc) {
      alpha_core += he1.reducedME(Ybeta, Fc);
    }
  }
  return f * alpha_core;
}

//------------------------------------------------------------------------------
double valence_tdhf(const DiracSpinor &Fv, const DiracOperator::E1 &he1,
                    const ExternalField::TDHF &dVE1, double omega,
                    const MBPT::NewSigma *const Sigma,
                    const std::vector<DiracSpinor> &force_orthog) {

  const auto f = (-1.0 / 3.0) / (Fv.twoj() + 1);
  auto Xv = dVE1.solve_dPsis(Fv, omega, ExternalField::dPsiType::X, Sigma);
  auto Yv = omega == 0.0 ?
                Xv :
                dVE1.solve_dPsis(Fv, omega, ExternalField::dPsiType::Y, Sigma);

  // Force orthogonality:
  if (!force_orthog.empty()) {
    for (const auto &Fx : force_orthog) {
      for (auto &Xbeta : Xv) {
        if (Fx.kappa() == Xbeta.kappa())
          Xbeta -= (Xbeta * Fx) * Fx;
      }
      for (auto &Ybeta : Yv) {
        if (Fx.kappa() == Ybeta.kappa())
          Ybeta -= (Ybeta * Fx) * Fx;
      }
    }
  }

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
std::pair<double, double> valence_SRN(
    const DiracSpinor &Fv, const std::vector<DiracSpinor> &spectrum,
    const DiracOperator::E1 &he1, const ExternalField::CorePolarisation *dVE1,
    double omega, bool do_tensor,
    // SR+N part:
    int max_n_SOS, int n_min_core, const std::vector<DiracSpinor> &hf_basis,
    const double en_core, const std::string &Qk_fname) {

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

//==============================================================================
std::pair<double, double> transition_SRN(
    const DiracSpinor &Fv, const DiracSpinor &Fw,
    const std::vector<DiracSpinor> &spectrum, const DiracOperator::E1 &he1,
    const ExternalField::CorePolarisation *dVE1,
    // SR+N part:
    int max_n_SOS, int n_min_core, const std::vector<DiracSpinor> &hf_basis,
    const double en_core, const std::string &Qk_fname) {
  // NOTE: Basis should be HF basis (used for MBPT), NOT spectrum

  // const auto f = 1.0 / 6.0;

  if (Fv.kappa() != -1 || Fw.kappa() != -1) {
    std::cout << "\nWARNING 578: transition_SRN formula only valid for "
                 "s-states (for now)\n";
  }

  std::cout << "Structure Radiation + Normalisation contribution: transition "
            << Fv.symbol() << "-" << Fw.symbol() << "\n";

  auto sr = MBPT::StructureRad(hf_basis, en_core, {n_min_core, 99}, Qk_fname);
  std::cout << "Including core states from n>=" << n_min_core
            << " in diagrams\n";
  std::cout << "Calculating SR+N for terms up to n=" << max_n_SOS
            << " in the sum-over-states\n";

  //std::cout << " n       da(n)        dSR(n)       a(sum)      SR(n)%\n";
  // std::cout << " n    d_vn      srn_vn    d_nw      srn_nw    da\n";
  std::cout << "n    d_vn   d_nw   alpha_n   beta_n    SRNvn  SRNnw  dAlpha   "
               "%     dBeta    %\n";
  double alpha0_ss = 0.0;
  double alpha_ss = 0.0;
  double beta0_ss = 0.0;
  double beta_ss = 0.0;
  for (const auto &n : spectrum) {
    if (n.en() < en_core)
      continue;
    if (n.n() > max_n_SOS)
      continue;
    if (he1.isZero(Fv.kappa(), n.kappa()) || he1.isZero(Fw.kappa(), n.kappa()))
      continue;

    const auto d0_vn = he1.reducedME(Fv, n) + (dVE1 ? dVE1->dV(Fv, n) : 0.0);
    const auto d0_nw = he1.reducedME(n, Fw) + (dVE1 ? dVE1->dV(n, Fw) : 0.0);

    // don't include RPA into SR+N (simpler)
    // use HF states for "legs"
    const auto v0 = std::find(hf_basis.cbegin(), hf_basis.cend(), Fv);
    const auto w0 = std::find(hf_basis.cbegin(), hf_basis.cend(), Fw);
    const auto n0 = std::find(hf_basis.cbegin(), hf_basis.cend(), n);
    assert(v0 != hf_basis.cend() && w0 != hf_basis.cend() &&
           n0 != hf_basis.cend() && "Missing relevant states from HF basis?");
    const auto [srn_vn, x1] = sr.srn(&he1, *v0, *n0, 0.0, nullptr);
    const auto [srn_nw, x2] = sr.srn(&he1, *n0, *w0, 0.0, nullptr);

    const auto d_vn = d0_vn + srn_vn;
    const auto d_nw = d0_nw + srn_nw;

    const auto f_de_a = 1.0 / (Fv.en() - n.en()) + 1.0 / (Fw.en() - n.en());
    const auto f_de_B = 1.0 / (Fv.en() - n.en()) - 1.0 / (Fw.en() - n.en());
    const auto f =
        he1.rme3js(Fv.twoj(), n.twoj(), 1) * he1.rme3js(n.twoj(), Fw.twoj(), 1);

    const auto da0 = f * d0_vn * d0_nw * f_de_a;
    const auto da = f * d_vn * d_nw * f_de_a;
    alpha0_ss += da0;
    alpha_ss += da;

    const auto fB = 1.0 / (3.0 * n.twojp1());
    const auto dB0 = fB * d0_vn * d0_nw * f_de_B;
    const auto dB = fB * d_vn * d_nw * f_de_B;
    beta0_ss += dB0;
    beta_ss += dB;

    const auto svn = std::abs(d0_vn) / d0_vn;
    const auto snw = std::abs(d0_nw) / d0_nw;

    printf("%3s %6.3f %6.3f %9.2e %9.2e %6.3f %6.3f %8.1e %4.1f%% %8.1e "
           "%4.1f%%\n",
           n.shortSymbol().c_str(), std::abs(d0_vn), std::abs(d0_nw), da0, dB0,
           srn_vn * svn, srn_nw * snw, da - da0, 100.0 * (da - da0) / da0,
           dB - dB0, 100.0 * (dB - dB0) / dB0);
  }

  const auto a_SRN = alpha_ss - alpha0_ss;
  const auto B_SRN = beta_ss - beta0_ss;
  std::cout << "\n";
  std::cout << "StrucRad+Norm: " << Fv.symbol() << " - " << Fw.symbol() << "\n";
  std::cout << "alpha(main): " << alpha0_ss << "\n";
  std::cout << "delta(SR)  :  " << a_SRN << " " << a_SRN / alpha0_ss * 100.0
            << "%\n";

  std::cout << "beta(main): " << beta0_ss << "\n";
  std::cout << "delta(SR)  :  " << B_SRN << " " << B_SRN / beta0_ss * 100.0
            << "%\n";

  std::cout << std::flush;

  return {a_SRN, B_SRN};
}
} // namespace alphaD

} // namespace Module
