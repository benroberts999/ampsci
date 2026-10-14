#include "Amplitudes/SecondOrder.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Modules/Modules.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

//==============================================================================
namespace Module {

void dynamicPolarisability(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_dynamicPolarisability{"dynamicPolarisability",
                                       "Calculates dynamic polarisabilities",
                                       &dynamicPolarisability};
} // namespace

//==============================================================================
void dynamicPolarisability(const IO::InputBlock &input,
                           const Wavefunction &wf) {
  IO::ChronoTimer t1("dynamicPolarisability");

  std::cout << "\n----------------------------------------------------------\n";
  std::cout << "Calculate atomic dynamic polarisabilities\n";

  input.check(
    {{"states",
      "Which states to calculate? (e.g., '7sp6d'). Must be a subset of "
      "valence. By default, all valence states are calculated."},
     {"tensor", "Do tensor polarisability a2(w) (as well as a0) [false]"},
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
     {"StructureRadiation{}",
      "Options for structure radiation and normalisation. If this block is "
      "included, SR+N is added to the single-particle matrix elements of the "
      "sum-over-states (evaluated at w = 0). Sum-over-states method only"}});

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
              "structure radiation is included [true]"},
     {"RPA_in_SR", "Include RPA in the SR+N diagrams [false]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto v_string =
    input.get("states", DiracSpinor::state_config(wf.valence()));
  const auto states = DiracSpinor::subset(wf.valence(), v_string);

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
      if (it != spectrum.end()) {
        *it = Fv;
      }
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

    // w_list = qip::uniform_range(w_min, w_max, num_w_steps);
    // logarithmic in omega is roughly ~linear in lambda
    w_list = qip::logarithmic_range(w_min, w_max, num_w_steps);

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
  // dV evaluated live from dVE1: if not rpa_omegaQ, it stays at w = 0
  auto *const dVptr = rpaQ ? &dVE1 : nullptr;

  // Mixed-states polarisabilities: alpha = A^0 / sqrt(3(2J+1)), with A^0 the
  // K = 0 second-order amplitude for t = s = E1 at (w, -w). The core part is
  // independent of the valence J (the sqrt cancels): use 2J = 1
  const auto alpha_ms_valence = [&](const DiracSpinor &Fv, double w) {
    return Amplitudes::ms_valence(0, Fv, Fv, &he1, &he1, w, -w, &dVE1, &dVE1,
                                  wf.Sigma())
             .first /
           std::sqrt(3.0 * (Fv.twoj() + 1));
  };
  const auto alpha_ms_core = [&](double w) {
    return Amplitudes::ms_core(0, 1, &he1, &he1, w, -w, wf.core(), &dVE1, &dVE1,
                               wf.Sigma()) /
           std::sqrt(6.0);
  };

  // static (w=0) core part.
  const auto ac0 = core_omegaQ ? 0.0 : alpha_ms_core(0.0);

  std::optional<MBPT::StructureRad> sr{std::nullopt};
  const auto sr_n_max =
    SR_input.get("n_max_legs", DiracSpinor::max_n(wf.core()) + 3);
  const auto sr_norm = SR_input.get("norm", true);
  const auto rpa_in_SR = rpaQ && SR_input.get("RPA_in_SR", false);
  if (t_SR_input) {
    const auto n_minmax = SR_input.get("n_minmax", std::vector{1});
    const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
    const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
    const auto Qk_file_t = SR_input.get("Qk_file", std::string{"false"});
    const std::string Qk_file =
      Qk_file_t != "false" ?
        Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t :
        "";
    std::cout << "\nIncluding structure radiation"
              << (sr_norm ? " and normalisation:\n" : " (no normalisation):\n");
    fmt::print("Applied to matrix elements with n <= {}\n", sr_n_max);
    if (rpa_in_SR) {
      std::cout << "Including RPA in SR+N\n";
    } else {
      std::cout << "Not including RPA in SR+N\n";
    }
    std::cout << std::flush;
    sr =
      MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max}, Qk_file);
    sr->solve_core(&he1, rpa_in_SR ? &dVE1 : nullptr);
  }

  // SR (+Norm) corrections: expensive, so calculated once (at w = 0,
  // ignoring their frequency dependence) and cached in a table. Only between
  // physical (low-n) states, as Amplitudes::me_table
  Coulomb::meTable<double> sr_tab{};
  if (sr) {
    IO::ChronoTimer t("Build SR table");
    std::cout << "Building table of SR matrix elements.." << std::flush;
    for (const auto &Fv : states) {
      if (Fv.n() > sr_n_max)
        continue;
      for (const auto &Fn : spectrum) {
        if (Fn.n() > sr_n_max || he1.isZero(Fn, Fv))
          continue;
        const auto *sr_dV = rpa_in_SR ? &dVE1 : nullptr;
        sr_tab.add(Fn, Fv,
                   sr_norm ? sr->srn(Fn, Fv, &he1, sr_dV, 0.0) :
                             sr->SR(Fn, Fv, 0.0));
      }
    }
    std::cout << " Done.\n" << std::flush;
  }

  // Full single-particle ME table: t0 + dV + SR, with the dV component at
  // +w (absorption, conj = false) or -w (emission, conj = true). Cheap to
  // rebuild each frequency: only the dV part varies (and the SR is cached)
  const auto build_me_table = [&](bool conj) {
    Coulomb::meTable<double> t_me{};
    for (const auto &Fv : states) {
      for (const auto &Fn : spectrum) {
        if (he1.isZero(Fn, Fv))
          continue;
        const auto dv_nv = dVptr ? dVptr->dV(Fn, Fv, conj) : 0.0;
        const auto dv_vn = dVptr ? dVptr->dV(Fv, Fn, conj) : 0.0;
        const auto sr_nv = sr_tab.getv(Fn, Fv);
        t_me.add(Fn, Fv, he1.reducedME(Fn, Fv) + dv_nv + sr_nv);
        t_me.add(Fv, Fn,
                 he1.reducedME(Fv, Fn) + dv_vn + he1.symm_sign(Fn, Fv) * sr_nv);
      }
    }
    return t_me;
  };

  // Excited (above-Fermi) part of the spectrum, for the core sum
  const auto excited =
    DiracSpinor::split_by_energy(spectrum, wf.FermiLevel()).second;

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
  for (auto &Fv : states) {
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
      if (dVE1.last_eps() > 1.0e-2) {
        // if tdhf didn't converge well last time, start from scratch
        // (otherwise, start from where we left off, since much faster)
        dVE1.clear();
        dVE1.solve_core(ww, 128, false);
      } else {
        dVE1.solve_core(ww, 5, false);
      }
    }
    // Full ME tables at this frequency (SOS paths only): the two legs of
    // the amplitude carry the field at +w and at -w
    const auto sos_tables = method != "MS" || do_tensor;
    const auto t_me =
      sos_tables ? build_me_table(false) : Coulomb::meTable<double>{};
    const auto s_me =
      sos_tables ? build_me_table(true) : Coulomb::meTable<double>{};

    // MS method is fine for the core, and _much_ faster, and core contributes
    // negligably..so fine.
    const auto ac =
      !core_omegaQ ? ac0 :
                     (method == "MS" ?
                        alpha_ms_core(ww) :
                        Amplitudes::sos_core(0, 1, &he1, &he1, ww, -ww,
                                             wf.core(), excited, dVptr, dVptr) /
                          std::sqrt(6.0));

    if (print)
      printf("%9.2e %9.2e %9.2e ", ww, lambda, ac);
    ofile << ww << " " << lambda << " " << ac << " ";
    // no core contrib to a2, but write zero so columns align
    o2file << ww << " " << lambda << " " << 0.0 << " ";
    std::vector<double> avs(states.size());
    std::vector<double> a2s(states.size());
#pragma omp parallel for if (method == "MS")
    for (auto iv = 0ul; iv < states.size(); ++iv) {
      const auto &Fv = states.at(iv);
      // alpha = A^0 / sqrt(3(2j+1)); the table already holds dV and SR
      const auto av =
        method == "MS" ?
          ac + alpha_ms_valence(Fv, ww) :
          ac + Amplitudes::sos_valence(0, Fv, Fv, &he1, &he1, ww, -ww, spectrum,
                                       nullptr, nullptr, t_me, s_me) /
                 std::sqrt(3.0 * (Fv.twoj() + 1));
      avs.at(iv) = av;
      if (do_tensor) {
        // alpha_2 = factor * A^2 (zero automatically for j < 1)
        const auto twoJ = double(Fv.twoj());
        const auto factor2 =
          twoJ < 2.0 ?
            0.0 :
            -std::sqrt(2.0 * twoJ * (twoJ - 1.0) /
                       (3.0 * (twoJ + 1.0) * (twoJ + 2.0) * (twoJ + 3.0)));
        a2s.at(iv) = factor2 * Amplitudes::sos_valence(
                                 2, Fv, Fv, &he1, &he1, ww, -ww, spectrum,
                                 nullptr, nullptr, t_me, s_me);
      }
    }
    for (auto &av : avs) {
      if (print)
        printf("%9.2e ", av);
      ofile << av << " ";
    }
    if (rpaQ && rpa_omegaQ) {
      if (print)
        printf("[%.0e]", dVE1.last_eps());
      ofile << dVE1.last_eps();
    }
    if (print)
      std::cout << "\n";
    ofile << "\n";

    if (do_tensor) {
      if (print) {
        std::cout << "          " << "          " << " a2(w)    ";
      }
      for (auto &a2 : a2s) {
        if (print)
          printf("%9.2e ", a2);
        o2file << a2 << " ";
      }
      if (rpaQ && rpa_omegaQ)
        o2file << dVE1.last_eps();
      if (print)
        std::cout << "\n";
      o2file << "\n";
    }
  }
}

} // namespace Module
