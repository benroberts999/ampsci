#include "Modules/thirdOrder.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/Sigma2.hpp"
#include "MBPT/StructureRad.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/ostream.hpp"

namespace Module {

void thirdOrderME(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"operator", "e.g., E1, hfs (see ampsci -o for available operators)"},
       {"options{}", "options specific to operator (see ampsci -o 'operator')"},
       {"rpa", "True/false (always uses diagram method) [true]"},
       {"scale",
        "Scaling factors for BO contributions (same as lambda in "
        "Correlations), comma-sepparated list. There must be one factor for "
        "each valence state, in the same order as the valence states."},
       {"omega",
        "Text or number. Freq. for RPA (and freq. dependent operators). Put "
        "'each' to solve at correct frequency for each transition. [0.0]"},
       {"Tderiv", "Include T_deriv term. Only if omega=each, and "
                  "frequency-dependent operator. [true]"},
       {"printBoth", "print <a|h|b> and <b|h|a> [false]"},
       {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
       {"off-diagonal",
        "Calculate off-diagonal matrix elements (if non-zero) [true]"},
       {"Qk_file",
        "true/false/filename - filename for QkTable file. If blank will "
        "not use QkTable; if exists, will read it in; if doesn't exist, will "
        "create it and write to disk. If 'true' will use default filename. "
        "Save time (10x) at cost of memory."},
       {"n_minmax", "list; min,max n for core/excited in internal diagram "
                    "lines: [1,inf]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  std::cout << "Calculating Third-Order matrix elements\n\n";
  IO::ChronoTimer timer("thirdOrderME");

  const auto oper = input.get<std::string>("operator", "");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  const auto h = DiracOperator::generate(oper, h_options, wf);

  // treat hyperfine operator differently: A constants instead of RME
  const bool hf_AB =
      qip::ci_compare(oper, "hfs") || qip::ci_compare(oper, "MLVP");

  const bool diagonal = input.get("diagonal", true);
  const bool off_diagonal = input.get("off-diagonal", true);

  std::cout << "\n"
            << "Matrix Elements - Operator: " << h->name() << "\n";
  if (hf_AB && h->rank() % 2 != 0) {
    std::cout << "Hyperfine A constants (magnetic type), K=" << h->rank()
              << "\n";
  } else if (hf_AB && h->rank() % 2 == 0) {
    std::cout << "Hyperfine B constants (electric type), K=" << h->rank()
              << "\n";
  } else {
    std::cout << "Reduced matrix elements\n";
  }
  std::cout << "Units: " << h->units() << "\n";

  // RPA: only use diagram method
  auto do_rpa = input.get("rpa", true);
  if (wf.core().empty())
    do_rpa = false;
  auto rpa = ExternalField::make_rpa(do_rpa ? "diagram" : "none", h.get(),
                                     wf.vHF(), true, wf.basis(), wf.identity());

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = qip::ci_compare(str_om, "each");
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  const auto do_Tderiv = input.get("Tderiv", true);

  if (h->freqDependantQ()) {
    std::cout << "Frequency-dependent operator; at omega = ";
    if (eachFreqQ)
      std::cout << "each transition frequency\n";
    else
      std::cout << omega << "\n";
  }

  std::map<std::string, double> scale_factors;
  bool do_scale = false;
  const auto lambda = input.get("scale", std::vector<double>{});
  if (!lambda.empty()) {
    do_scale = true;
    std::cout << "Using scaling factors for BO:\n";
    for (std::size_t i = 0; i < wf.valence().size(); ++i) {
      const auto lam = lambda.size() > i ? lambda.at(i) : 1.0;
      std::cout << wf.valence().at(i) << " : " << lam << "\n";
      scale_factors[wf.valence().at(i).shortSymbol()] = lam;
    }
  }

  const auto n_minmax = input.get("n_minmax", std::vector{1});
  const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
  const auto Qk_file_t = input.get("Qk_file", std::string{"false"});
  std::string Qk_file =
      Qk_file_t != "false" ?
          Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t :
          "";

  std::cout
      << "\nIncluding Structure radiation, normalisation of states and BO:\n";
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";
  if (!Qk_file.empty()) {
    std::cout << "Will read/write Qk integrals to file: " << Qk_file << "\n";
  } else {
    std::cout << "Will calculate Qk integrals on-the-fly\n";
  }
  std::cout << std::flush;

  auto sr =
      MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max}, Qk_file);

  // use basis states, not valence states
  std::vector<DiracSpinor> orbs;
  for (const auto &v : wf.valence()) {
    orbs.push_back(v);
  }

  if (h->freqDependantQ() && !eachFreqQ) {
    h->updateFrequency(omega);
  }

  if (rpa && !eachFreqQ) {
    rpa->solve_core(omega);
  }

  std::stringstream summary;

  // Diagonal elements:
  if (diagonal && h->parity() == 1) {

    if (eachFreqQ && h->freqDependantQ()) {
      h->updateFrequency(0.0);
    }
    if (eachFreqQ && rpa) {
      rpa->solve_core(0.0);
    }

    for (const auto &a : orbs) {
      const auto &b = a;

      if (h->isZero(a.kappa(), b.kappa()))
        continue;

      const auto ff = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                  h->rank(), a.kappa(), b.kappa()) :
                              1.0;

      const auto h0 = h->reducedME(a, b);
      const auto dv = rpa ? rpa->dV(a, b) : 0.0;

      fmt::print("\n{:4s}\n", a.shortSymbol());
      fmt::print("   h0: {:11.4e} + {:11.4e}\n", ff * h0, ff * dv);

      const auto [tb, dvtb] = sr.srTB(h.get(), a, b, 0.0, rpa.get());
      const auto [c, dvc] = sr.srC(h.get(), a, b, rpa.get());

      fmt::print("   SR: {:11.4e} + {:11.4e}\n", ff * (tb + c),
                 ff * (dvtb + dvc - tb - c));

      const auto [n, dvn] = sr.norm(h.get(), a, b, rpa.get());

      fmt::print("   Nm: {:11.4e} + {:11.4e}\n", ff * n, ff * (dvn - n));

      const auto fa = do_scale ? scale_factors[a.shortSymbol()] : 1.0;
      const auto fb = do_scale ? scale_factors[b.shortSymbol()] : 1.0;

      const auto [bo, dvbo] = sr.BO(h.get(), a, b, rpa.get(), fa, fb);

      fmt::print("   BO: {:11.4e} + {:11.4e}\n", ff * bo, ff * (dvbo - bo));

      const auto total = ff * (h0 + dv + dvtb + dvc + dvn + dvbo);
      fmt::print("Total: {:11.4e}\n", total);

      fmt::print(summary,
                 "{:3s} {:3s} {:6.3f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} "
                 "{:11.4e} {:11.4e}\n",
                 a.shortSymbol(), b.shortSymbol(), 0.0, ff * h0, ff * dv,
                 ff * (dvtb + dvc), ff * dvn, ff * dvbo, total);
    }
  }

  // Off-Diagonal elements:
  if (off_diagonal) {
    const auto pi = h->parity();

    for (std::size_t ib = 0; ib < orbs.size(); ib++) {
      const auto &b = orbs.at(ib);
      for (std::size_t ia = 0; ia < orbs.size(); ia++) {
        const auto &a = orbs.at(ia);

        if (a == b)
          continue;
        if (h->isZero(a.kappa(), b.kappa()))
          continue;

        // Ensure even-parity state on right for odd-parity operators
        if (pi == -1) {
          if (b.parity() == -1)
            continue;
        } else {
          if (ib > ia)
            continue;
        }

        const auto ww = eachFreqQ ? a.en() - b.en() : omega;

        fmt::print("\n{:4s} - {:4s}: omega = {:.6f}\n", a.shortSymbol(),
                   b.shortSymbol(), ww);

        if (eachFreqQ && h->freqDependantQ()) {
          h->updateFrequency(ww);
        }
        if (eachFreqQ && rpa) {
          rpa->solve_core(ww);
        }

        if (h->isZero(a.kappa(), b.kappa()))
          continue;

        const auto ff = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                    h->rank(), a.kappa(), b.kappa()) :
                                1.0;

        const auto h0 = h->reducedME(a, b);
        const auto dv = rpa ? rpa->dV(a, b) : 0.0;
        fmt::print("   h0: {:11.4e} + {:11.4e}\n", ff * h0, ff * dv);

        const auto [tb, dvtb] = sr.srTB(h.get(), a, b, ww, rpa.get());
        const auto [c, dvc] = sr.srC(h.get(), a, b, rpa.get());

        fmt::print("   SR: {:11.4e} + {:11.4e}\n", ff * (tb + c),
                   ff * (dvtb + dvc - tb - c));

        const auto [n, dvn] = sr.norm(h.get(), a, b, rpa.get());

        fmt::print("   Nm: {:11.4e} + {:11.4e}\n", ff * n, ff * (dvn - n));

        // Calculate T_deriv:
        double T_deriv = 0.0;
        if (eachFreqQ && h->freqDependantQ() && do_Tderiv) {
          std::cout << "   T_deriv (added to BO):\n";
          const auto de2_a =
              MBPT::Sigma_vw(a, a, sr.Yk(), sr.core(), sr.excited());
          const auto de2_b =
              MBPT::Sigma_vw(b, b, sr.Yk(), sr.core(), sr.excited());
          const auto dw = de2_a - de2_b;
          std::cout << "        da(2): " << de2_a << "\n";
          std::cout << "        db(2): " << de2_b << "\n";
          std::cout << "        dw(2): " << dw << "\n";

          h->updateFrequency(ww + 0.005);
          const auto h_plus = h->reducedME(a, b);
          h->updateFrequency(ww - 0.005);
          const auto h_minus = h->reducedME(a, b);
          const auto d_h = (h_plus - h_minus) / 0.01;
          std::cout << "        dT/dw: " << ff * d_h << "\n";
          std::cout << "        T_dw : " << ff * d_h * dw << "\n";
          T_deriv = d_h * dw;
          h->updateFrequency(ww);
        }

        const auto fa = do_scale ? scale_factors[a.shortSymbol()] : 1.0;
        const auto fb = do_scale ? scale_factors[b.shortSymbol()] : 1.0;

        auto [bo, dvbo] = sr.BO(h.get(), a, b, rpa.get(), fa, fb);
        bo += T_deriv;
        dvbo += T_deriv;

        fmt::print("   BO: {:11.4e} + {:11.4e}\n", ff * bo, ff * (dvbo - bo));

        const auto total = ff * (h0 + dv + dvtb + dvc + dvn + dvbo);

        fmt::print("Total: {:11.4e}\n", total);

        fmt::print(summary,
                   "{:3s} {:3s} {:6.3f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} "
                   "{:11.4e} {:11.4e}\n",
                   a.shortSymbol(), b.shortSymbol(), ww, ff * h0, ff * dv,
                   ff * (dvtb + dvc), ff * dvn, ff * dvbo, total);
      }
    }
  }

  std::cout << "\n" << h->name() << " (" << h->units() << ")\n";
  fmt::print("{:3s} {:3s} {:>6s} {:>11s} {:>11s} {:>11s} {:>11s} "
             "{:>11s} {:>11s}\n",
             "w", "v", "omega", "t0", "RPA", "SR", "Norm", "BO", "Total");
  std::cout << summary.str() << "\n";
}

} // namespace Module
