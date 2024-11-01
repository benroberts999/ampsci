#include "Modules/isotopeShift.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/ostream.hpp"
#include "qip/Array.hpp"
#include "qip/Vector.hpp"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>

#include <iostream>

namespace Module {

void isotopeShift(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({{"", "Calculates field shift using MBPT, including 2nd-order "
                    "(quadratic) field shift, via the TDHF/RPA(basis) method"},
               {"RPA", "Include RPA? true/Ffalse [true]"},
               {"num_iter", "Number of iterations to compute changes for [20]"},
               {"dr", "Smallest change in rms charge radius (fm) [0.0001]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  using namespace qip::overloads;

  IO::ChronoTimer timer("isotopeShift");

  const auto num_iter = std::max(input.get("num_iter", 20ul), 2ul);
  const auto delta_r = input.get("dr", 0.0001);
  const auto rpaQ = input.get("RPA", true);

  // Initial (reference) nuclear parameters:
  const auto &nuc0 = wf.nucleus();
  const auto r_rms0 = nuc0.r_rms();
  const auto r20 = r_rms0 * r_rms0;
  const auto r40 = DiracOperator::fieldshift::r4(wf.grid(), wf.nucleus());

  std::cout
      << "\nCalculating shift in valence state energies using.\n"
      << "dE^(1) = F d<r^2> + G^2 d<r^2>^2 + G^4 d<r^4>\n"
      << "[nb: E is binding energy (different sign from someother works)]\n";

  // For second-order G^2 correction,
  // See, e.g., Eq.(8) PhysRevA.103.L030801 (2021)

  std::cout << "\nInitial nuclear parameters:\n"
            << "r_rms_0 = " << r_rms0 << " fm, "
            << "<r^2>_0 = " << r20 << " fm^2, "
            << "<r^4>_0 = " << r40 << " fm^4\n";

  // Grid for delta_rrms:
  const auto gtype = GridType::logarithmic;
  // const auto gtype = GridType::linear;
  const auto range = gtype == GridType::logarithmic ?
                         qip::logarithmic_range<double, std::size_t> :
                         qip::uniform_range<double, std::size_t>;

  const auto max_dr = delta_r * int(num_iter) / 2.0;
  const auto drs_t = range(delta_r, max_dr, num_iter / 2);
  const auto drs = qip::merge(-1.0 * qip::reverse(drs_t), drs_t);

  // Arrays to store data for fitting:
  qip::Array F1_data(wf.valence().size(), num_iter);
  qip::Array G2_data(wf.valence().size(), num_iter);

  std::cout << "\ndE(1) (GHz)\n";
  fmt::print("{:6s} {:8s} {:8s}", "r_rms", " d<r2>", " d<r4>");
  for (const auto &Fv : wf.valence()) {
    fmt::print("  {:10s}", Fv.shortSymbol());
  }
  std::cout << std::endl;
  std::stringstream gout, gout2;

  // Loop over nuclear charge radii:
  auto nuc = nuc0;
  for (std::size_t i = 0ul; i < drs.size(); ++i) {
    if (std::abs(drs.at(i)) < 1.0e-10)
      continue;

    nuc.set_rrms(r_rms0 + drs.at(i));

    DiracOperator::fieldshift fis(wf.grid(), nuc0, nuc);

    const auto r_rms = nuc.r_rms();
    const auto delta_r2 = fis.dr2();
    const auto delta_r4 = fis.dr4();

    // Always use TDHF now
    ExternalField::TDHF dVfis(&fis, wf.vHF());

    if (rpaQ)
      dVfis.solve_core(0.0, 100, false);

    fmt::print("{:6.4f} {:8.5f} {:8.5f}", r_rms, delta_r2, delta_r4);
    fmt::print(gout, "{:6.4f} {:8.5f} {:8.2e}", r_rms, delta_r2,
               delta_r2 * delta_r2);
    fmt::print(gout2, "{:6.4f} {:8.5f} {:8.2e}", r_rms, delta_r2,
               delta_r2 * delta_r2);

    //-------------------------------------------------------
    // Loop over valence states:
    for (std::size_t j = 0ul; j < wf.valence().size(); ++j) {
      const auto &Fv = wf.valence().at(j);

      const auto k = fis.rme3js(Fv.twoj(), Fv.twoj()); // RME -> ME

      // First-order correction:
      auto F_rme = fis.reducedME(Fv, Fv) + dVfis.dV(Fv, Fv);
      const auto dE = F_rme * k * PhysConst::Hartree_GHz;

      fmt::print(" {:11.4e}", dE);

      F1_data(j, i) = dE;

      // Second-order G^2 correction:
      {

        // Use TDHF method:
        const auto hFv =
            fis.reduced_rhs(Fv.kappa(), Fv) + dVfis.dV_rhs(Fv.kappa(), Fv);

        const auto dFv1 = dVfis.solve_dPsi(Fv, 0.0, ExternalField::dPsiType::X,
                                           Fv.kappa(), wf.Sigma());

        // .. and using SOS method:
        // If including correlations, should have a spectrum for MBPT
        const auto &basis = wf.spectrum().empty() ? wf.basis() : wf.spectrum();
        // Note: negative energy states are required here!!
        const auto dFv2 =
            ExternalField::solveMixedState_basis(Fv, hFv, 0.0, basis);

        const auto G2_rme = fis.reducedME(Fv, dFv1) + dVfis.dV(Fv, dFv1);
        const auto G2_rme2 = fis.reducedME(Fv, dFv2) + dVfis.dV(Fv, dFv2);
        const auto dE2 = G2_rme * k * k * PhysConst::Hartree_GHz;

        // Just ude TDHF for the fit: more stable
        G2_data(j, i) = dE2;

        fmt::print(gout, " {:11.4e}", G2_rme);
        fmt::print(gout2, " {:11.4e}", G2_rme2);
      }
    }
    if (dVfis.last_eps() > 1.0e-8)
      std::cout << " ***";
    std::cout << std::endl;
    fmt::print(gout, "\n");
    fmt::print(gout2, "\n");
  }

  // Print the G2 results below:
  {
    std::cout << "\ndE(2) (GHz) - TDHF method\n";
    fmt::print("{:6s} {:8s} {:8s}", "r_rms", " d<r2>", " d<r2>^2");
    for (const auto &Fv : wf.valence()) {
      fmt::print("  {:10s}", Fv.shortSymbol());
    }
    std::cout << std::endl;
    std::cout << gout.str();

    std::cout << "\ndE(2) (GHz) - SOS method\n"
                 " * note: Requires negative energy states!\n";
    fmt::print("{:6s} {:8s} {:8s}", "r_rms", " d<r2>", " d<r2>^2");
    for (const auto &Fv : wf.valence()) {
      fmt::print("  {:10s}", Fv.shortSymbol());
    }
    std::cout << std::endl;
    std::cout << gout2.str();
  }

  // Used for the linear fits:
  const auto dr2s = (r_rms0 + drs) * (r_rms0 + drs) - r20;
  const auto dr2s2 = dr2s * dr2s;

  std::cout << "\nv     F (GHz/fm^2)   err*       G2 (GHz/fm^4)   err\n";

  // For the transitions:
  std::stringstream out_trans;
  std::optional<double> F0, G20;
  std::string v0;
  for (auto i = 0ul; i < wf.valence().size(); ++i) {
    const auto &Fv = wf.valence().at(i);
    auto G2_data_v = G2_data.row(i);
    auto F1_data_v = F1_data.row(i);

    // Linear fit, use GSL
    // Not strictly necisary
    // https://www.gnu.org/software/gsl/doc/html/lls.html

    double F{0.0}, G{0.0}, covF{0.0}, covG{0.0}, sumsq{0.0};

    gsl_fit_mul(dr2s.data(), 1, F1_data_v.data(), 1, dr2s.size(), &F, &covF,
                &sumsq);

    gsl_fit_mul(dr2s2.data(), 1, G2_data_v.data(), 1, dr2s2.size(), &G, &covG,
                &sumsq);

    fmt::print("{:4s} {:13.6e}   {:7.1e}   {:13.6e}    {:7.1e}\n",
               Fv.shortSymbol(), F, std::sqrt(covF), G, std::sqrt(covG));

    // Store first F0, G0
    if (!F0) {
      F0 = F;
      G20 = G;
      v0 = Fv.shortSymbol();
    }
    fmt::print(out_trans, "{:4s} - {:4s} {:13.6e}   {:13.6e}\n",
               Fv.shortSymbol(), v0, F - *F0, G - *G20);
  }
  std::cout << "\n";
  std::cout << "(* err is just from the fit: "
            << "it is a minimum numerical error only)\n";

  std::cout << "\nTransitions:\n";
  std::cout << "           F (GHz/fm^2)   G2 (GHz/fm^4)\n";
  std::cout << out_trans.str() << "\n";
}

//==============================================================================
void fieldShift(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"", "Calculates field shift: F = d(E)/d(<r^2>) by direct calculation. "
            "Note: copies the same correlation potential; this is OK, but not "
            "exact (i.e., neglects the SR contribution)"},
       {"core_relaxation", "Include Core relaxation (equiv to RPA)? [true]"},
       {"print", "Print each step to screen? [true]"},
       {"write", "Write dE(r^2) to file? [false]"},
       {"minmax_delta", "Minimum relative shift in r [1.0e-5, 1.0e-3]"},
       {"num_steps", "Number of steps for fit (for each sign)? [5]"},
       {"grid", "Logarithmic or linear grid for dr2 [logarithmic]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto core_relax = input.get("core_relaxation", true);
  const auto print = input.get("print", true);
  const auto write = input.get("write", false);

  const auto [min_d, max_d] =
      input.get("minmax_delta", std::array{1.0e-5, 1.0e-3});

  const auto num_steps = input.get<unsigned long>("num_steps", 5);

  const auto grid_type = input.get("grid", std::string{"logarithmic"});

  Wavefunction wfB(wf.grid_sptr(), wf.nucleus(), wf.alpha() / PhysConst::alpha);

  std::cout << "Calculating field shift corrections for \n"
            << wf.atom() << ", " << wf.nucleus() << "\n"
            << "By fitting de = F<dr^2> for small delta r\n"
            << "Directly re-solves Hartree-Fock at each step\n"
            << "(Note: does not re-calculate Sigma!)\n";

  wfB.copySigma(wf.Sigma());
  const auto core_string = wf.coreConfiguration();
  const auto val_string = DiracSpinor::state_config(wf.valence());
  const auto r0 = wf.get_rrms();

  std::cout << "\n";
  wfB.solve_core("HartreeFock", 0.0, core_string, 0.0, true);
  if (!core_relax) {
    std::cout << "Not including Core relaxtion\n";
  } else {
    std::cout << "Including Core relaxtion\n";
  }

  // Optionally write dE(R) to file for plotting
  std::ofstream of;
  if (write) {
    of.open(wf.identity() + "_FS.txt");
    of << "dr dr2 ";
    for (auto &v : wf.valence()) {
      of << " " << v;
    }
    of << "\n";
  }

  const auto gtype = qip::ci_wc_compare(grid_type, "log*") ?
                         GridType::logarithmic :
                         GridType::linear;

  const auto drs_t =
      gtype == GridType::logarithmic ?
          qip::logarithmic_range(r0 * min_d, r0 * max_d, num_steps) :
          qip::uniform_range(r0 * min_d, r0 * max_d, num_steps);
  using namespace qip::overloads;
  const auto drs = qip::merge(-1.0 * qip::reverse(drs_t), drs_t);

  if (print) {
    std::cout << "\n   r_rms (fm)     del(r)     del(r^2)     dE (GHz)   F "
                 "(GHz/fm^2)\n";
  } else {
    std::cout << "\nRunning...\n";
  }

  // Store data - used for linear fit to get F
  std::vector<std::vector<std::pair<double, double>>> data(wf.valence().size());

  for (const auto del : drs) {
    const auto rB = r0 + del;
    const auto dr2 = rB * rB - r0 * r0;

    auto nuc_b = wf.nucleus();
    nuc_b.set_rrms(rB);

    wfB.update_Vnuc(Nuclear::formPotential(nuc_b, wf.grid().r()));

    if (write) {
      of << rB - r0 << " " << dr2;
    }

    if (core_relax)
      wfB.solve_core("HartreeFock", 0.0, core_string, 0.0, false);
    wfB.solve_valence(val_string, false);
    wfB.hartreeFockBrueckner(false);

    for (auto i = 0ul; i < wfB.valence().size(); ++i) {
      const auto &Fv = wfB.valence().at(i);
      const auto &Fv0 = *wf.getState(Fv.n(), Fv.kappa());
      const auto dE = (Fv.en() - Fv0.en()) * PhysConst::Hartree_GHz;
      const auto tF = dE / dr2;
      if (print)
        printf("%4s  %7.5f  %+8.6f  %11.4e  %11.4e  %10.3e\n",
               Fv.shortSymbol().c_str(), rB, rB - r0, dr2, dE, tF);
      auto &data_v = data[i];
      data_v.emplace_back(dr2, dE);

      if (write) {
        of << " " << dE;
      }
    }

    if (write) {
      of << "\n";
    }
    if (print) {
      std::cout << "\n";
    }
  }

  std::cout << "\n";

  // Fit straight line to data
  std::cout << "\nv     F (GHz/fm^2)   err*\n";
  for (auto i = 0ul; i < wfB.valence().size(); ++i) {
    const auto &Fv = wfB.valence().at(i);
    auto &data_v = data[i];

    [[maybe_unused]] double c0, c1, cov00, cov01, cov11, sumsq;

    // Fit, without c0
    // https://www.gnu.org/software/gsl/doc/html/lls.html
    gsl_fit_mul(&data_v[0].first, 2, &data_v[0].second, 2, data_v.size(), &c1,
                &cov11, &sumsq);

    // std::cout << Fv.symbol() << " "
    //           << "F = " << c1 << " GHz/fm^2, sd = " << std::sqrt(sumsq) << "\n";

    fmt::print("{:4s} {:13.6e}   {:7.1e}\n", Fv.shortSymbol(), c1,
               std::sqrt(cov11));
  }
  std::cout << "\n";
  std::cout << "(* err is just from the fit: "
            << "it is a minimum numerical error only)\n";
}

} // namespace Module