#include "Modules/isotopeShift.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

#include <gsl/gsl_fit.h>

namespace Module {

void fieldShift(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"", "Calculates field shift: F = d(E)/d(<r^2>) by direct calculation"},
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
            << "By fitting de = F<dr^2> for small delta r\n";

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
      of << rB - r0 << " " << " " << dr2;
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
  for (auto i = 0ul; i < wfB.valence().size(); ++i) {
    const auto &Fv = wfB.valence().at(i);
    auto &data_v = data[i];

    [[maybe_unused]] double c0, c1, cov00, cov01, cov11, sumsq;

    // Fit, without c0
    // https://www.gnu.org/software/gsl/doc/html/lls.html
    gsl_fit_mul(&data_v[0].first, 2, &data_v[0].second, 2, data_v.size(), &c1,
                &cov11, &sumsq);

    std::cout << Fv.symbol() << " "
              << "F = " << c1 << " GHz/fm^2, sd = " << std::sqrt(sumsq) << "\n";
  }
}

//------------------------------------------------------------------------------
void fieldShift2(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"", "Calculates field shift via matrix element, including 2nd-order "
            "Field shift via TDHF method"},
       {"drrms", "Effective shift in r_rms in fm; must be small; answer "
                 "should be independent? [0.0005]"},
       {"rpa", "Include RPA (uses TDHF method)? [true]"},
       {"use_basis", "Use basis (i.e., MBPT) for 2nd order term, rather than "
                     "TDHF? [false]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto drrms = input.get("drrms", 0.0005);
  const auto rpa = input.get("rpa", true);
  const auto use_basis = input.get("use_basis", false);

  const auto nuc1 = wf.nucleus();
  auto nuc2 = nuc1;
  nuc2.set_rrms(nuc1.r_rms() + drrms);

  DiracOperator::F_fs Ffs(wf.grid(), nuc1, nuc2);
  std::cout << "F1 should be in: " << Ffs.units() << "\n";
  std::cout
      << "\n* Note: second-order term unstable, and not properly checked!\n\n";

  ExternalField::TDHF tdhf(&Ffs, wf.vHF());
  if (rpa)
    tdhf.solve_core(0.0);

  std::cout << "\n";

  fmt::print("{:4s} {:13s}   {:13s}\n", "", "F1 (GHz/fm^2)", "F2 (MHz/fm^4)");
  for (const auto &v : wf.valence()) {
    const auto k = Ffs.rme3js(v.twoj(), v.twoj()); // RME to ME
    const auto f1 = k * Ffs.reducedME(v, v);
    const auto dv = k * tdhf.dV(v, v);

    // If including correlations, should have a spectrum for MBPT
    const auto &basis = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

    const auto hFv = Ffs.reduced_rhs(v.kappa(), v) + tdhf.dV_rhs(v.kappa(), v);
    const auto dFv =
        use_basis ? ExternalField::solveMixedState_basis(v, hFv, 0.0, basis) :
                    tdhf.solve_dPsi(v, 0.0, ExternalField::dPsiType::X,
                                    v.kappa(), wf.Sigma());

    // Two Three J's = (-1 * k * k)
    const auto k_so = 2.0 * k * k / PhysConst::Hartree_GHz * 1.0e3;
    const auto f_so = k_so * (Ffs.reducedME(v, dFv) + tdhf.dV(v, dFv));

    // std::cout << v << " " << f1 + dv << " " << f_so << "\n";

    fmt::print("{:4s} {:+13.5e}   {:+13.5e}\n", v.shortSymbol(), f1 + dv, f_so);
  }
}

} // namespace Module
