#include "Modules/HFAnomaly.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/matrixElements.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/Methods.hpp"
#include <fstream>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <iostream>
#include <string>
#include <vector>

namespace Module {

//==============================================================================
void BW_effect(const std::vector<DiracSpinor> &valence,
               const DiracOperator::TensorOperator *const h0,
               const ExternalField::DiagramRPA *const rpa0,
               const DiracOperator::TensorOperator *const h,
               const ExternalField::DiagramRPA *const rpa) {

  fmt::print("\n{:4s} {:11s} {:11s}  {:6s}", "", "A0_HF (MHz)", "A_HF",
             "eBW(%)");
  fmt::print("  {:11s} {:11s}  {:6s}", "A0_RPA", "A_RPA", "eBW(%)");
  std::cout << "\n";
  for (const auto &Fv : valence) {
    const auto rme_to_A = DiracOperator::Hyperfine::convert_RME_to_AB(
      h0->rank(), Fv.kappa(), Fv.kappa());
    const auto A0 = h0->reducedME(Fv, Fv) * rme_to_A;
    const auto dA0 = rpa0 ? rpa0->dV(Fv, Fv) * rme_to_A : 0.0;
    const auto A = h->reducedME(Fv, Fv) * rme_to_A;
    const auto dA = rpa ? rpa->dV(Fv, Fv) * rme_to_A : 0.0;

    const auto eps_HF = 100.0 * (A - A0) / A0;
    const auto eps_RPA = 100.0 * (A + dA - A0 - dA0) / (A0 + dA0);

    fmt::print("{:4s} {:.5e} {:.5e} {:7.4f}", Fv.shortSymbol(), A0, A, eps_HF);
    if (rpa)
      fmt::print("  {:.5e} {:.5e} {:7.4f}", A0 + dA0, A + dA, eps_RPA);
    std::cout << "\n";
  }
}

//==============================================================================
void tune_Rmag(const DiracSpinor &Fv, const double eps_target,
               const std::optional<IO::InputBlock> &hfs_options,
               const Wavefunction &wf,
               const DiracOperator::TensorOperator *const h0,
               const ExternalField::DiagramRPA *const rpa0) {

  std::cout << "State: " << Fv << ", BW target: " << eps_target << "%\n";
  if (rpa0) {
    std::cout << "Including RPA\n";
  }

  // initial guess for magnetic radius: charge radius
  double t_rrms0 = wf.get_rrms();

  std::unique_ptr<ExternalField::DiagramRPA> rpa_t{nullptr};
  if (rpa0) {
    rpa_t = std::make_unique<ExternalField::DiagramRPA>(h0, rpa0);
  }

  // Define function to be minimised:
  const auto delta_eps = [&](double t_rrms) {
    // const auto &Fv = *pFv;
    auto t_hfs_options = hfs_options ? *hfs_options : IO::InputBlock{};
    // t_hfs_options.add("rrms=0.1;");
    t_hfs_options.add(IO::Option{"rrms", std::to_string(t_rrms)});
    t_hfs_options.add(IO::Option{"print", "false"});
    const auto ht = DiracOperator::generate("hfs", t_hfs_options, wf);

    if (rpa_t) {
      rpa_t->update_t0s(ht.get());
      rpa_t->solve_core(0.0, 50, false);
      if (rpa_t->last_eps() > 1.0e-4) {
        std::cout << "* WARNING: RPA didn't converge?: " << rpa_t->last_eps()
                  << "\n";
      }
    }

    const auto rme_to_A = DiracOperator::Hyperfine::convert_RME_to_AB(
      h0->rank(), Fv.kappa(), Fv.kappa());
    const auto A =
      (ht->reducedME(Fv, Fv) + (rpa_t ? rpa_t->dV(Fv, Fv) : 0.0)) * rme_to_A;
    const auto A0 =
      (h0->reducedME(Fv, Fv) + (rpa0 ? rpa0->dV(Fv, Fv) : 0.0)) * rme_to_A;

    const auto eps = 100.0 * (A - A0) / A0;
    const auto delta = eps - eps_target;
    fmt::print("{:.4f} {:+.4f} [{:+.4f}]\n", t_rrms, eps, std::abs(delta));
    return delta;
  };

  std::cout << "Use Newton's method:\n";
  std::cout << "Rm_rms eps(%)  [delta]\n";
  // Use Newtons method to do fit"
  const auto [rrms_fitted, error] = qip::Newtons(
    delta_eps, t_rrms0, {0.75 * t_rrms0, 1.5 * t_rrms0}, 1.0e-3, 0.05, 100);

  // to determine error function: d(eps)/d(rrms):
  const auto err_factor = 1.002;
  const auto del_eps = delta_eps(rrms_fitted * err_factor);
  const auto del_rms = rrms_fitted * err_factor - rrms_fitted;
  const auto dRrms_deps = del_rms / del_eps;

  std::cout << "\nFitted Rmag_rms = " << rrms_fitted
            << " fm  [fit error: " << error << "]\n";
  std::cout << "Fitted Rmag     = " << std::sqrt(5.0 / 3) * rrms_fitted
            << " fm\n";

  std::cout << "\nError function, d(Rrms)/d(eps) = " << dRrms_deps << " fm/%\n";
}

//==============================================================================
void BW_screening_factor(const Wavefunction &wf,
                         const DiracOperator::TensorOperator *const h0,
                         const ExternalField::DiagramRPA *const rpa0,
                         const DiracOperator::TensorOperator *const h,
                         const ExternalField::DiagramRPA *const rpa) {

  std::cout << "\n-------------------------------\n";
  std::cout << "Screening factors, as per PhysRevA.105.052802:\n";
  std::cout << "\nCalculate H-like wavefunctions (for screening):\n";
  const auto grid_params = wf.grid().params();
  auto wfH =
    Wavefunction(grid_params, wf.nucleus(), wf.alpha() / PhysConst::alpha);
  std::cout << wfH.grid().gridParameters() << "\n";
  wfH.solve_valence("1s2p");
  wfH.printValence();

  fmt::print("\n{:3s} {:11s} {:7s} {:4s} {:11s} {:7s} {:7s}", "", "H-like",
             "eps(%)", "", "HF", "eps(%)", " x(HF)");
  fmt::print("  {:11s} {:7s} {:7s}", "RPA", "eps(%)", " x(RPA)");
  std::cout << "\n";

  for (const auto &Fv : wf.valence()) {
    if (Fv.twoj() > 1)
      continue;

    const auto rme_to_A = DiracOperator::Hyperfine::convert_RME_to_AB(
      h0->rank(), Fv.kappa(), Fv.kappa());

    const auto A0 = h0->reducedME(Fv, Fv) * rme_to_A;
    const auto A = h->reducedME(Fv, Fv) * rme_to_A;

    const auto dA0 = rpa0 ? rpa0->dV(Fv, Fv) * rme_to_A : 0.0;
    const auto dA = rpa ? rpa->dV(Fv, Fv) * rme_to_A : 0.0;

    const auto eps_hf = 100.0 * (A - A0) / A0;
    const auto eps_rpa = 100.0 * (A + dA - A0 - dA0) / (A0 + dA0);

    const auto F1s = *wfH.getState(Fv.l() + 1, Fv.kappa());
    const auto A01s = h0->reducedME(F1s, F1s) * rme_to_A;
    const auto A1s = h->reducedME(F1s, F1s) * rme_to_A;
    const auto eps_1s = 100.0 * (A1s - A01s) / A01s;

    const auto x_HF = eps_hf / eps_1s;
    const auto x_RPA = eps_rpa / eps_1s;

    fmt::print("{:3s} {:.5e} {:7.4f} {:4s} {:.5e} {:7.4f} {:7.4f}",
               F1s.shortSymbol(), A01s, eps_1s, Fv.shortSymbol(), A0, eps_hf,
               x_HF);
    if (rpa0)
      fmt::print("  {:.5e} {:7.4f} {:7.4f}", A, eps_rpa, x_RPA);
    std::cout << "\n";
  }

  const auto F1s = *wfH.getState(1, -1);
  const auto F2p = *wfH.getState(2, 1);

  const auto rme_to_A_s =
    DiracOperator::Hyperfine::convert_RME_to_AB(h0->rank(), -1, -1);
  const auto rme_to_A_p =
    DiracOperator::Hyperfine::convert_RME_to_AB(h0->rank(), 1, 1);

  const auto A01s = h0->reducedME(F1s, F1s) * rme_to_A_s;
  const auto A1s = h->reducedME(F1s, F1s) * rme_to_A_s;
  const auto eps_1s = 100.0 * (A1s - A01s) / A01s;

  const auto A02p = h0->reducedME(F2p, F2p) * rme_to_A_p;
  const auto A2p = h->reducedME(F2p, F2p) * rme_to_A_p;
  const auto eps_2p = 100.0 * (A2p - A02p) / A02p;

  const auto eta_sp = eps_1s / eps_2p;
  std::cout << "eta_sp = " << eta_sp << "\n";
}

//==============================================================================
std::array<double, 3> fit(std::vector<double> &xd, std::vector<double> &yd,
                          std::size_t terms) {
  std::array<double, 3> out;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;
  const auto n = xd.size();
  X = gsl_matrix_alloc(n, terms);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);

  c = gsl_vector_alloc(terms);
  cov = gsl_matrix_alloc(terms, terms);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < terms; ++j) {
      gsl_matrix_set(X, i, j, std::pow(xd.at(i), 2 * (j + 1)));
    }

    gsl_vector_set(y, i, yd.at(i));
    gsl_vector_set(w, i, 1.0);
  }

  double chisq;
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, terms);
  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
  gsl_multifit_linear_free(work);

  auto C = [&](std::size_t i) { return gsl_vector_get(c, i); };

  out[0] = C(0);
  out[1] = C(1);
  out[2] = C(2);

  // fmt::print("K{}(Rm) = ", 2 * terms);
  // for (std::size_t j = 0; j < terms; ++j) {
  //   const auto p = 2 * (j + 1);
  //   fmt::print("{:+.3e} *Rm**{}", C(j), p);
  //   if (j + 1 != terms)
  //     std::cout << " ";
  // }
  // std::cout << "\n";
  // std::cout << "#         [";
  fmt::print("k{} = [", 2 * terms);
  for (std::size_t j = 0; j < terms; ++j) {
    const auto p = 2 * (j + 1);
    fmt::print("{:.4e}", C(j), p);
    if (j + 1 != terms)
      std::cout << ", ";
  }
  std::cout << "]\n";

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  return out;
}

//==============================================================================
std::pair<std::array<double, 3>, std::array<double, 3>>
b_moments(const std::string &iso, const DiracSpinor &v, double R0_fm,
          int max_power) {

  std::cout << "\nb moments: " << v << "\n";

  const auto &grid = v.grid();

  const auto R0 = R0_fm / PhysConst::aB_fm;
  const auto i0 = grid.getIndex(R0);
  std::cout << "Fit up to: Rmax = " << R0_fm << " fm  --  ";
  std::cout << "[points within Rmax: " << i0 << "]\n";

  // Why?? Just use atomig grid!
  // const auto ri = 0.001 * R0;
  // const auto rf = 1.1 * R0;
  // const int num_steps = 10000;
  // const bool log = true;
  // const auto R_pts = log ? qip::logarithmic_range(ri, rf, num_steps) :
  //                          qip::uniform_range(ri, rf, num_steps);

  const auto r3 = grid.rpow(3.0);
  const auto r2inv = grid.rpow(-2.0);

  // std::cout << ri * PhysConst::aB_fm << " - " << rf * PhysConst::aB_fm << " in "
  //           << num_steps << " (" << (log ? "logarithmic" : "linear")
  //           << ") steps\n";

  const auto fg_inf =
    NumCalc::integrate(1.0, 0, v.max_pt(), v.f(), v.g(), r2inv, grid.drdu());

  const auto fg_rm = [&v, &grid, &r2inv](std::size_t i_rm) {
    // const auto im = grid.getIndex(rm, true);
    // if (im < 20) {
    //   std::cout << "Warning: not enough points within nuclear radius: " << rm
    //             << "\n";
    // }
    return NumCalc::integrate(1.0, 0, i_rm, v.f(), v.g(), r2inv, grid.drdu());
  };

  const auto fg_r3_rm = [&v, &grid, &r3, &r2inv](std::size_t i_rm) {
    // const auto im = grid.getIndex(rm, true);
    return NumCalc::integrate(1.0, 0, i_rm, v.f(), v.g(), r2inv, r3,
                              grid.drdu());
  };

  std::ofstream of;
  std::ofstream of2;
  if (!iso.empty()) {
    // Only write output if input string given
    of.open("KS_" + iso + "_" + v.shortSymbol() + ".txt");
    of2.open("KL_" + iso + "_" + v.shortSymbol() + ".txt");
  }

  std::vector<double> xd;
  std::vector<double> yd;
  std::vector<double> zd;

  for (std::size_t i_rm = 10; i_rm < i0; ++i_rm) {
    const auto rm = grid.r(i_rm);
    const auto rm3 = rm * rm * rm;
    const auto F_rm = fg_rm(i_rm) / fg_inf;
    const auto F_r3_rm = (fg_rm(i_rm) - fg_r3_rm(i_rm) / rm3) / fg_inf;
    of << rm * PhysConst::aB_fm << " " << F_rm << "\n";
    of2 << rm * PhysConst::aB_fm << " " << F_r3_rm << "\n";
    xd.push_back(rm * PhysConst::aB_fm);
    yd.push_back(F_rm);
    zd.push_back(F_r3_rm);
  }

  const int max_terms = max_power / 2;

  const auto Ks = fit(xd, yd, std::size_t(max_terms));
  const auto Kl = fit(xd, zd, std::size_t(max_terms));

  return {Ks, Kl};

  // std::cout << "\nKS(R):\n";
  // for (int i = 1; i <= max_terms; ++i)
  //   fit(xd, yd, std::size_t(i));

  // std::cout << "\nKL(R):\n";
  // for (int i = 1; i <= max_terms; ++i)
  //   fit(xd, zd, std::size_t(i));
}

//==============================================================================
//==============================================================================
void HFAnomaly(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"hfs_options{}", "Options for HFS operator (see -o hfs)"},
     {"rpa", "Include RPA (diagram method)? [false]"},
     {"screening", "Calculate screening parameters (x and eta) [false]"},
     {"b_power", "Maximum power for which to caculate b moments in "
                 "KS(R),KL(R) expansion. Must be >=2 to calculate any [0]"},
     {"b_max_ratio_Rnuc",
      "Maximum radius to include in fit for b (as ratio of R_nuc) [1]"},
     {"eps_target",
      "Tune magnetic radius to match experimental BW effect "
      "(eps). Two inputs, comma separated: state (in 'short "
      "symbol' form), and eps (in %). E.g.: '6s+, -0.05' [optional]"},
     {"A2", "Second isotope (for differential anomaly) [optional]"},
     {"Nucleus2{}", "Nuclear (charge) parameters for isotope 2 (see -a "
                    "Nucleus); uses default for A2 if blank. [optional]"},
     {"hfs2_options{}", "Options for HFS operator for isotope 2 (see -o hfs)"},
     {"1D2_target",
      "Fits magnetic radius for both nuclei to reproduce given target for "
      "differential hyperfine anomaly (1D2). Two inputs, comma separated: "
      "state (in 'short symbol' form), and 1D2 (in %). E.g.: '6s+, 0.5' "
      "[optional]"},
     {"min_max_steps",
      "List: For 1D2_target only: minimum and maximum magnetic "
      "radii, as a fraction of charge radius, and number of "
      "steps in the fit. Default=[0.9,1.5,10]"}});
  if (input.has_option("help")) {
    return;
  }

  const auto hfs_options = input.getBlock("hfs_options");

  // Generate hyperfine operator (including BW effect from input)
  const auto h = DiracOperator::generate(
    "hfs", hfs_options ? *hfs_options : IO::InputBlock{}, wf);

  // Build pointlike HFS operator
  IO::InputBlock h0_options{"h0", "F = pointlike;"};
  // Ensure same parameters as other operator
  if (hfs_options) {
    if (hfs_options->has_option("mu"))
      h0_options.add(*hfs_options->getOption("mu"));
    if (hfs_options->has_option("I"))
      h0_options.add(*hfs_options->getOption("I"));
    if (hfs_options->has_option("Q"))
      h0_options.add(*hfs_options->getOption("Q"));
    if (hfs_options->has_option("k"))
      h0_options.add(*hfs_options->getOption("k"));
    if (hfs_options->has_option("print"))
      h0_options.add(*hfs_options->getOption("print"));
  }
  const auto h0 = DiracOperator::generate("hfs", h0_options, wf);

  // RPA:
  const auto rpaQ = input.get("rpa", false);
  std::unique_ptr<ExternalField::DiagramRPA> rpa0{nullptr}, rpa{nullptr};
  if (rpaQ) {
    std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
    rpa0 = std::make_unique<ExternalField::DiagramRPA>(h0.get(), wf.basis(),
                                                       wf.vHF(), wf.identity());
    rpa = std::make_unique<ExternalField::DiagramRPA>(h.get(), rpa0.get());

    std::cout << "Solving RPA:\n"
              << "For pointlike:\n";
    rpa0->solve_core(0.0, 100, true);
    std::cout << "For finite-magnetisation:\n";
    rpa->solve_core(0.0, 100, true);
  }

  std::cout << "\n-------------------------------\n";
  fmt::print("Bohr-Weisskopf effect, Z={}, A={}\n", wf.Znuc(), wf.Anuc());
  fmt::print("R_charge = {:.4f} fm [Rrms_charge = {:.4f} fm]\n",
             std::sqrt(5.0 / 3) * wf.get_rrms(), wf.get_rrms());

  // 1. Calculate BW effect
  BW_effect(wf.valence(), h0.get(), rpa0.get(), h.get(), rpa.get());

  // calculate the B moments in KS(R) and KL(R) expansions
  // maximum moment to consider:
  const int b_power = input.get("b_power", 0);
  // Maximum radius to include in fit (as ratio of R_nuc)
  const double max_ratio = input.get("b_max_ratio_Rnuc", 1.0);
  if (b_power >= 2) {
    const auto iso = wf.identity() + "-" + std::to_string(wf.Anuc());

    for (const auto &v : wf.valence()) {
      b_moments(iso, v, std::sqrt(5.0 / 3) * wf.get_rrms() * max_ratio,
                b_power);
    }
  }

  //----------------------------------------------------------------------------
  // 1.5: Tune magnetic radius to reproduce experimental BW effect:
  if (input.get("eps_target")) {
    std::cout << "\n-------------------------------\n";
    std::cout << "Tuning Rmag to reproduce Bohr-Weiskopf effect:\n";
    const auto [a, eps_s] =
      *input.get<std::array<std::string, 2>>("eps_target");
    // fails if eps_s invalid 'double':
    const auto eps_target = eps_s == "" ? 0.0 : std::stod(eps_s);
    const auto pFv = wf.getState(a);
    if (!pFv) {
      std::cout << "Invalid state: " << a << "\n";
    } else {
      tune_Rmag(*pFv, eps_target, hfs_options, wf, h0.get(), rpa0.get());
    }
  }

  //----------------------------------------------------------------------------

  // 2. BW screening
  if (input.get("screening", false)) {
    BW_screening_factor(wf, h0.get(), rpa0.get(), h.get(), rpa.get());
  }

  // no longer need RPA for pointlike; release memory
  rpa0.reset(nullptr);

  //----------------------------------------------------------------------------

  // 3. Differential HF anomaly
  const auto A2 = input.get<int>("A2");
  if (A2) {
    std::cout << "\n-------------------------------\n";
    std::cout << "Differential hyperfine anomaly:\n";
    fmt::print("A1 = {}\n", wf.Anuc());
    fmt::print("A2 = {}\n", *A2);

    std::cout << "\nRe-solving HF for second isotope:\n";
    // Set nuclear type for second nucleus. If given
    const auto nuc_type = input.get({"Nucleus2"}, "type", std::string{"Fermi"});
    // Get default nucleus:
    auto nucleus = Nuclear::Nucleus{wf.Znuc(), *A2, nuc_type};
    // over-ride default options
    const auto rrms = input.get<double>({"Nucleus2"}, "rrms");
    const auto t = input.get<double>({"Nucleus2"}, "t");
    const auto c_hdr = input.get<double>({"Nucleus2"}, "c");
    if (t) {
      nucleus.t() = *t;
    }
    if (rrms) {
      nucleus.set_rrms(*rrms);
    }
    if (c_hdr) {
      // this will over-ride given rms
      nucleus.set_rrms(Nuclear::rrms_formula_c_t(*c_hdr, nucleus.t()));
    }
    // If A or given rrms are zero, explicitely set to pointlike nucleus
    // This isn't required, but makes output more explicit
    if (nucleus.a() == 0.0 || nucleus.r_rms() == 0.0) {
      nucleus.t() = 0.0;
      nucleus.set_rrms(0.0);
      nucleus.type() = Nuclear::ChargeDistro::point;
    }

    auto wf2 = Wavefunction(wf.grid().params(), std::move(nucleus),
                            wf.alpha() / PhysConst::alpha);

    std::cout << "A1 = " << wf.nucleus() << "\n";
    std::cout << "A2 = " << wf2.nucleus() << "\n";

    wf2.solve_core("HartreeFock", wf.vHF()->x_Breit(),
                   wf.coreConfiguration_nice());
    wf2.solve_valence(DiracSpinor::state_config(wf.valence()));
    if (wf.Sigma()) {
      wf2.copySigma(wf.Sigma());
      wf2.hartreeFockBrueckner();
    }
    wf2.printValence();
    wf2.basis() = wf.basis(); // OK??

    std::cout << "\nNote: Assuming correlation potential and basis (for RPA) "
                 "is the "
                 "same accross isotopes! Test by swapping A1 and A2\n";

    auto hfs2_options = input.getBlock("hfs2_options") ?
                          *input.getBlock("hfs2_options") :
                          IO::InputBlock{};

    // Generate hyperfine operator (including BW effect from input)
    std::cout << "\nIsotope 2:";
    const auto h2 = DiracOperator::generate("hfs", hfs2_options, wf2);

    std::unique_ptr<ExternalField::DiagramRPA> rpa2{nullptr};
    if (rpaQ) {
      std::cout << "\nIncluding RPA (diagram method) - must have basis\n";

      rpa2 = std::make_unique<ExternalField::DiagramRPA>(h2.get(), rpa.get());

      std::cout << "Solving RPA:\n";
      std::cout << "For finite-magnetisation (isotope 2):\n";
      rpa2->solve_core(0.0, 100, true);
    }

    // std::cout << "\n-------------------------------\n";
    // fmt::print("Bohr-Weisskopf effect, Z={}, A={}\n", wf2.Znuc(), wf2.Anuc());
    // fmt::print("R_charge = {:.4f} fm [Rrms_charge = {:.4f} fm]\n",
    //            std::sqrt(5.0 / 3) * wf2.get_rrms(), wf2.get_rrms());

    // // Build pointlike HFS operator
    // IO::InputBlock h2_0_options{"h0", "F(r)=pointlike;"};
    // // Ensure same parameters as other operator
    // {
    //   if (hfs2_options.has_option("mu"))
    //     h2_0_options.add(*hfs_options->getOption("mu"));
    //   if (hfs2_options.has_option("I"))
    //     h2_0_options.add(*hfs_options->getOption("I"));
    //   if (hfs2_options.has_option("Q"))
    //     h2_0_options.add(*hfs_options->getOption("Q"));
    //   if (hfs2_options.has_option("k"))
    //     h2_0_options.add(*hfs_options->getOption("k"));
    //   if (hfs2_options.has_option("print"))
    //     h2_0_options.add(*hfs_options->getOption("print"));
    // }
    // const auto h2_0 = DiracOperator::generate("hfs", h2_0_options, wf);

    // // 1. Calculate BW effect
    // BW_effect(wf2.valence(), h2_0.get(), nullptr, h2.get(), nullptr);

    std::cout << "\nDifferential hyperfine anomaly:\n";
    fmt::print("\n{:4s} {:11s} {:11s} {:8s}  {:11s} {:11s} {:8s}\n", "",
               "A1_HF/g", "A2_HF/g", "1D2_HF(%)", "A1_RPA/g", "A2_RPA/g",
               "1D2_RPA(%)");
    for (const auto &Fv : wf.valence()) {
      const auto rme_to_A = DiracOperator::Hyperfine::convert_RME_to_AB(
        h0->rank(), Fv.kappa(), Fv.kappa());

      const auto Ahf1 = (h->reducedME(Fv, Fv) + (rpa ? rpa->dV(Fv, Fv) : 0.0)) *
                        rme_to_A / h->getc();
      const auto A0hf1 = h->reducedME(Fv, Fv) * rme_to_A / h->getc();

      const auto Fv2 = *wf2.getState(Fv.n(), Fv.kappa());

      const auto Ahf2 =
        (h2->reducedME(Fv2, Fv2) + (rpa2 ? rpa2->dV(Fv2, Fv2) : 0.0)) *
        rme_to_A / h2->getc();
      const auto A0hf2 = h2->reducedME(Fv2, Fv2) * rme_to_A / h2->getc();

      const auto D12 = 100.0 * (Ahf1 / Ahf2 - 1.0);
      const auto D012 = 100.0 * (A0hf1 / A0hf2 - 1.0);

      fmt::print("{:4s} {:.5e} {:.5e}  {:+.5f}", Fv.shortSymbol(), A0hf1, A0hf2,
                 D012);
      if (rpaQ)
        fmt::print("  {:.5e} {:.5e}  {:+.5f}", Ahf1, Ahf2, D12);
      std::cout << "\n";
    }
    // no longer use rpa2; release memory
    rpa2.reset(nullptr);

    // Tune _both_ magnetic radii to fit target differential anomaly, 1D2
    if (input.get("1D2_target")) {
      std::cout << "\n-------------------------------\n";
      std::cout << "Tune magnetic radii of both isotopes to match 1D2:\n";

      const auto [a, str_1D2] =
        *input.get<std::array<std::string, 2>>("1D2_target");
      // fails if eps_s invalid 'double':
      const auto target_1D2 = str_1D2 == "" ? 0.0 : std::stod(str_1D2);
      const auto pFv = wf.getState(a);
      std::cout << "1D2 target: " << target_1D2 << "\n";
      if (!pFv) {
        std::cout << "Invalid state: " << a << "\n";
        return;
      }
      const auto &Fv1 = *pFv;
      const auto &Fv2 = *wf2.getState(a);

      const auto min_max_steps =
        input.get<std::vector<double>>("min_max_steps", {0.9, 1.5, 10});
      const int num_steps = int(min_max_steps.at(2));
      const double r1_i = min_max_steps.at(0) * wf.get_rrms();
      const double r1_f = min_max_steps.at(1) * wf.get_rrms();
      const double dr = (r1_f - r1_i) / (num_steps - 1);

      std::unique_ptr<ExternalField::DiagramRPA> rpa1{nullptr}, rpa22{nullptr};

      std::cout << "\nRm_rms(1) Rm_rms(2) 1D2\n";
      for (int i = 0; i < num_steps; ++i) {
        auto r1 = r1_i + i * dr;

        // 1. Calculate A1(r1)
        // 2. Fit r2 : 1D2 = target

        auto h1_options = hfs_options ? *hfs_options : IO::InputBlock{};
        h1_options.add(IO::Option{"rrms", std::to_string(r1)});
        h1_options.add(IO::Option{"print", "false"});
        const auto h1 = DiracOperator::generate("hfs", h1_options, wf);

        if (rpaQ) {
          if (rpa1) {
            rpa1->update_t0s(h1.get());
          } else {
            rpa1 =
              std::make_unique<ExternalField::DiagramRPA>(h1.get(), rpa.get());
            rpa22 =
              std::make_unique<ExternalField::DiagramRPA>(h1.get(), rpa1.get());
          }
          rpa1->solve_core(0.0, 40, false);
        }

        const auto rme_to_A = DiracOperator::Hyperfine::convert_RME_to_AB(
          h0->rank(), Fv1.kappa(), Fv1.kappa());

        const auto A1 =
          (h1->reducedME(Fv1, Fv1) + (rpa1 ? rpa1->dV(Fv1, Fv1) : 0.0)) *
          rme_to_A / h1->getc();

        // double r2 = r1;
        double hf_1D2 = -1.0;
        const auto delta_1D2 = [&](double r2) {
          auto h2_options = hfs2_options;
          h2_options.add(IO::Option{"rrms", std::to_string(r2)});
          h2_options.add(IO::Option{"print", "false"});
          auto h22 = DiracOperator::generate("hfs", h2_options, wf2);

          if (rpaQ) {
            rpa22->update_t0s(h22.get());
            rpa22->solve_core(0.0, 60, false);
          }

          const auto A22 =
            (h22->reducedME(Fv2, Fv2) + (rpa22 ? rpa22->dV(Fv2, Fv2) : 0.0)) *
            rme_to_A / h22->getc();
          hf_1D2 = 100.0 * (A1 / A22 - 1.0);
          return hf_1D2 - target_1D2;
        };

        const auto r2_0 = /*(r1 / wf.get_rrms())**/ wf2.get_rrms();
        const auto [r2_fitted, error] = qip::Newtons(
          delta_1D2, r2_0, {0.5 * r2_0, 2.0 * r2_0}, 1.0e-5, 0.001, 150);

        const auto D12_error = std::abs(hf_1D2 / target_1D2 - 1.0);
        const auto flag = D12_error > 0.01 ? "***" : "";
        fmt::print("{:.5f}   {:.5f}   {:.4f}  {}\n", r1, r2_fitted, hf_1D2,
                   flag);
      }
    }
  }
  //
}

//==============================================================================
void b_plot(const IO::InputBlock &input, const Wavefunction &wf) {
  //

  input.check({
    {"r1", "minimum r_rms in fm [0.5*rrms]"},
    {"r2", "maximum r_rms in fm [1.5*rrms]"},
    {"t", "skin thickness [2.3]"},
    {"num_steps", "Number of steps along (r1,r2) [100]"},
  });
  if (input.has_option("help")) {
    return;
  }

  const auto r0 = wf.get_rrms();

  const auto r1 = input.get("r1", 0.5 * r0);
  const auto r2 = input.get("r2", 1.5 * r0);
  const auto t = input.get("t", wf.nucleus().t());

  const auto num_steps = input.get("num_steps", 100);

  const auto core_str = wf.coreConfiguration_nice();
  const auto valence_str = DiracSpinor::state_config(wf.valence());
  const auto grid = wf.grid_sptr();
  auto nucleus = wf.nucleus();

  // For grid uncertainty
  const std::size_t small_points =
    std::max(2000ul, std::size_t(0.2 * (double)grid->num_points()));
  const double small_r0 = std::min(1.0e-6, 10.0 * grid->r0());

  const auto small_grid1 = std::make_shared<const Grid>(Grid(
    grid->r0(), grid->rmax(), small_points, grid->type(), grid->loglin_b()));

  const auto small_grid2 = std::make_shared<const Grid>(
    Grid(small_r0, grid->rmax(), small_points, grid->type(), grid->loglin_b()));

  using Data = std::vector<std::array<double, 3>>;

  std::vector<double> rs;
  Data KSs, KLs, dKSs, dKLs;

  const auto dr = (r2 - r1) / (num_steps - 1);
  for (int ii = 0; ii < num_steps; ++ii) {

    double r_rms = r1 + ii * dr;

    nucleus.set_rrms(r_rms);
    nucleus.t() = t;

    const auto R0 = std::sqrt(5.0 / 3) * r_rms;

    Wavefunction wf_t(grid, nucleus);
    Wavefunction wf_s(small_grid1, nucleus);
    Wavefunction wf_s2(small_grid2, nucleus);

    std::cout << "\n--------------------------\n\n";
    std::cout << wf_t.nucleus() << "\n";

    wf_t.solve_core("HartreeFock", 0.0, core_str, 1.0e-13, true);
    wf_t.solve_valence(valence_str);
    wf_t.printValence();
    const auto &v = wf_t.valence().front();

    wf_s.solve_core("HartreeFock", 0.0, core_str, 1.0e-13, true);
    wf_s.solve_valence(valence_str);
    wf_s.printValence();
    const auto &v_s = wf_s.valence().front();

    wf_s2.solve_core("HartreeFock", 0.0, core_str, 1.0e-13, true);
    wf_s2.solve_valence(valence_str);
    wf_s2.printValence();
    const auto &v_s2 = wf_s2.valence().front();

    std::cout << "\nMain:\n";
    const auto [KS, KL] = b_moments("", v, R0, 20);

    rs.push_back(r_rms);
    KSs.push_back(KS);
    KLs.push_back(KL);

    // Bunch: for uncertainty
    std::cout << "\nFewer terms:\n";
    const auto [KS1, KL1] = b_moments("", v, R0, 14);

    std::cout << "\nSmaller grid (same r0):\n";         // This dominates!
    const auto [KS2, KL2] = b_moments("", v_s, R0, 20); // small grid

    std::cout << "\nSmaller grid (larger r0):\n";
    const auto [KS5, KL5] = b_moments("", v_s2, R0, 20); // small grid

    // different maximum fit radii
    std::cout << "\nFit to 0.9R:\n";
    const auto [KS3, KL3] = b_moments("", v, 0.9 * R0, 20);
    std::cout << "\nFit to 1.1R:\n";
    const auto [KS4, KL4] = b_moments("", v, 1.1 * R0, 20);

    // uncerainty:
    auto KSx = KS; // just for array shape
    auto KLx = KL;
    for (std::size_t i = 0; i < KSx.size(); ++i) {
      KSx[i] = std::abs(0.5 * qip::max_difference(KS[i], KS1[i], KS2[i], KS3[i],
                                                  KS4[i], KS5[i]));
      KLx[i] = std::abs(0.5 * qip::max_difference(KL[i], KL1[i], KL2[i], KL3[i],
                                                  KL4[i], KL5[i]));
    }
    dKSs.push_back(KSx);
    dKLs.push_back(KLx);
  }

  std::cout << "\n\n";
  for (std::size_t i = 0; i < std::size_t(num_steps); ++i) {
    fmt::print("{:.6f} {:+.6e} {:+.6e} {:+.6e} {:+.2e} {:+.2e} {:+.2e} {:+.6e} "
               "{:+.6e} {:+.6e} {:+.2e} {:+.2e} {:+.2e}\n",
               rs[i],                              //
               KSs[i][0], KSs[i][1], KSs[i][2],    //
               dKSs[i][0], dKSs[i][1], dKSs[i][2], //
               KLs[i][0], KLs[i][1], KLs[i][2],    //
               dKLs[i][0], dKLs[i][1], dKLs[i][2]);
  }
  std::cout << "\n";
}

} // namespace Module
