#include "Modules/fitParametric.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <tuple>
#include <vector>

namespace Module {
//******************************************************************************
void fitParametric(const IO::InputBlock &input, const Wavefunction &wf) {

  input.checkBlock_old({"statesToFit", "method", "fitWorst"});

  auto which_states = input.get<std::string>("statesToFit", "core");
  auto which_method = input.get<std::string>("method", "Green");

  auto fit_worst = input.get("fitWorst", true);

  std::vector<AtomData::DiracSEnken> states;
  // bool fit_worst = true; // XXX update to be input?

  if (which_states == "core" || which_states == "both") {
    for (const auto &phi : wf.core)
      states.emplace_back(phi.n, phi.k, phi.en());
  }
  if (which_states == "valence" || which_states == "both") {
    for (const auto &phi : wf.valence)
      states.emplace_back(phi.n, phi.k, phi.en());
  }
  // XXX Can easily update later to take a user-given list!

  bool green = true;
  if (which_method != "Green") {
    green = false;
    which_method = "Tietz";
  }

  if (!states.empty()) {
    double H, d;
    auto Z = wf.Znuc();
    // auto A = wf.Anuc();
    auto num_points = wf.rgrid->num_points();
    auto r0 = wf.rgrid->r0();
    auto rmax = wf.rgrid->rmax();
    GridParameters gp(num_points, r0, rmax, wf.rgrid->loglin_b(),
                      wf.rgrid->type());
    auto nuc_params = wf.get_nuclearParameters();

    std::tie(H, d) =
        FitParametric::performFit(states, Z, gp, nuc_params, green, fit_worst);

    std::cout << "\nBest fit parameters (core) for ";
    if (green)
      printf("Green: \n  H=%7.5f  d=%7.5f\n", H, d);
    else
      printf("Tietz: \n  t=%7.5f  g=%7.5f\n", H, d);

    // Now, solve using the above-found best-fit parameters:
    Wavefunction wf_prm(gp, nuc_params);
    if (green)
      for (auto r : wf_prm.rgrid->r())
        wf_prm.vdir.push_back(Parametric::green(Z, r, H, d));
    else
      for (auto r : wf_prm.rgrid->r())
        wf_prm.vdir.push_back(Parametric::tietz(Z, r, H, d));
    for (auto &nk : states) {
      wf_prm.solveNewValence(nk.n, nk.k, nk.en);
    }

    std::cout << "\nUsing " << which_method
              << " Parametric potential: " << wf.atom() << "\n";

    double en0 = wf_prm.valence.front().en();
    int i = 0;
    for (auto &phi : wf_prm.valence) {
      // const auto njl = phi.symbol().c_str();
      // double rinf = wf.rinf(phi);
      double rinf = wf.rgrid->r()[phi.max_pt() - 1];
      double eni = phi.en();
      double enT = states[std::size_t(i++)].en;
      printf("%7s %2i  %3.0f %3i  %5.0e  %13.7f  %11.4f %8.2f%%\n",
             phi.symbol().c_str(), phi.k, rinf, phi.its(), phi.eps(), eni,
             (eni - en0) * PhysConst::Hartree_invcm, 100. * (enT - eni) / enT);
    }

    printf("\nParams: {%i, %.3f, %.3f}\n", Z, H, d);
  }

  return;
}

//******************************************************************************
std::tuple<double, double>
FitParametric::performFit(const std::vector<AtomData::DiracSEnken> &states,
                          int Z, const GridParameters &gp,
                          const Nuclear::Parameters &nuc_params, bool green,
                          bool fit_worst) {

  std::cout << "\nPerforming fit (for ";
  if (green)
    std::cout << "Green ";
  else
    std::cout << "Teitz ";
  std::cout << "potential).\n";
  if (fit_worst)
    std::cout << "Fitting for worst state.\n";
  else
    std::cout << "Fitting by sum of (relative) squares.\n";

  // convergence parameters for finding best-fit H and d (or t and g)
  double eps = 1.e-6;
  int max_its = 100;

  double GHmin = 0.05, GHmax = 15.;
  double Gdmin = 0.01, Gdmax = 5.;

  double Hmin = GHmin, Hmax = GHmax;
  double dmin = Gdmin, dmax = Gdmax;

  double best_H = 0, best_d = 0;

  for (int nit = 0; nit < max_its; nit++) {
    const int n_array = 32;
    const int n_params = 3;
    double array[n_array][n_array][n_params];
    double dH = (Hmax - Hmin) / n_array;
    double dd = (dmax - dmin) / n_array;

// Solve equation for each (H,d), store values
#pragma omp parallel for
    for (int n = 0; n < n_array; n++) {
      double H = Hmin + n * dH;
      for (int m = 0; m < n_array; m++) {
        double d = dmin + m * dd;
        Wavefunction wf(gp, nuc_params);
        if (green)
          for (auto r : wf.rgrid->r())
            wf.vdir.push_back(Parametric::green(Z, r, H, d));
        else
          for (auto r : wf.rgrid->r())
            wf.vdir.push_back(Parametric::tietz(Z, r, H, d));
        // fits for the worst state
        double fx = 0;
        if (fit_worst) {
          for (std::size_t ns = 0; ns < states.size(); ns++) {
            wf.solveNewValence(states[ns].n, states[ns].k, states[ns].en);
            auto fx2 = fabs((wf.valence[ns].en() - states[ns].en) /
                            (wf.valence[ns].en() + states[ns].en));
            if (fx2 > fx)
              fx = fx2;
          }
        } else {
          // sum-of-squares
          for (std::size_t ns = 0; ns < states.size(); ns++) {
            wf.solveNewValence(states[ns].n, states[ns].k, states[ns].en);
            // fx += std::pow(wf.orbitals[ns].en() - states[ns].en, 2);
            fx += std::pow((wf.valence[ns].en() - states[ns].en) /
                               (wf.valence[ns].en() + states[ns].en),
                           2);
          }
        }
        array[n][m][0] = fx;
        array[n][m][1] = H;
        array[n][m][2] = d;
      }
    }

    // Find the "so-far" best-fit
    double bH = 0, bd = 0, bfx = 99999999.;
    for (int n = 0; n < n_array; n++) {
      for (int m = 0; m < n_array; m++) {
        if (array[n][m][0] < bfx) {
          bfx = array[n][m][0];
          bH = array[n][m][1];
          bd = array[n][m][2];
        }
      }
    }

    // Adjust H and d search range (hone-in):
    Hmin = bH - (0.2 * n_array) * dH;
    Hmax = bH + (0.2 * n_array) * dH;
    dmin = bd - (0.2 * n_array) * dd;
    dmax = bd + (0.2 * n_array) * dd;
    if (Hmax > GHmax)
      Hmax = GHmax;
    if (dmax > Gdmax)
      dmax = Gdmax;
    if (Hmin < GHmin)
      Hmin = GHmin;
    if (dmin < Gdmin)
      dmin = Gdmin;
    best_H = bH;
    best_d = bd;

    // Print progress to screen; quit if convergence OK
    printf("%2i %6.4f %6.4f  %.1e\n", nit, bH, bd, fmax(dH, dd));
    if (dH < eps && dd < eps)
      break;
  }

  return std::make_tuple(best_H, best_d);
}

} // namespace Module
