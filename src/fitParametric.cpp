#include "ATI_atomInfo.hpp"
#include "ChronoTimer.hpp"
#include "ElectronOrbitals.hpp"
#include "FPC_physicalConstants.hpp"
#include "FileIO_fileReadWrite.hpp"

#include "HartreeFockClass.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "Parametric_potentials.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
/*
Finds the best-fit parameter values for the Green or Tietz potentials
*/

// struct EOnken {
//   int n;
//   int k;
//   double en;
//   EOnken(int in_n, int in_k, double in_en) : n(in_n), k(in_k), en(in_en){};
//   EOnken(){};
// };

std::tuple<double, double> performFit(const std::vector<EOnken> &states, int Z,
                                      int A, int ngp, double r0, double rmax,
                                      bool green, bool fit_worst);

//******************************************************************************
int main(int argc, char *argv[]) {
  ChronoTimer sw; // start the overall timer

  std::string input_file = (argc > 1) ? argv[1] : "fitParametric.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input parameters:
  std::string Z_str, str_core;
  int A;
  int ngp;
  double r0, rmax;
  int igreen;
  std::vector<EOnken> states;
  std::vector<EOnken> val_states;
  bool fit_worst;

  auto in_str_list = FileIO::readInputFile_byEntry(input_file);
  {
    int fit_type;
    auto tp = std::forward_as_tuple(Z_str, A, str_core, r0, rmax, ngp, igreen,
                                    fit_type);
    FileIO::stringstreamVectorIntoTuple(in_str_list, tp);
    fit_worst = (fit_type == 0) ? true : false;

    auto n_els = std::tuple_size<decltype(tp)>::value; // fuck you c++
    auto n_val_els = in_str_list.size() - n_els;
    if (n_val_els % 3 != 0)
      std::cout << "Error with valence input. Wrong number of input elements\n";
    for (std::size_t i = n_els; i < in_str_list.size(); i += 3) {
      auto tv = std::vector<std::string>(in_str_list.begin() + i,
                                         in_str_list.begin() + i + 3);
      EOnken tmp;
      auto tp2 = std::forward_as_tuple(tmp.n, tmp.k, tmp.en);
      FileIO::stringstreamVectorIntoTuple(tv, tp2);
      val_states.push_back(tmp);
    }
  }

  bool do_HF = true;
  if (str_core == "na")
    do_HF = false;

  int Z = ATI::get_z(Z_str);
  if (Z == 0)
    return 2;

  bool green = true;
  std::string which = "Green";
  if (igreen == 1) {
    green = false;
    which = "Tietz";
  }

  printf("\n Finding best-fit parameters for %s potential, %s Z=%i\n",
         which.c_str(), Z_str.c_str(), Z);
  std::cout << "*********************************************************\n";

  if (do_HF) {
    ElectronOrbitals hfwf(Z, A, ngp, r0, rmax);
    HartreeFock hf(hfwf, str_core, 1.e-9);
    for (auto &phi : hfwf.core_orbitals) {
      // // don't fit for both j=l+/-1/2, just one!
      // if (phi.k < 0 && phi.k != -1)
      //   continue;
      states.emplace_back(phi.n, phi.k, phi.en);
    }
  }

  if (states.size() > 0) {
    double H, d;
    std::tie(H, d) = performFit(states, Z, A, ngp, r0, rmax, green, fit_worst);

    std::cout << "\nBest fit parameters (core) for ";
    if (green)
      printf("Green: \n  H=%7.5f  d=%7.5f\n\n", H, d);
    else
      printf("Tietz: \n  t=%7.5f  g=%7.5f\n\n", H, d);

    // Now, solve using the above-found best-fit parameters:
    ElectronOrbitals wf(Z, A, ngp, r0, rmax);
    if (green)
      for (auto r : wf.rgrid.r)
        wf.vdir.push_back(Parametric::green(Z, r, H, d));
    else
      for (auto r : wf.rgrid.r)
        wf.vdir.push_back(Parametric::tietz(Z, r, H, d));
    for (auto &nk : states) {
      wf.solveNewValence(nk.n, nk.k, nk.en);
    }

    printf(" nl_j    k  Rinf its   eps     En (au)    \n");
    // double en0 = wf.orbitals.front().en;
    int i = 0;
    for (auto &phi : wf.core_orbitals) {
      auto njl = phi.symbol().c_str();
      double rinf = wf.rinf(phi);
      double eni = phi.en;
      double enT = states[i++].en;
      printf("%7s %2i  %3.0f %3i  %5.0e  %10.4f  %8.2f%%\n", njl, phi.k, rinf,
             phi.its, phi.eps, eni, 100. * (enT - eni) / enT);
    }

    printf("\nCore: {%i, %.3f, %.3f}\n", Z, H, d);
  }

  if (val_states.size() > 0) {
    std::cout << "\n";
    double H, d;
    std::tie(H, d) =
        performFit(val_states, Z, A, ngp, r0, rmax, green, fit_worst);

    std::cout << "\nBest fit parameters (valence) for ";
    if (green)
      printf("Green: \n  H=%7.5f  d=%7.5f\n\n", H, d);
    else
      printf("Tietz: \n  t=%7.5f  g=%7.5f\n\n", H, d);

    // Now, solve using the above-found best-fit parameters:
    ElectronOrbitals wf(Z, A, ngp, r0, rmax);
    if (green)
      for (auto r : wf.rgrid.r)
        wf.vdir.push_back(Parametric::green(Z, r, H, d));
    else
      for (auto r : wf.rgrid.r)
        wf.vdir.push_back(Parametric::tietz(Z, r, H, d));
    for (auto &nk : val_states) {
      wf.solveNewValence(nk.n, nk.k, nk.en);
    }

    printf(" nl_j    k  Rinf its   eps     En (au)        En (/cm)\n");
    double en0 = wf.valence_orbitals.front().en;
    int i = 0;
    for (auto &phi : wf.valence_orbitals) {
      auto njl = phi.symbol().c_str();
      double rinf = wf.rinf(phi);
      double eni = phi.en;
      double enT = val_states[i++].en;
      printf("%7s %2i  %3.0f %3i  %5.0e  %13.7f  %11.4f %8.2f%%\n", njl, phi.k,
             rinf, phi.its, phi.eps, eni, (eni - en0) * FPC::Hartree_invcm,
             100. * (enT - eni) / enT);
    }

    printf("\nValence: {%i, %.3f, %.3f}\n", Z, H, d);
  }

  std::cout << "\nTime: " << sw.reading_str() << "\n";

  return 0;
}

//******************************************************************************
std::tuple<double, double> performFit(const std::vector<EOnken> &states, int Z,
                                      int A, int ngp, double r0, double rmax,
                                      bool green, bool fit_worst) {

  std::cout << "Performing fit (for ";
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
        ElectronOrbitals wf(Z, A, ngp, r0, rmax);
        if (green)
          for (auto r : wf.rgrid.r)
            wf.vdir.push_back(Parametric::green(Z, r, H, d));
        else
          for (auto r : wf.rgrid.r)
            wf.vdir.push_back(Parametric::tietz(Z, r, H, d));
        // fits for the worst state
        double fx = 0;
        if (fit_worst) {
          for (std::size_t ns = 0; ns < states.size(); ns++) {
            wf.solveNewValence(states[ns].n, states[ns].k, states[ns].en);
            auto fx2 = fabs((wf.valence_orbitals[ns].en - states[ns].en) /
                            (wf.valence_orbitals[ns].en + states[ns].en));
            if (fx2 > fx)
              fx = fx2;
          }
        } else {
          // sum-of-squares
          for (std::size_t ns = 0; ns < states.size(); ns++) {
            wf.solveNewValence(states[ns].n, states[ns].k, states[ns].en);
            // fx += pow(wf.orbitals[ns].en - states[ns].en, 2);
            fx += pow((wf.valence_orbitals[ns].en - states[ns].en) /
                          (wf.valence_orbitals[ns].en + states[ns].en),
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
