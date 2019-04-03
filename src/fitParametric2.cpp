#include "ATI_atomInfo.h"
#include "ChronoTimer.h"
#include "ElectronOrbitals.h"
#include "FPC_physicalConstants.h"
#include "HartreeFockClass.h"
#include "NumCalc_quadIntegrate.h"
#include "PRM_parametricPotentials.h"
#include <cmath>
#include <fstream>
#include <iostream>
/*
Finds the best-fit parameter values for the Green or Tietz potentials
*/

int main(int argc, char *argv[]) {
  ChronoTimer sw; // start the overall timer

  std::string input_file = (argc > 1) ? argv[1] : "fitParametric.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // bool sphere = true; // Finite nucleus? makes no difference??

  // Input parameters:
  std::string Z_str, str_core;
  int num_val, l_max; // valence states to calc
  int A;
  int ngp;
  double r0, rmax;
  int igreen;
  std::vector<int> in_n, in_k;
  std::vector<double> in_en;
  std::ifstream ifile;
  ifile.open("fitParametric2.in"); // input file
  {
    std::string junk;
    ifile >> Z_str >> A;
    getline(ifile, junk);
    ifile >> str_core;
    getline(ifile, junk);
    ifile >> num_val >> l_max;
    getline(ifile, junk);
    ifile >> r0 >> rmax >> ngp;
    getline(ifile, junk);
    ifile >> igreen;
    getline(ifile, junk);
    // int nstates;
    // ifile >> nstates;
    // getline(ifile, junk);
    // for (int ns = 0; ns < nstates; ns++) {
    //   int tn, tk;
    //   double te;
    //   ifile >> tn >> tk >> te;
    //   getline(ifile, junk);
    //   in_n.push_back(tn);
    //   in_k.push_back(tk);
    //   in_en.push_back(te);
    // }
  }
  ifile.close();

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

  // convergence parameters for finding best-fit H and d (or t and g)
  double eps = 1.e-5;
  int max_its = 100;

  printf("\n Finding best-fit parameters for %s potential, %s Z=%i\n",
         which.c_str(), Z_str.c_str(), Z);
  printf("*********************************************************\n");

  if (do_HF) {
    ElectronOrbitals hfwf(Z, A, ngp, r0, rmax);
    HartreeFock hf(hfwf, str_core, 1.e-6);
    // for (auto &phi : hfwf.orbitals) {
    //   if (phi.k < 0 && phi.k != -1)
    //     continue;
    //   in_n.push_back(phi.n);
    //   in_k.push_back(phi.k);
    //   in_en.push_back(phi.en);
    // }

    // Create list of valence states to solve for
    // Solves for lowest num_val states with given l
    std::vector<std::vector<int>> lst;
    if ((int)hfwf.Ncore() >= hfwf.Znuc())
      num_val = 0;
    for (int l = 0; l <= l_max; l++) {
      int n0 = hfwf.maxCore_n(l) + 1;
      if (n0 == 1)
        n0 += l;
      for (int nv = 0; nv < num_val; nv++) {
        for (int tk = 0; tk < 2; tk++) {     // loop over k
          int ka = (tk == 0) ? l : -(l + 1); // kappa
          if (ka == 0)
            continue; // no j=l-1/2 for l=s
          lst.push_back({n0 + nv, ka});
        }
      }
    }
    for (const auto &nk : lst) {
      int n = nk[0];
      int k = nk[1];
      hf.solveValence(n, k);
    }

    for (auto &phi : hfwf.orbitals) {
      if (phi.k < 0 && phi.k != -1)
        continue;
      in_n.push_back(phi.n);
      in_k.push_back(phi.k);
      in_en.push_back(phi.en);
    }
  }

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
        // if (sphere)
        // wf.formNuclearPotential(NucleusType::spherical);
        if (green)
          for (auto r : wf.rgrid.r)
            wf.vdir.push_back(PRM::green(Z, r, H, d));
        else
          for (auto r : wf.rgrid.r)
            wf.vdir.push_back(PRM::tietz(Z, r, H, d));
        double fx = 0;
        for (std::size_t ns = 0; ns < in_n.size(); ns++) {
          wf.solveLocalDirac(in_n[ns], in_k[ns], in_en[ns]);
          fx += pow((wf.orbitals[ns].en - in_en[ns]), 2);
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

  double H = best_H;
  double d = best_d;
  std::cout << "\nBest fit parameters for ";
  if (green)
    printf("Green: \n  H=%7.5f  d=%7.5f\n\n", H, d);
  else
    printf("Tietz: \n  t=%7.5f  g=%7.5f\n\n", H, d);

  // Now, solve using the above-found best-fit parameters:
  ElectronOrbitals wf(Z, A, ngp, r0, rmax);
  //  if (sphere)
  //    wf.formNuclearPotential(NucleusType::spherical);
  if (green)
    for (auto r : wf.rgrid.r)
      wf.vdir.push_back(PRM::green(Z, r, H, d));
  else
    for (auto r : wf.rgrid.r)
      wf.vdir.push_back(PRM::tietz(Z, r, H, d));
  for (std::size_t ns = 0; ns < in_n.size(); ns++)
    wf.solveLocalDirac(in_n[ns], in_k[ns], in_en[ns]);

  printf(" n l_j    k Rinf its  eps     En (au)            En (/cm)\n");
  double en0 = wf.orbitals.front().en;
  int i = 0;
  for (auto &phi : wf.orbitals) {
    auto njl = phi.symbol().c_str();
    double rinf = wf.rinf(phi);
    double eni = phi.en;
    double enT = in_en[i++];
    printf("%7s %2i  %3.0f %3i  %5.0e  %10.4f  %13.3f  %8.2f%%\n", njl, phi.k,
           rinf, phi.its, wf.orbitals[i].eps, eni,
           (eni - en0) * FPC::Hartree_invcm, 100. * (enT - eni) / enT);
  }

  printf("{%i, %.3f, %.3f}\n", Z, H, d);

  std::cout << "\nTime: " << sw.reading_str() << "\n";

  return 0;
}
