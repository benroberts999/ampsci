#include "ATI_atomInfo.hpp"
#include "ChronoTimer.hpp"
#include "ElectronOrbitals.hpp"
#include "FPC_physicalConstants.hpp"
#include "FileIO_fileReadWrite.hpp"
#include "PRM_localPotentials.hpp"
#include <iostream>
#include <tuple>

int main(int argc, char *argv[]) {
  ChronoTimer sw; // start the overall timer

  std::string input_file = (argc > 1) ? argv[1] : "parametricPotential.in";
  std::cout << "Reading input from: " << input_file << "\n";

  double varalpha = 1; // need same number as used for the fitting!

  // Input options
  std::string Z_str;
  int A;
  double r0, rmax;
  int ngp;
  double Tf, Tt, Tg; // Teitz potential parameters
  double Gf, Gh, Gd; // Green potential parameters
  int n_max, l_max;
  std::string str_core;

  { // Open and read the input file:
    auto tp = std::forward_as_tuple(Z_str, A, str_core, r0, rmax, ngp, Gf, Gh,
                                    Gd, Tf, Tt, Tg, n_max, l_max);
    FileIO::setInputParameters("parametricPotential.in", tp);
  }

  int Z = ATI::get_z(Z_str);
  if (A == 0)
    A = ATI::defaultA(Z); // if none given, get default A

  // Normalise the Teitz/Green weights:
  if (Gf != 0 || Tf != 0) {
    double TG_norm = Gf + Tf;
    Gf /= TG_norm;
    Tf /= TG_norm;
  }

  // If H,d etc are zero, use default values
  if (Gf != 0 && Gh == 0)
    PRM::defaultGreen(Z, Gh, Gd);
  if (Tf != 0 && Tt == 0)
    PRM::defaultTietz(Z, Tt, Tg);

  std::cout << "\nRunning parametric potential for " << Z_str << ", Z=" << Z
            << ", A=" << A << "\n";
  std::cout << "*************************************************\n";
  if (Gf != 0)
    printf("%3.0f%% Green potential: H=%.4f  d=%.4f\n", Gf * 100., Gh, Gd);
  if (Tf != 0)
    printf("%3.0f%% Tietz potential: T=%.4f  g=%.4f\n", Tf * 100., Tt, Tg);

  // Generate the orbitals object:
  ElectronOrbitals wf(Z, A, ngp, r0, rmax, varalpha);

  wf.rgrid.print();

  // Fill the electron part of the potential
  wf.vdir.clear();
  wf.vdir.reserve(wf.rgrid.ngp);
  // for (int i = 0; i < wf.rgrid.ngp; i++) {
  for (auto r : wf.rgrid.r) {
    double tmp = 0;
    if (Gf != 0)
      tmp += Gf * PRM::green(Z, r, Gh, Gd);
    if (Tf != 0)
      tmp += Tf * PRM::tietz(Z, r, Tt, Tg);
    wf.vdir.push_back(tmp);
  }

  // Solve for core states
  wf.solveInitialCore(str_core);

  // Calculate the valence (and excited) states
  for (int n = 1; n <= n_max; n++) {
    for (int l = 0; l <= l_max; l++) {
      if (l + 1 > n)
        continue;
      for (int tk = 0; tk < 2; tk++) {
        int k;
        if (tk == 0)
          k = l;
        else
          k = -(l + 1);
        if (k == 0)
          continue;
        if (wf.isInCore(n, k))
          continue;
        wf.solveLocalDirac(n, k);
      }
    }
  }

  // make list of energy indices in sorted order:
  // std::vector<int>
  auto sorted_list = wf.sortedEnergyList();

  // Output results:
  std::cout << "\n n l_j    k Rinf its    eps      En (au)        En (/cm)\n";
  int m = 0;
  for (auto i : sorted_list) {
    auto &psi = wf.orbitals[i];
    if (m++ == (int)wf.coreIndexList.size()) {
      std::cout << " ######### Valence: ######\n";
      std::cout
          << " n l_j    k Rinf its    eps      En (au)        En (/ cm)\n";
    }
    double rinf = wf.rinf(psi);
    printf("%7s %2i  %3.0f %3i  %5.0e  %11.5f %15.3f\n", psi.symbol().c_str(),
           psi.k, rinf, psi.its, wf.orbitals[i].eps, psi.en,
           psi.en * FPC::Hartree_invcm);
  }

  std::cout << "\n " << sw.reading_str() << "\n";
  return 0;
}
