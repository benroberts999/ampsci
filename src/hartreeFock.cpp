#include "ATI_atomInfo.h"
#include "ChronoTimer.h"
#include "ElectronOrbitals.h"
#include "FPC_physicalConstants.h"
#include "FileIO_fileReadWrite.h"
#include "HartreeFockClass.h"
#include "Nucleus.h"
#include "NumCalc_quadIntegrate.h"
#include "PRM_parametricPotentials.h"
#include <cmath>
#include <iostream>
#include <tuple>

int main(int argc, char *argv[]) {
  ChronoTimer timer; // start the overall timer

  std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  std::string Z_str;
  int A;
  std::string str_core;
  double r0, rmax;
  int ngp;
  double eps_HF;      // HF convergance
  int num_val, l_max; // valence states to calc
  double varalpha, varalpha2;
  bool exclude_exchange;

  { // Open and read the input file:
    int i_excl_ex;
    auto tp = std::forward_as_tuple(Z_str, A, str_core, r0, rmax, ngp, eps_HF,
                                    num_val, l_max, varalpha2, i_excl_ex);
    FileIO::setInputParameters(input_file, tp);
    exclude_exchange = i_excl_ex == 1 ? true : false;
  }

  // Change varAlph^2 to varalph
  if (varalpha2 == 0)
    varalpha2 = 1.e-10;
  varalpha = sqrt(varalpha2);

  int Z = ATI::get_z(Z_str);
  if (A == -1)
    A = ATI::defaultA(Z); // if none given, get default A

  if (exclude_exchange)
    std::cout << "\nRunning Hartree (excluding exchange) for ";
  else
    std::cout << "\nRunning Hartree-Fock for ";
  std::cout << Z_str << "; Z=" << Z << " A=" << A << "\n"
            << "*************************************************\n";

  // Generate the orbitals object:
  ElectronOrbitals wf(Z, A, ngp, r0, rmax, varalpha);
  wf.rgrid.print();
  std::cout << "\n";

  // Solve Hartree equations for the core:
  timer.start(); // start the timer for HF
  HartreeFock hf(wf, str_core, eps_HF, exclude_exchange);
  double core_energy = hf.calculateCoreEnergy();
  std::cout << "core: " << timer.lap_reading_str() << "\n";

  // Create list of valence states to solve for
  // Solves for lowest num_val states with given l
  std::vector<std::vector<int>> lst;
  if ((int)wf.Ncore() >= wf.Znuc())
    num_val = 0;
  for (int l = 0; l <= l_max; l++) {
    int n0 = wf.maxCore_n(l) + 1;
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

  // Solve for the valence states:
  timer.start();
  for (const auto &nk : lst) {
    int n = nk[0];
    int k = nk[1];
    hf.solveValence(n, k);
  }
  if (lst.size() > 0)
    std::cout << "Valence: " << timer.lap_reading_str() << "\n";

  // make list of energy indices in sorted order:
  std::vector<int> sorted_by_energy_list;
  wf.sortedEnergyList(sorted_by_energy_list, true);

  // Output results:
  std::cout << "\nCore: " << Z_str << ", Z=" << Z << " A=" << A << "\n";
  std::cout << "     state   k   Rinf its    eps       En (au)      En (/cm)\n";
  bool val = false;
  double en_lim = 0;
  int count = 0;
  for (int i : sorted_by_energy_list) {
    auto &phi = wf.orbitals[i];
    ++count;
    if (val && en_lim == 0)
      en_lim = fabs(phi.en); // give energies wrt core
    double rinf = wf.rinf(phi);
    // auto xxxx = phi.symbol().c_str();
    printf("%2i) %7s %2i  %5.1f %3i  %5.0e %13.7f %13.1f", i,
           phi.symbol().c_str(), phi.k, rinf, phi.its, phi.eps, phi.en,
           phi.en * FPC::Hartree_invcm);
    if (val)
      printf(" %10.2f\n", (phi.en + en_lim) * FPC::Hartree_invcm);
    else
      std::cout << "\n";
    if (count == (int)wf.coreIndexList.size()) {
      // printf("E_core = %.5f au\n", core_energy);
      std::cout << "E_core = " << core_energy << " au\n";
      if (wf.coreIndexList.size() == wf.orbitals.size())
        break;
      std::cout
          << "Val: state   "
          << "k   Rinf its    eps       En (au)      En (/cm)   En (/cm)\n";
      val = true;
    }
  }

  std::cout << "\n Total time: " << timer.reading_str() << "\n";

  bool run_test = false;
  if (run_test) {
    std::cout << "Test orthonormality [should all read 0]:\n";
    std::cout << "       ";
    for (auto &psi_b : wf.orbitals)
      printf("   %1i %2i  ", psi_b.n, psi_b.k);
    std::cout << "\n";
    for (auto &psi_a : wf.orbitals) {
      printf("%1i %2i  ", psi_a.n, psi_a.k);
      for (auto &psi_b : wf.orbitals) {
        if (psi_b > psi_a)
          continue;
        if (psi_a.k != psi_b.k) {
          std::cout << "         ";
          continue;
        }
        double xo = wf.radialIntegral(psi_a, psi_b);
        if (psi_a.n == psi_b.n)
          xo -= 1;
        printf(" %7.0e ", xo);
      }
      std::cout << "\n";
    }
  }

  bool print_wfs = false;
  if (print_wfs) {
    std::ofstream of("hf-orbitals.txt");
    of << "r ";
    for (auto &psi : wf.orbitals)
      of << "\"" << psi.symbol(true) << "\" ";
    of << "\n";
    for (size_t i = 0; i < wf.rgrid.ngp; i++) {
      of << wf.rgrid.r[i] << " ";
      for (auto &psi : wf.orbitals)
        of << psi.f[i] << " ";
      of << "\n";
    }
  }

  // for (size_t i = 0; i < 20; i++) {
  //   printf("%.2e:  %.2e %.2e  :  %.2e %.2e\n", wf.rgrid.r[i],
  //          wf.orbitals[17].f[i], wf.orbitals[17].g[i], wf.orbitals[19].f[i],
  //          wf.orbitals[19].g[i]);
  // }

  bool testpnc = false;
  if (testpnc) {
    double t = 2.3;
    double c = Nucleus::approximate_rc(wf.Anuc());

    auto rho = Nucleus::fermiNuclearDensity(1, t, c, wf.rgrid);

    double Gf = FPC::GFe11;
    double Cc = (Gf / sqrt(8.)) * (-wf.Nnuc()); // Qw/(-N)

    double Ac = 2. / 6.; // angular coef
    auto a6s_i = wf.getStateIndex(6, -1);
    auto a7s_i = wf.getStateIndex(7, -1);
    std::cout << a6s_i << "," << a7s_i << "\n";
    auto &a6s = wf.orbitals[a6s_i];
    auto &a7s = wf.orbitals[a7s_i];
    double pnc = 0;
    for (auto np : wf.orbitals) {
      if (np.k != 1)
        continue; // p_1/2 only
      int n = np.n;
      // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
      double d7s = wf.radialIntegral(a7s, np, wf.rgrid.r);
      double w6s = wf.radialIntegral(np, a6s, rho, Operator::gamma5);
      double dE6s = a6s.en - np.en;
      double d6s = wf.radialIntegral(np, a6s, wf.rgrid.r);
      double w7s = wf.radialIntegral(a7s, np, rho, Operator::gamma5);
      double dE7s = a7s.en - np.en;
      double pnc1 = Cc * Ac * d7s * w6s / dE6s;
      double pnc2 = Cc * Ac * d6s * w7s / dE7s;
      std::cout << "n=" << n << " pnc= " << pnc1 << " + " << pnc2 << " = "
                << pnc1 + pnc2 << "\n";
      pnc += pnc1 + pnc2;
    }
    std::cout << "Total= " << pnc << "\n";
    std::cout << "\n Total time: " << timer.reading_str() << "\n";
  }

  return 0;
}
