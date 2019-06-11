#include "ATI_atomInfo.hpp"
#include "ChronoTimer.hpp"
#include "DiracOperator.hpp"
#include "ElectronOrbitals.hpp"
#include "FPC_physicalConstants.hpp"
#include "FileIO_fileReadWrite.hpp"
#include "HartreeFockClass.hpp"
#include "Nucleus.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "Operators.hpp"
#include "PRM_localPotentials.hpp"
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

  // Generate the orbitals object:
  ElectronOrbitals wf(Z, A, ngp, r0, rmax, varalpha);
  // Print info to screen:
  if (exclude_exchange)
    std::cout << "\nRunning Hartree (excluding exchange) for ";
  else
    std::cout << "\nRunning Hartree-Fock for ";
  wf.printAtom();
  wf.printNuclearParams();
  wf.rgrid.print();
  std::cout << "********************************************************\n";

  // Solve Hartree equations for the core:
  timer.start(); // start the timer for HF
  HartreeFock hf(wf, str_core, eps_HF, exclude_exchange);
  double core_energy = hf.calculateCoreEnergy();
  std::cout << "core: " << timer.lap_reading_str() << "\n";

  // Create list of valence states to solve for
  if ((int)wf.Ncore() >= wf.Znuc())
    num_val = 0;
  auto val_lst = wf.listOfStates_nk(num_val, l_max);

  // Solve for the valence states:
  timer.start();
  for (const auto &nk : val_lst) {
    int n = nk[0];
    int k = nk[1];
    hf.solveValence(n, k);
  }
  wf.orthonormaliseOrbitals(wf.valence_orbitals, 2);
  if (val_lst.size() > 0)
    std::cout << "Valence: " << timer.lap_reading_str() << "\n";

  // Output results:
  std::cout << "\nHartree Fock: ";
  wf.printAtom();
  bool sorted = true;
  wf.printCore(sorted);
  std::cout << "E_core = " << core_energy
            << " au;  = " << core_energy * FPC::Hartree_invcm << "/cm\n";
  wf.printValence(sorted);

  std::cout << "\n Total time: " << timer.reading_str() << "\n";

  //*********************************************************
  //               TESTS
  //*********************************************************

  bool run_test = true;
  if (run_test) {
    std::cout << "Test orthonormality [log-scale, should all read 0]:\n";

    for (int i = 0; i < 3; i++) {
      auto tmp_b = (i == 2) ? wf.valence_orbitals : wf.core_orbitals;
      auto tmp_a = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
      // Core-Core:
      if (i == 0)
        std::cout << "\nCore-Core\n    ";
      else if (i == 1)
        std::cout << "\nValence-Core\n    ";
      else
        std::cout << "\nValence-Valence\n    ";
      for (auto &psi_b : tmp_b)
        printf("%2i%2i", psi_b.n, psi_b.k);
      std::cout << "\n";
      for (auto &psi_a : tmp_a) {
        printf("%2i%2i", psi_a.n, psi_a.k);
        for (auto &psi_b : tmp_b) {
          if (psi_b > psi_a)
            continue;
          if (psi_a.k != psi_b.k) {
            std::cout << "    ";
            continue;
          }
          double xo = (psi_a * psi_b);
          if (psi_a.n == psi_b.n)
            xo -= 1;
          // xo = xo == 0 ? 0 : ;
          // printf(" e%+3.0f", xo);
          if (xo == 0)
            printf("   0");
          else
            printf(" %+3.0f", log10(fabs(xo)));
        }
        std::cout << "\n";
      }
    }

    std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";
    double c = 1. / wf.get_alpha();
    DiracOperator w(c, GammaMatrix::g5, 1, true);
    RadialOperator x_a(wf.rgrid, -1);
    DiracOperator y(c * c, DiracMatrix(0, 0, 0, -2));
    DiracOperator z1(wf.vnuc);
    DiracOperator z2(wf.vdir);
    for (int i = 0; i < 2; i++) {
      auto &tmp_orbs = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
      bool valence = (i == 0) ? false : true;
      for (auto &psi : tmp_orbs) {
        auto k = psi.k;
        DiracOperator z3(hf.get_vex(psi, valence));
        DiracOperator x_b(c, DiracMatrix(0, 1 - k, 1 + k, 0), 0, true);
        auto rhs = (w * psi) + (x_a * (x_b * psi)) + (y * psi) + (z1 * psi) +
                   (z2 * psi) + (z3 * psi);
        double R = psi * rhs;
        double ens = psi.en;
        double fracdiff = (R - ens) / ens;
        printf("<%i% i|H|%i% i> = %17.11f, E = %17.11f; % .0e\n", psi.n, psi.k,
               psi.n, psi.k, R, ens, fracdiff);
      }
    }
  }

  bool print_wfs = false;
  if (print_wfs) {
    std::ofstream of("hf-orbitals.txt");
    of << "r ";
    for (auto &psi : wf.core_orbitals)
      of << "\"" << psi.symbol(true) << "\" ";
    for (auto &psi : wf.valence_orbitals)
      of << "\"" << psi.symbol(true) << "\" ";
    of << "\n";
    for (std::size_t i = 0; i < wf.rgrid.ngp; i++) {
      of << wf.rgrid.r[i] << " ";
      for (auto &psi : wf.core_orbitals)
        of << psi.f[i] << " ";
      for (auto &psi : wf.valence_orbitals)
        of << psi.f[i] << " ";
      of << "\n";
    }
  }

  bool testpnc = false;
  if (testpnc) {
    double t = 2.3;
    double c = Nucleus::approximate_c_hdr(wf.Anuc());
    PNCnsiOperator hpnc(c, t, wf.rgrid, -wf.Nnuc());
    E1Operator he1(wf.rgrid);

    double Ac = 2. / 6.; // angular coef
    auto a6s_i = wf.getStateIndex(6, -1);
    auto a7s_i = wf.getStateIndex(7, -1);
    auto &a6s = wf.valence_orbitals[a6s_i];
    auto &a7s = wf.valence_orbitals[a7s_i];
    std::cout << "E_pnc: " << wf.Anuc() << "-" << ATI::atomicSymbol(wf.Znuc())
              << " " << a6s.symbol() << " -> " << a7s.symbol() << "\n";

    double pnc = 0;
    for (int i = 0; i < 2; i++) {
      auto &tmp_orbs = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
      for (auto &np : tmp_orbs) {
        if (np.k != 1)
          continue; // p_1/2 only
        // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
        double pnc1 =
            Ac * (a7s * (he1 * np)) * (np * (hpnc * a6s)) / (a6s.en - np.en);
        double pnc2 =
            Ac * (a7s * (hpnc * np)) * (np * (he1 * a6s)) / (a7s.en - np.en);
        std::cout << "n=" << np.n << " pnc= " << pnc1 << " + " << pnc2 << " = "
                  << pnc1 + pnc2 << "\n";
        pnc += pnc1 + pnc2;
      }
    }
    std::cout << "Total= " << pnc << "\n";
    std::cout << "\n Total time: " << timer.reading_str() << "\n";
  }

  bool test_hfs = false;
  if (test_hfs) {
    // Test hfs and Operator
    double muN = 2.751818;            // XXX Rb
    double IN = (3. / 2.);            // XXX Rb
    auto r_rms = 4.1989 / FPC::aB_fm; // XXX Rb
    // auto r_rms = Nucleus::approximate_r_rms(wf.Anuc());
    std::cout << "Gridpoints below Rrms: " << wf.rgrid.getIndex(r_rms) << "\n";

    // example for using lambda
    auto l1 = [](double r, double) { return 1. / (r * r); };
    // auto l2 = [](double r, double rN) { return r > rN ? 1. / (r * r) : 0.;
    // };
    HyperfineOperator vhfs(muN, IN, r_rms, wf.rgrid, l1);

    // for (auto i : wf.valenceIndexList) {
    //   auto &phi = wf.orbitals[i];
    for (auto phi : wf.valence_orbitals) {
      auto A_tmp = phi * (vhfs * phi);

      double j = phi.j();
      auto factor = FPC::Hartree_MHz * phi.k / (j * (j + 1.));

      std::cout << phi.symbol() << ": ";
      std::cout << A_tmp * factor << "\n";
    }
  }

  return 0;
}
