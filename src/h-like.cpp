#include "ATI_atomInfo.hpp"
#include "ChronoTimer.hpp"
#include "DiracOperator.hpp"
// #include "DiracSpinor.hpp"
#include "ElectronOrbitals.hpp"
#include "FileIO_fileReadWrite.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <cmath>
#include <iostream>
#include <tuple>

int main(int argc, char *argv[]) {
  ChronoTimer sw; // start the overall timer

  std::string input_file = (argc > 1) ? argv[1] : "h-like.in";
  std::cout << "Reading input from: " << input_file << "\n";

  int Z, A, n_max, l_max, ngp;
  double r0, rmax, varalpha;
  bool extra;
  {
    int iextra;
    auto tp = std::forward_as_tuple(Z, A, n_max, l_max, r0, rmax, ngp, varalpha,
                                    iextra);
    FileIO::setInputParameters(input_file, tp);
    extra = (iextra == 1) ? true : false;
  }

  printf("\nRunning SolveDBS for Local H-like potential, Z=%i\n", Z);
  printf("*************************************************\n");

  // Generate the orbitals object:
  ElectronOrbitals wf(Z, A, ngp, r0, rmax, varalpha);

  // Solve the Dirac equation for H-like ions:
  for (int n = 1; n <= n_max; n++) {
    for (int i = 1; i < 2 * n; i++) { // loop through each kappa state
      int k = int(pow(-1, i) * ceil(0.5 * i));
      int l = ATI::l_k(k);
      if (l > l_max)
        continue;
      double eng = ATI::diracen(Z, n, k, wf.get_alpha());
      wf.solveLocalDirac(n, k, eng);
    }
  }

  wf.rgrid.print();
  if (varalpha != 1)
    std::cout << "varalpha = c/c_eff = " << varalpha << " ";
  if (varalpha < 1)
    std::cout << "(non-relativistic scenario)\n";
  if (varalpha > 1)
    std::cout << "(hyper-relativistic scenario)\n";

  std::cout << "\n";

  std::cout << " n l_j    k  R_inf its eps     En (au)            Error (au)\n";
  for (auto &psi : wf.orbitals) {
    double del = psi.en - ATI::diracen(wf.Znuc(), psi.n, psi.k, wf.get_alpha());
    // wf.diracen(wf.Znuc(), psi.n, psi.k);
    double rinf = wf.rinf(psi);
    printf("%7s (%2i)  %3.0f %3i  %5.0e  %.15f  %7.0e\n", psi.symbol().c_str(),
           psi.k, rinf, psi.its, psi.eps, psi.en, del);
  }
  // wf.orthonormaliseOrbitals(2);

  if (extra) {
    // Calculate the expectation value of r^rpow for each state in list:
    std::cout << "\nExpectation value of r^n (radial integral)\n";
    std::cout << "          ";
    for (int in = -2; in <= 2; in++) {
      if (in == 0)
        continue;
      std::cout << " <nk|r^" << in << "|nk>   ";
    }
    std::cout << "\n";
    for (auto &phi : wf.orbitals) {
      printf("%7s : ", phi.symbol().c_str());
      for (int in = -2; in <= 2; in++) {
        if (in == 0)
          continue;
        std::vector<double> rton;
        rton.reserve(wf.rgrid.ngp);
        for (auto r : wf.rgrid.r)
          rton.push_back(pow(r, in));
        DiracOperator rp(1, rton);
        // double R1 = wf.radialIntegral(phi, phi, rton);
        double R1 = phi * (rp * phi);
        printf("%13.8f, ", R1);
      }
      std::cout << "\n";
    }

    // Testing Dirac Eq. by evaluating <a|H|a> - ME of Hamiltonian
    std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";
    double alpha = wf.get_alpha();
    double a2 = pow(alpha, 2);
    for (auto &psi : wf.orbitals) {
      // std::vector<double> dQ(wf.rgrid.ngp);
      // NumCalc::diff(psi.g, wf.rgrid.drdu, wf.rgrid.du, dQ);
      std::vector<double> dQ =
          NumCalc::derivative(psi.g, wf.rgrid.drdu, wf.rgrid.du);
      std::vector<double> rad;
      for (std::size_t i = 0; i < wf.rgrid.ngp; i++) {
        double x1 = -2 * psi.f[i] * dQ[i] / alpha;
        double x2 = 2 * psi.k * psi.f[i] * psi.g[i] / (wf.rgrid.r[i] * alpha);
        double x3 = -2 * pow(psi.g[i], 2) / a2;
        double x4 = wf.vnuc[i] * (pow(psi.f[i], 2) + pow(psi.g[i], 2));
        rad.push_back(x1 + x3 + x2 + x4);
      }
      double R = NumCalc::integrate(rad, wf.rgrid.drdu) * wf.rgrid.du;
      double ens = psi.en;
      double fracdiff = (R - ens) / ens;
      printf("<%i% i|H|%i% i> = % .15f, E(%i% i) = % .15f; % .0e\n", psi.n,
             psi.k, psi.n, psi.k, R, psi.n, psi.k, ens, fracdiff);
    }
  }

  std::cout << "\n Total time: " << sw.reading_str() << "\n";

  return 0;
}
