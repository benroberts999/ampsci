#include "ATI_atomInfo.h"
#include "ChronoTimer.h"
#include "ElectronOrbitals.h"
#include "FileIO_fileReadWrite.h"
#include "INT_quadratureIntegration.h"
#include <cmath>
#include <iostream>
#include <tuple>

int main(int argc, char *argv[]) {
  ChronoTimer sw(true); // start the overall timer

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
      double eng = wf.diracen(Z, n, k);
      wf.solveLocalDirac(n, k, eng);
    }
  }

  printf("Grid: pts=%i h=%7.5f; r0=%.1e, Rmax=%5.1f\n", wf.ngp, wf.h,
         wf.r.front(), wf.r.back());
  if (varalpha != 1)
    std::cout << "varalpha = c/c_eff = " << varalpha << " ";
  if (varalpha < 1)
    std::cout << "(non-relativistic scenario)\n";
  if (varalpha > 1)
    std::cout << "(hyper-relativistic scenario)\n";

  std::cout << "\n";

  std::cout << " n l_j    k  R_inf its eps     En (au)            Error (au)\n";
  for (auto i : wf.stateIndexList) {
    int n = wf.n_pqn(i);
    int k = wf.ka(i);
    double del = wf.en[i] - wf.diracen(wf.Znuc(), n, k);
    double rinf = wf.rinf(i);
    printf("%7s (%2i)  %3.0f %3i  %5.0e  %.15f  %7.0e\n",
           wf.seTermSymbol(i).c_str(), k, rinf, wf.itslist[i], wf.epslist[i],
           wf.en[i], del);
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
    for (auto s : wf.stateIndexList) {
      printf("%7s : ", wf.seTermSymbol(s).c_str());
      for (int in = -2; in <= 2; in++) {
        if (in == 0)
          continue;
        std::vector<double> rton;
        rton.reserve(wf.ngp);
        for (auto r : wf.r)
          rton.push_back(r);
        double R1 = wf.radialIntegral(s, s, rton);
        printf("%13.8f, ", R1);
      }
      std::cout << "\n";
    }

    // Testing Dirac Eq. by evaluating <a|H|a> - ME of Hamiltonian
    std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";
    double alpha = wf.get_alpha();
    double a2 = pow(alpha, 2);
    for (auto s : wf.stateIndexList) {
      std::vector<double> dQ(wf.ngp);
      INT::diff(wf.g[s], wf.drdt, wf.h, dQ);
      std::vector<double> rad;
      for (int i = 0; i < wf.ngp; i++) {
        double x1 = -2 * wf.f[s][i] * dQ[i] / alpha;
        double x2 =
            2 * wf.kappa[s] * wf.f[s][i] * wf.g[s][i] / (wf.r[i] * alpha);
        double x3 = -2 * pow(wf.g[s][i], 2) / a2;
        double x4 = wf.vnuc[i] * (pow(wf.f[s][i], 2) + pow(wf.g[s][i], 2));
        rad.push_back(x1 + x3 + x2 + x4);
      }
      double R = INT::integrate2(rad, wf.drdt) * wf.h;
      double fracdiff = (R - wf.en[s]) / wf.en[s];
      printf("<%i% i|H|%i% i> = % .15f, E(%i% i) = % .15f; % .0e\n",
             wf.n_pqn(s), wf.ka(s), wf.n_pqn(s), wf.ka(s), R, wf.n_pqn(s),
             wf.ka(s), wf.en[s], fracdiff);
    }
  }

  std::cout << "\n Total time: " << sw.reading_str() << "\n";

  return 0;
}
