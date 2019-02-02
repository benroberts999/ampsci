#include "ATI_atomInfo.h"
#include "ChronoTimer.h"
#include "ElectronOrbitals.h"
#include "FileIO_fileReadWrite.h"
#include "INT_quadratureIntegration.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <tuple>

int main() {

  ChronoTimer sw(true); // start timer

  int Z, A, n_max, l_max, ngp;
  double r0, rmax, varalpha;
  bool extra;
  {
    int iextra;
    auto tp = std::forward_as_tuple(Z, A, n_max, l_max, r0, rmax, ngp, varalpha,
                                    iextra);
    FileIO::setInputParameters("h-like.in", tp);
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

  printf("Grid: pts=%i h=%7.5f Rmax=%5.1f\n", wf.ngp, wf.h, wf.r[wf.ngp - 1]);
  if (varalpha != 1)
    printf("varalpha = c/c_eff = %.1e  ", varalpha);
  if (varalpha < 1)
    std::cout << "(non-relativistic scenario)\n";
  if (varalpha > 1)
    std::cout << "(hyper-relativistic scenario)\n";

  std::cout << "\n";

  printf(" n l_j    k  R_inf its eps     En (au)            Error (au)\n");
  auto num_states = wf.nlist.size();
  for (auto i = 0ul; i < num_states; i++) {
    int n = wf.nlist[i];
    int k = wf.kappa[i];
    int twoj = 2 * abs(k) - 1;
    int l = (abs(2 * k + 1) - 1) / 2;
    double del = wf.en[i] - wf.diracen(wf.Znuc(), n, k);
    double rinf = wf.r[wf.pinflist[i]];
    printf("%2i %s_%i/2 (%2i)  %3.0f %3i  %5.0e  %.15f  %7.0e\n", n,
           ATI::l_symbol(l).c_str(), twoj, k, rinf, wf.itslist[i],
           wf.epslist[i], wf.en[i], del);
  }

  // wf.orthonormaliseOrbitals(2);

  if (extra) {
    // Calculate the expectation value of r^rpow for each state in list:
    printf("\nExpectation value of r^n (radial integral)\n");
    std::cout << "          ";
    for (int in = -2; in <= 2; in++) {
      if (in == 0)
        continue;
      printf(" <nk|r^%i|nk>   ", in);
    }
    std::cout << "\n";
    for (auto s : wf.stateIndexList) {
      printf("%7s : ", wf.seTermSymbol(s).c_str());
      for (int in = -2; in <= 2; in++) {
        if (in == 0)
          continue;
        std::vector<double> rad1;
        for (int i = 0; i < wf.ngp; i++) {
          double x1 = (wf.f[s][i] * wf.f[s][i] + wf.g[s][i] * wf.g[s][i]) *
                      pow(wf.r[i], in);
          rad1.push_back(x1);
        }
        double R1 = INT::integrate2(rad1, wf.drdt) * wf.h;
        printf("%13.8f, ", R1);
      }
      std::cout << "\n";
    }

    // Testing Dirac Eq. by evaluating <a|H|a> - ME of Hamiltonian
    printf("\nTesting wavefunctions: <n|H|n>  (numerical error)\n");
    double alpha = wf.get_alpha();
    double a2 = pow(alpha, 2);
    // for(auto s=0ul; s<num_states; s++){
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
             wf.nlist[s], wf.kappa[s], wf.nlist[s], wf.kappa[s], R, wf.nlist[s],
             wf.kappa[s], wf.en[s], fracdiff);
    }
  }

  std::cout << "\n Total time: " << sw.reading_str() << "\n";

  return 0;
}
