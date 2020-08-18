#include "Modules/Module_testFeynman.hpp"
#include "Coulomb/Coulomb.hpp"
#include "IO/UserInput.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <vector>

namespace Module {
void testFeynman(const IO::UserInputBlock &input, const Wavefunction &wf) {
  std::cout << "\ntestFeynman:\n";

  input.checkBlock({"real_omega", "max_l", "pol_basis", "screening", "rmin",
                    "rmax", "stride", "include_G"});

  const double omre = input.get("real_omega", -0.2);
  const auto method = MBPT::Method::Feynman;
  const auto min_n_core = 1;                       // not used?
  const int max_l_excited = input.get("max_l", 4); // up to g

  const auto include_G = input.get("include_G", false); // up to g

  const MBPT::Sigma_params sigp{method, min_n_core, max_l_excited, false,
                                true,   omre,       false,         include_G};

  const auto rmin = input.get("rmin", 1.0e-4);
  const auto rmax = input.get("rmax", 30.0);
  const auto stride = input.get("stride", std::size_t(4));
  const MBPT::rgrid_params gridp{rmin, rmax, stride};

  const MBPT::FeynmanSigma Sigma(wf.getHF(), wf.basis, sigp, gridp, {}, "NA");

  std::cout << "\n***********************************\n\n";

  // this is just to check Pol.
  // Should be zero?
  std::cout << "Testing polarisation operator, Int[ pi(1,2) d1 d2]:\n";
  std::cout << "ka kA   Int(pi)      Expected     eps\n";
  for (const auto ka : {-1, 1, -2, 2, -3}) {
    for (const auto kA : {-1, 1, -2, 2, -3}) {
      // std::cout << "ka=" << ka << " kalpha=" << kA << "\n";
      auto pi_aA = Sigma.Polarisation(ka, kA, omre, 0.0, MBPT::GrMethod::basis);
      // auto rePi = pi_aA.get_real();
      auto imPi = pi_aA.get_imaginary();
      // double sum_re = 0.0;
      double sum_im = 0.0;
      for (auto i = 0ul; i < pi_aA.size; ++i) {
        const auto dri = Sigma.dr_subToFull(i);
        for (auto j = 0ul; j < pi_aA.size; ++j) {
          const auto drj = Sigma.dr_subToFull(j);
          // sum_re += rePi.ff[i][j] * dri * drj;
          sum_im += imPi.ff[i][j] * dri * drj;
        }
      }
      // std::cout << "Int(Pi) = " << sum_im << "\n";
      // std::cout << sum_re << " + " << sum_im << "i\n";

      double expected = 0.0;
      for (const auto &Fa : Sigma.m_holes) {
        if (Fa.k != ka)
          continue;
        for (const auto &FA : Sigma.m_excited) {
          if (FA.k != kA)
            continue;
          const auto me = (Fa * FA);
          const auto eaA = (Fa.en - FA.en);
          const auto denom = eaA * eaA - omre * omre; // add imag!
          expected += me * me * 2.0 * eaA / denom;
        }
      }
      const auto eps = std::abs((sum_im - expected) / (expected + sum_im));

      printf("%2i %2i| %11.4e  %11.4e  %7.1e", ka, kA, sum_im, expected, eps);
      if ((eps > 1.0e-3 && std::abs(expected) > 1.0e-6) ||
          (std::abs(expected) < 1.0e-6 && std::abs(sum_im) > 1.0e-3))
        std::cout << " ** ";
      std::cout << "\n";
    }
  }
  std::cin.get();

  //----------------------------------------------------------------------------
  // Just for testing Vx matrix.
  std::cout << "Test Vx matirx (also, q^k):\n";
  {
    double worst = 0.0;
    const DiracSpinor *Fworst = &wf.core.front();
    for (const auto orbs : {&wf.core, &wf.valence}) {
      for (const auto &a : *orbs) {
        const auto &Vx = Sigma.get_Vx_kappa(a.k);
        const auto vxmat = Sigma.Sigma_G_Fv_2(a, Vx, a);
        const auto vxhf = a * (wf.getHF()->calc_vexFa(a));
        const auto eps = std::abs((vxmat - vxhf) / vxhf);
        if (eps > worst) {
          worst = eps;
          Fworst = &a;
        }
        printf("%4s %.4e %.4e | %.1e\n", a.shortSymbol().c_str(), vxmat, vxhf,
               eps);
      }
    }
    std::cout << "Worst: " << Fworst->symbol() << " eps=" << worst << "\n";
  }
  std::cout << "\n";

  //----------------------------------------------------------------------------
  // Just for testing Vx matrix.
  std::cout << "Test Q matirx: <aa|q^k|aa> = R^k_aaaa\n";
  {
    for (const auto orbs : {/*&wf.core,*/ &wf.valence}) {
      for (const auto &a : *orbs) {
        double worst = 0.0;
        int worstk = -1;
        for (int k = 0; k <= Sigma.maxk(); ++k) {

          const auto &qk = Sigma.get_qk(k);
          const auto aQa =
              (Sigma.G_single(a, a, 1.0)).mult_elements_by(qk).get_real();
          const auto q_aa = Sigma.Sigma_G_Fv_2(a, aQa, a);
          const auto Rk = Coulomb::Rk_abcd(a, a, a, a, k);
          const auto eps = std::abs((q_aa - Rk) / Rk);
          if (eps > worst) {
            worst = eps;
            worstk = k;
          }

          printf("k=%2i %4s %.4e %.4e | %.1e\n", k, a.shortSymbol().c_str(),
                 q_aa, Rk, eps);
        }

        std::cout << "Worst (" << a.shortSymbol() << ") : k=" << worstk
                  << " eps=" << worst << "\n\n";
      }
    }
  }
  std::cout << "\n";
  std::cin.get();

  //----------------------------------------------------------------------------
  std::cout << "Testing Gr(ev + w)\n";
  std::cout << "Comparing <v|G(e)|v> = sum_n <v|n><n|v>/(e-en) = 1/(e-en)\n";
  const double env = -0.127;
  for (const auto &omim : {0.0, 0.2, 10.0}) {

    std::cout << "\nev=" << env << ", w=" << omre << " + " << omim << "i\n";
    std::cout << "       Gr           Basis        Expected   |  Gr/Bs  Gr/Ex  "
                 "Bs/Ex\n";
    for (auto kappa : {-1, 1, -2 /*, 2, -3*/}) {
      std::cout << "kappa: " << kappa << "\n";

      const auto Gr1 = Sigma.Green(kappa, env + omre, omim, MBPT::States::both,
                                   MBPT::GrMethod::Green);
      const auto Gr2 = Sigma.Green(kappa, env + omre, omim, MBPT::States::both,
                                   MBPT::GrMethod::basis);

      // ????????????????
      auto zero1 = (Gr1 * Gr1.inverse()).plusIdent(-1.0);
      auto zero2 = (Gr2 * Gr2.inverse()).plusIdent(-1.0);
      // needs dr ? no ?
      auto [a, b] = zero1.max_el();
      auto [c, d] = zero2.max_el();
      std::cout << "Inverse (Gr): " << a << " " << b << "\n";
      std::cout << "Inverse (Basis): " << c << " " << d << "\n";
      // ????????????????

      for (const auto orbs : {&wf.core, &wf.valence}) {
        for (const auto &Fv : *orbs) {
          if (Fv.k != kappa)
            continue;
          auto vGv1 = Fv * Sigma.Sigma_G_Fv(Gr1.get_real(), Fv);
          auto vGv2 = Fv * Sigma.Sigma_G_Fv(Gr2.get_real(), Fv);
          // auto expect = 1.0 / (omre - Fv.en);

          auto denom = (omre - Fv.en) * (omre - Fv.en) - omim * omim;
          auto expect_re = (omre - Fv.en) / denom;
          auto expect_im = -omim / denom;

          auto error0 = std::abs((vGv1 - vGv2) / (vGv1 + vGv2));
          auto error1 = std::abs((vGv1 - expect_re) / expect_re);
          auto error2 = std::abs((vGv2 - expect_re) / expect_re);
          printf("%4s  %11.4e  %11.4e  %11.4e |  %.0e  %.0e  %.0e\n",
                 Fv.shortSymbol().c_str(), vGv1, vGv2, expect_re, error0,
                 error1, error2);

          if (omim != 0.0) {
            auto vGv1_i = Fv * Sigma.Sigma_G_Fv(Gr1.get_imaginary(), Fv);
            auto vGv2_i = Fv * Sigma.Sigma_G_Fv(Gr2.get_imaginary(), Fv);
            auto error0_i = std::abs((vGv1_i - vGv2_i) / (vGv1_i + vGv2_i));
            auto error1_i = std::abs((vGv1_i - expect_im) / expect_im);
            auto error2_i = std::abs((vGv2_i - expect_im) / expect_im);
            printf("      %+11.4ei %+11.4ei %+11.4ei|  %.0e  %.0e  %.0e\n",
                   vGv1_i, vGv2_i, expect_im, error0_i, error1_i, error2_i);
          }
        } // Fv
      }   // orbs (core/val)
    }     // kappa
  }       // omim

  //----------------------------------------------------------------------------

  /*

    1. Test Q vs. yk
    2. Test Vexchange
    3. Test inversion of real and imaginary Gmatrix
    4. Test Gr at real real + Imaginary omega
    5. Test Gr_core and Gr_excited
    6. Test Pol operator somehow?
    7. Test w integration somehow
        7.1 Omega grid: just integrate some function over grid
        7.2 Maybe pi(e) = \int (dw/2pi) G(w)*G(w+e)

  */

  // //------------------------------------------------------------------------------
  // void FeynmanSigma::sumPol(const ComplexGMatrix &pi_aA) const {
  //   // this is just to check Pol
  //   auto rePi = pi_aA.get_real();
  //   auto imPi = pi_aA.get_imaginary();
  //   double sum_re = 0.0;
  //   double sum_im = 0.0;
  //   for (auto i = 0ul; i < m_subgrid_points; ++i) {
  //     const auto dri = dr_subToFull(i);
  //     for (auto j = 0ul; j < m_subgrid_points; ++j) {
  //       const auto drj = dr_subToFull(j);
  //       sum_re += rePi.ff[i][j] * dri * drj;
  //       sum_im += imPi.ff[i][j] * dri * drj;
  //     }
  //   }
  //   std::cout << sum_re << "+" << sum_im << "i\n";
  // }

} // namespace Module
} // namespace Module
