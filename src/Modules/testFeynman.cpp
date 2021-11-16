#include "Modules/testFeynman.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <complex>
#include <iostream>
#include <vector>

using ComplexDouble = std::complex<double>;

namespace Module {

void testFeynman(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\ntestFeynman:\n";

  input.checkBlock_old({"real_omega", "max_l", "screening", "rmin", "rmax",
                        "stride", "include_G", "testQ", "testGreen", "testGQ",
                        "testPol", "omim"});

  const double omre = input.get("real_omega", -0.2);
  const auto method = MBPT::Method::Feynman;
  const auto min_n_core = 3;                       // Used in Polarisation
  const int max_l_excited = input.get("max_l", 4); // up to g

  const auto include_G = input.get("include_G", false); // up to g

  double w0 = 0.01;
  double wratio = 1.5;

  const MBPT::Sigma_params sigp{method, min_n_core, include_G, max_l_excited,
                                false,  true,       omre,      w0,
                                wratio, false,      false};

  const auto rmin = input.get("rmin", 1.0e-4);
  const auto rmax = input.get("rmax", 30.0);
  const auto stride = input.get("stride", std::size_t(4));
  const MBPT::rgrid_params gridp{rmin, rmax, stride};

  const MBPT::FeynmanSigma Sigma(wf.getHF(), wf.basis, sigp, gridp, "NA");

  Sigma.print_subGrid();

  const auto testQ = input.get("testQ", true);
  const auto testGreen = input.get("testGreen", false);
  const auto testGQ = input.get("testGQ", true);
  const auto testPol = input.get("testPol", false);

  const auto omim_v = input.get("omim", std::vector{0.0, 1.0, 10.0, 100.0});

  //----------------------------------------------------------------------------
  // 1. Test Q vs. yk
  // 2. Test Vexchange
  if (testQ)
    Feyn::test_Q(wf, Sigma);

  // 3. Test inversion of real and imaginary Gmatrix
  // 4. Test Gr at real real + Imaginary omega
  // 5. Test Gr_core and Gr_excited
  if (testGreen)
    Feyn::test_green(wf, Sigma, omre, omim_v);

  if (testGQ)
    Feyn::test_GQ(wf, Sigma, omre, omim_v);

  // 6. Test Pol operator
  if (testPol)
    Feyn::test_pol(wf, Sigma, omre, omim_v, max_l_excited);

  //----------------------------------------------------------------------------

  /*
    7. Test w integration somehow
        7.1 Omega grid: just integrate some function over grid
        7.2 Maybe pi(e) = \int (dw/2pi) G(w)*G(w+e)
  */
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
void Feyn::test_Q(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma) {

  // Just for testing Vx matrix.
  std::cout << "\n***********************************\n\n";
  std::cout << "Test Vx matirx (also, q^k):\n";

  {
    std::cout << "      Vx-Matrix   <v|vx|v>   | eps\n";
    double worst = 0.0;
    const DiracSpinor *Fworst = &wf.core.front();
    for (const auto orbs : {&wf.core, &wf.valence}) {
      for (const auto &a : *orbs) {
        const auto &Vx = Sigma.get_Vx_kappa(a.k);
        const auto vxmat = Sigma.act_G_Fv_2(a, Vx, a);
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
  // For testing Q matrix.
  std::cout << "\n***********************************\n\n";
  std::cout << "Test Q matirx: <aa|q^k|aa> = R^k_aaaa\n";
  {
    std::cout << "           <aa|q^k|aa>  R^k_aaaa   | eps\n";
    for (const auto orbs : {/*&wf.core,*/ &wf.valence}) {
      for (const auto &a : *orbs) {
        double worst = 0.0;
        int worstk = -1;
        for (int k = 0; k <= Sigma.maxk(); ++k) {

          const auto &qk = Sigma.get_qk(k);
          const auto aQa =
              (Sigma.G_single(a, a, 1.0)).mult_elements_by(qk).get_real();
          const auto q_aa = Sigma.act_G_Fv_2(a, aQa, a);
          const auto Rk = Coulomb::Rk_abcd(a, a, a, a, k);
          const auto eps = std::abs((q_aa - Rk) / Rk);
          if (eps > worst) {
            worst = eps;
            worstk = k;
          }

          printf("k=%2i %4s  %.4e   %.4e | %.1e\n", k, a.shortSymbol().c_str(),
                 q_aa, Rk, eps);
        }

        std::cout << "Worst (" << a.shortSymbol() << ") : k=" << worstk
                  << " eps=" << worst << "\n\n";
      }
    }
  }
  std::cout << "\n";
}

//******************************************************************************
void Feyn::test_green(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma,
                      double omre, const std::vector<double> &omim_v) {
  //----------------------------------------------------------------------------
  std::cout << "Testing Gr(ev + w)\n";
  std::cout << "Comparing <v|G(e)|v> = sum_n <v|n><n|v>/(e-en) = 1/(e-ev)\n";
  const double env = wf.valence[0].en();

  const auto max_ki = DiracSpinor::max_kindex(wf.valence);

  for (const auto &omim : omim_v) {

    ComplexDouble omega{omre, omim};

    std::cout << "\nev=" << env << ", w=" << omre << " + " << omim << "i\n";
    for (auto ik = 0; ik <= max_ki; ik++) {
      auto kappa = Angular::kappaFromIndex(ik);
      std::cout << "kappa: " << kappa << "\n";

      const auto Gr1 = Sigma.Green(kappa, {env + omre, omim},
                                   MBPT::States::both, MBPT::GrMethod::Green);
      const auto Gr2 = Sigma.Green(kappa, {env + omre, omim},
                                   MBPT::States::both, MBPT::GrMethod::basis);

      // ????????????????
      auto zero1 = (Gr1 * Gr1.inverse()).plusIdent(-1.0);
      auto zero2 = (Gr2 * Gr2.inverse()).plusIdent(-1.0);
      // needs dr ? no ?
      auto [a, b] = zero1.max_el();
      auto [c, d] = zero2.max_el();
      std::cout << "Inverse (Gr): " << a << " " << b << "\n";
      std::cout << "Inverse (Basis): " << c << " " << d << "\n";
      // ????????????????

      std::cout
          << "       Gr           Basis        Expected   |  Gr/Bs  Gr/Ex  "
             "Bs/Ex\n";

      for (const auto orbs : {&wf.core, &wf.valence}) {
        for (const auto &Fv : *orbs) {
          if (Fv.k != kappa)
            continue;
          auto vGv1 = Fv * Sigma.act_G_Fv(Gr1.get_real(), Fv);
          auto vGv2 = Fv * Sigma.act_G_Fv(Gr2.get_real(), Fv);

          auto denom = omega + ComplexDouble{env - Fv.en(), 0};
          auto expected = 1.0 / denom;

          auto expect_re = expected.real();
          auto expect_im = expected.imag();

          auto error0 = std::abs((vGv1 - vGv2) / (vGv1 + vGv2));
          auto error1 = std::abs((vGv1 - expect_re) / expect_re);
          auto error2 = std::abs((vGv2 - expect_re) / expect_re);
          printf("%4s  %11.4e  %11.4e  %11.4e |  %.0e  %.0e  %.0e\n",
                 Fv.shortSymbol().c_str(), vGv1, vGv2, expect_re, error0,
                 error1, error2);

          if (omim != 0.0) {
            auto vGv1_i = Fv * Sigma.act_G_Fv(Gr1.get_imaginary(), Fv);
            auto vGv2_i = Fv * Sigma.act_G_Fv(Gr2.get_imaginary(), Fv);
            auto error0_i = std::abs((vGv1_i - vGv2_i) / (vGv1_i + vGv2_i));
            auto error1_i = std::abs((vGv1_i - expect_im) / expect_im);
            auto error2_i = std::abs((vGv2_i - expect_im) / expect_im);
            printf("      %+11.4ei %+11.4ei %+11.4ei|  %.0e  %.0e  %.0e\n",
                   vGv1_i, vGv2_i, expect_im, error0_i, error1_i, error2_i);
          }
        } // Fv
      }   // orbs (core/val)
    }     // kappa
    std::cout << "\n";
  } // omim
}

//******************************************************************************
void Feyn::test_GQ(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma,
                   double omre, const std::vector<double> &omim_v) {
  //----------------------------------------------------------------------------
  std::cout << "Testing Gr(e + w)\n";
  std::cout << "Comparing <v| G*Q |v> = sum_i R^k_viiv / (e + w - ei)\n";
  const double env = wf.valence[0].en();
  const auto max_ki = DiracSpinor::max_kindex(wf.basis);

  const auto method = MBPT::GrMethod::Green;
  // const auto method = MBPT::GrMethod::basis;

  for (int k = 0; k <= Sigma.maxk(); ++k) {

    for (auto omim : omim_v) {
      const ComplexDouble en{env + omre, omim};

      std::cout << "\nk: " << k << ", Green/Expected\n";
      std::cout << "ev=" << env << ", w=" << omre << " + " << omim << "i\n";
      std::cout
          << "   :   Total G(e)       G_core(e)        G_ex(e)         eps\n";

      for (const auto &Fv : wf.valence) {

        ComplexDouble core = {0.0, 0.0}, ex = {0.0, 0.0};
        auto vgqv_r = 0.0, vgqCv_r = 0.0, vgqXv_r = 0.0;
        auto vgqv_i = 0.0, vgqCv_i = 0.0, vgqXv_i = 0.0;

        for (auto ik = 0; ik <= max_ki; ik++) {
          const auto kappa = Angular::kappaFromIndex(ik);

          const auto &qk = Sigma.get_qk(k);

          const auto gq =
              Sigma.Green(kappa, {env + omre, omim}, MBPT::States::both, method)
                  .mult_elements_by(qk);

          const auto gqC =
              Sigma.Green(kappa, {env + omre, omim}, MBPT::States::core, method)
                  .mult_elements_by(qk);

          const auto gqX = Sigma
                               .Green(kappa, {env + omre, omim},
                                      MBPT::States::excited, method)
                               .mult_elements_by(qk);

          vgqv_r += Sigma.act_G_Fv_2(Fv, gq.get_real(), Fv);
          vgqCv_r += Sigma.act_G_Fv_2(Fv, gqC.get_real(), Fv);
          vgqXv_r += Sigma.act_G_Fv_2(Fv, gqX.get_real(), Fv);

          vgqv_i += Sigma.act_G_Fv_2(Fv, gq.get_imaginary(), Fv);
          vgqCv_i += Sigma.act_G_Fv_2(Fv, gqC.get_imaginary(), Fv);
          vgqXv_i += Sigma.act_G_Fv_2(Fv, gqX.get_imaginary(), Fv);

          for (const auto &Fi : wf.basis) {
            if (Fi.k != kappa)
              continue;
            auto Rk = Coulomb::Rk_abcd(Fv, Fi, Fi, Fv, k);
            auto denom = (en - ComplexDouble{Fi.en()});
            ComplexDouble term = Rk / denom;
            if (wf.isInCore(Fi.n, Fi.k))
              core += term;
            else
              ex += term;
          }
        }
        auto eps_r = vgqv_r / (core + ex).real() - 1.0;
        auto eps_i = vgqv_i / (core + ex).imag() - 1.0;
        printf("%3s:  %7.4f/%7.4f  %7.4f/%7.4f  %7.4f/%7.4f  %.0e\n",
               Fv.shortSymbol().c_str(), vgqv_r, (core + ex).real(), vgqCv_r,
               (core).real(), vgqXv_r, (ex).real(), eps_r);
        if (omim != 0.0)
          printf("  i:  %7.4f/%7.4f  %7.4f/%7.4f  %7.4f/%7.4f  %.0e\n", vgqv_i,
                 (core + ex).imag(), vgqCv_i, (core).imag(), vgqXv_i,
                 (ex).imag(), eps_i);
      }
    }
  }
  std::cout << "\n";
}

//******************************************************************************
void Feyn::test_pol(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma,
                    double omre, const std::vector<double> &omim_v,
                    int max_l_excited) {

  for (auto test_pol_method : {MBPT::GrMethod::Green, MBPT::GrMethod::basis}) {

    const std::string str_meth = test_pol_method == MBPT::GrMethod::Green ?
                                     "Green fn method" :
                                     "basis method";

    std::cout << "\n***********************************\n\n";

    std::cout
        << "Testing polarisation operator: <v|G*Q*Pi*Q|v> = C * [R^k]^2\n";

    for (auto omim : omim_v) {
      const ComplexDouble om{omre, omim};
      std::cout << "\n" << str_meth << "\n";
      std::cout << "w = " << omre << " + " << omim << " i\n";
      std::cout << "       k   vGQPQv       [R^k]^2      eps\n";
      for (int k = 0; k <= Sigma.maxk(); ++k) {

        const auto &qk = Sigma.get_qk(k);
        const auto pik = Sigma.Polarisation_k(k, {omre, omim}, test_pol_method);
        const auto qpq = (qk * pik * qk);

        // for (const auto &Fv : wf.core) {
        for (const auto &Fv : wf.valence) {
          double sum_k_gqpq = 0.0;
          double sum_k_gqpq_i = 0.0;
          ComplexDouble sum_k_RR = {0.0, 0.0};

          for (auto kB : {-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7}) {
            if (Angular::l_k(kB) > max_l_excited)
              break;
            auto f2 = Angular::Ck_kk(k, Fv.k, kB);
            if (f2 == 0)
              continue;

            const auto gqpq =
                Sigma
                    .Green(kB, {Fv.en() + omre, omim}, MBPT::States::both,
                           MBPT::GrMethod::Green)
                    .mult_elements_by(qpq);

            sum_k_gqpq += Sigma.act_G_Fv_2(Fv, gqpq.get_imaginary(), Fv);
            if (omim != 0.0)
              sum_k_gqpq_i += Sigma.act_G_Fv_2(Fv, gqpq.get_real(), Fv);

            ComplexDouble sum1 = {0.0, 0.0};
            for (const auto &FB : wf.basis) {
              if (FB.k != kB)
                continue;
              for (const auto &Fa : wf.core) {
                for (const auto &FA : wf.basis) {
                  const auto f = Angular::Ck_kk(k, Fa.k, FA.k);
                  if (f == 0.0)
                    continue;
                  const auto ide1 =
                      1.0 / (ComplexDouble{Fv.en() - FB.en()} + om);
                  const auto ide2 =
                      1.0 / (ComplexDouble{Fa.en() - FA.en()} - om);
                  const auto ide3 =
                      1.0 / (ComplexDouble{Fa.en() - FA.en()} + om);
                  const auto Dinv = ide1 * (ide2 + ide3);
                  const auto Rk = Coulomb::Rk_abcd(Fv, FA, FB, Fa, k);
                  sum1 += Rk * Rk * f * f * Dinv * (1.0 / (2 * k + 1));
                }
              }
            }
            sum_k_RR += ComplexDouble{0.0, 1.0} * sum1; // i from pi
          }
          const auto eps =
              (sum_k_gqpq - sum_k_RR.imag()) / (sum_k_gqpq + sum_k_RR.imag());
          const auto epsi = (sum_k_gqpq_i - sum_k_RR.real()) /
                            (sum_k_gqpq_i + sum_k_RR.real());
          if (sum_k_RR.imag() != 0.0 || sum_k_gqpq != 0.0) {
            printf("%6s%2i| %11.4e  %11.4e  %6.0e\n", Fv.symbol().c_str(), k,
                   sum_k_gqpq, sum_k_RR.imag(), eps);
          }
          if (sum_k_RR.real() != 0.0 || sum_k_gqpq_i != 0.0)
            printf("%6s  | %11.4ei %11.4ei %6.0e\n", Fv.symbol().c_str(),
                   sum_k_gqpq_i, sum_k_RR.real(), epsi);
        }
        std::cout << "\n";
      }
    }

    //

    std::cout << "***********************************\n\n";
    std::cout << "Testing polarisation operator, Int[ pi(1,2) d1 d2]:\n";
    for (auto omim : omim_v) {
      std::cout << "\n" << str_meth << "\n";
      std::cout << "w = " << omre << " + " << omim << " i\n";
      std::cout << " k    Int(pi)      Expected     delta\n";
      for (int k = 0; k <= Sigma.maxk(); ++k) {
        auto pi = Sigma.Polarisation_k(k, {omre, omim}, test_pol_method);
        ComplexDouble om = {omre, omim};
        const auto rePi = pi.get_real();
        const auto imPi = pi.get_imaginary();

        double sum_re = 0.0;
        double sum_im = 0.0;
        for (auto i = 0ul; i < pi.size; ++i) {
          const auto dri = Sigma.dr_subToFull(i);
          for (auto j = 0ul; j < pi.size; ++j) {
            const auto drj = Sigma.dr_subToFull(j);
            sum_re += imPi.ff[i][j] * dri * drj;
            sum_im += rePi.ff[i][j] * dri * drj;
          }
        }

        ComplexDouble expected = {0.0, 0.0};
        for (const auto &Fa : Sigma.m_holes) {
          for (const auto &FA : Sigma.m_excited) {
            const auto f = Angular::Ck_kk(k, Fa.k, FA.k);
            if (f == 0.0)
              continue;
            const auto me = (Fa * FA);
            const auto ide2 = 1.0 / (ComplexDouble{Fa.en() - FA.en()} - om);
            const auto ide3 = 1.0 / (ComplexDouble{Fa.en() - FA.en()} + om);
            const auto Dinv = (ide2 + ide3);
            expected += me * me * f * f * Dinv * (1.0 / (2 * k + 1));
          }
        }
        expected *= ComplexDouble{0.0, 1.0}; // i from pi
        const auto del = (sum_re - expected.imag());
        const auto del2 = (sum_im - expected.real());

        if (sum_re != 0.0 || expected.imag() != 0.0)
          printf("%2i | %11.4e  %11.4e  %8.1e\n", k, sum_re, expected.imag(),
                 del);
        if (sum_im != 0.0 || expected.real() != 0.0)
          printf("   | %11.4ei %11.4ei %8.1e\n", sum_im, expected.real(), del2);
      }
    }
  }
}

} // namespace Module
