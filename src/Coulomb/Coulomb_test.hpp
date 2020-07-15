#pragma once
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <string>

namespace UnitTest {

namespace helper {

//******************************************************************************
std::vector<double> yk_naive(const DiracSpinor &Fa, const DiracSpinor &Fb,
                             int k) {
  // Naiive (slow, but simple) implementation of yk
  const auto &gr = *Fa.rgrid;
  std::vector<double> yk(gr.r.size());
#pragma omp parallel for
  for (auto i = 0ul; i < yk.size(); ++i) {
    auto r = gr.r[i];

    auto rtkr = [](double x, double y, int kk) {
      return x < y ? std::pow(x / y, kk) / y : std::pow(y / x, kk) / x;
    };

    std::vector<double> f;
    f.reserve(gr.r.size());
    for (auto j = 0ul; j < yk.size(); ++j) {
      auto rp = gr.r[j];
      f.push_back(rtkr(r, rp, k) * (Fa.f[j] * Fb.f[j] + Fa.g[j] * Fb.g[j]));
    }

    auto p0 = 0ul; // std::max(Fa.p0, Fb.p0);
    auto pi = 0ul; // std::min(Fa.pinf, Fb.pinf);
    yk[i] = NumCalc::integrate(gr.du, p0, pi, f, gr.drdu);
  }

  return yk;
}

//******************************************************************************
std::vector<double> yk_naive(const std::vector<double> &v, const Grid &gr,
                             int k) {
  // THIS WORKS - checked against MMA
  // Naiive (slow, but simple) implementation of yk
  std::vector<double> yk(gr.r.size());
#pragma omp parallel for
  for (auto i = 0ul; i < yk.size(); ++i) {
    auto r = gr.r[i];

    auto rtkr = [](double x, double y, int kk) {
      return x < y ? std::pow(x / y, kk) / y : std::pow(y / x, kk) / x;
    };

    std::vector<double> f;
    f.reserve(gr.r.size());
    for (auto j = 0ul; j < yk.size(); ++j) {
      auto rp = gr.r[j];
      f.push_back(rtkr(r, rp, k) * v[j]);
    }

    auto p0 = 0ul; // std::max(Fa.p0, Fb.p0);
    auto pi = 0ul; // std::min(Fa.pinf, Fb.pinf);
    yk[i] = NumCalc::integrate(gr.du, p0, pi, f, gr.drdu);
  }

  return yk;
}
//
// //------------------------------------------------------------------------------
// std::vector<double> yk_naive2(const DiracSpinor &Fa, const DiracSpinor &Fb,
//                               int k) {
//   // Naiive (slow, but simple) implementation of yk
//   const auto &gr = *Fa.rgrid;
//
//   std::vector<double> yk(gr.r.size());
//   // #pragma omp parallel for
//   for (auto i = 0ul; i < yk.size(); ++i) {
//     auto r = gr.r[i];
//
//     std::vector<double> rp1, rp2;
//     rp1.reserve(gr.r.size());
//     rp2.reserve(gr.r.size());
//     for (auto rp : gr.r) {
//       rp1.push_back(std::pow(rp / r, k) / r);
//       rp2.push_back(std::pow(r / rp, k) / rp);
//     }
//
//     yk[i] = (                                                             //
//                 NumCalc::integrate(1.0, 0, i, rp1, Fa.f, Fb.f, gr.drdu)   //
//                 + NumCalc::integrate(1.0, 0, i, rp1, Fa.g, Fb.g, gr.drdu) //
//                 + NumCalc::integrate(1.0, i, 0, rp2, Fa.f, Fb.f, gr.drdu) //
//                 + NumCalc::integrate(1.0, i, 0, rp2, Fa.g, Fb.g, gr.drdu)
//                 //
//                 ) *
//             gr.du;
//   }
//
//   return yk;
// }
// std::vector<double> yk_naive2(const std::vector<double> &v, const Grid &gr,
//                               int k) {
//   // Naiive (slow, but simple) implementation of yk
//   // const auto &gr = *Fa.rgrid;
//
//   std::vector<double> yk(gr.r.size());
//   // #pragma omp parallel for
//   for (auto i = 0ul; i < yk.size(); ++i) {
//     auto r = gr.r[i];
//
//     std::vector<double> rp1, rp2;
//     rp1.reserve(gr.r.size());
//     rp2.reserve(gr.r.size());
//     for (auto rp : gr.r) {
//       rp1.push_back(std::pow(rp / r, k) / r);
//       rp2.push_back(std::pow(r / rp, k) / rp);
//     }
//
//     yk[i] = (NumCalc::integrate(gr.du, 0, i + 1, rp1, v, gr.drdu) + // works!
//              NumCalc::integrate(gr.du, i + 1, gr.r.size(), rp2, v,
//                                 gr.drdu) // doesn't?
//     );
//   }
//
//   return yk;
// }

//******************************************************************************
double check_ykab_Tab(const std::vector<DiracSpinor> &a,
                      const std::vector<DiracSpinor> &b,
                      const Coulomb::YkTable &Yab) {
  //
  double worst = 0.0;
  for (const auto &Fa : a) {
    for (const auto &Fb : b) {
      const auto [kmin, kmax] = Coulomb::YkTable::k_minmax(Fa, Fb);
      for (int k = kmin; k <= kmax; ++k) {
        // Only check if Angular factor is non-zero (since Ykab only calc'd in
        // this case)
        if (!Angular::Ck_kk_SR(k, Fa.k, Fb.k))
          continue;
        const auto &y1 = Yab.get_yk_ab(k, Fa, Fb);
        const auto y2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto y3 = Coulomb::yk_ab(Fb, Fa, k); // symmetric
        const auto del = std::abs(qip::compare(y1, y2).first) +
                         std::abs(qip::compare(y2, y3).first);
        if (del > worst)
          worst = del;
      }
    }
  }

  return worst;
}

//******************************************************************************
std::vector<double> check_ykab(const std::vector<DiracSpinor> &a,
                               const std::vector<DiracSpinor> &b) {

  // Compared Yk as calculated by fast Coulomb::yk_ab routine (used in the code)
  // to helper::yk_naive, a very slow, but simple version. In theory, should be
  // exactly the same

  std::vector<double> worst; // = 0.0;
  // nb: slow, so only check sub-set
  for (auto ia = 0ul; ia < a.size(); ia += 4) {
    const auto &Fa = a[ia];
    for (auto ib = 1ul; ib < b.size(); ib += 4) {
      const auto &Fb = b[ib];
      const auto [kmin, kmax] = Coulomb::YkTable::k_minmax(Fa, Fb);
      for (int k = kmin; k <= kmax; ++k) {
        if (!Angular::Ck_kk_SR(k, Fa.k, Fb.k))
          continue;
        if (std::size_t(k + 1) > worst.size())
          worst.resize(std::size_t(k + 1));
        const auto y2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto y4 = helper::yk_naive(Fa, Fb, k); // slow
        // const auto y5 = helper::yk_naive2(Fa, Fb, k); // slow

        // for (int i = 0; i < 1000; i += 100) {
        //   std::cout << k << " " << i << " " << y4[i] << " " << y5[i] << "\n";
        // }

        const auto del = std::abs(qip::compare(y2, y4).first);
        // const auto del = std::abs(qip::compare(y4, y5).first);
        // const auto del = std::abs(qip::compare(y2, y5).first);

        if (del > worst[std::size_t(k)])
          worst[std::size_t(k)] = del;
      }
    }
  }

  return worst;
}

} // namespace helper

//******************************************************************************
//******************************************************************************
bool Coulomb(std::ostream &obuff) {
  bool pass = true;

  { // First, test quad int:
    // Define a wavefunction-like function, func,
    auto func = [](double x) {
      return x * std::exp(-0.2 * x) * (1.0 + 2.0 * std::sin(x));
    };
    // that has an exact integral, Intfunc:
    auto Intfunc = [](double x) {
      return -(5.0 / 169.0) * std::exp(-0.2 * x) *
             (169.0 * (5.0 + x) + 5.0 * (5.0 + 13.0 * x) * std::cos(x) +
              (13.0 * x - 60.0) * std::sin(x));
    };

    // const auto pts_lst = std::vector<std::size_t>{750, 1000, 2000};
    const auto pts_lst = std::vector<std::size_t>{1000};
    for (const auto pts : pts_lst) {
      const Grid grll(1.0e-6, 100.0, pts, GridType::loglinear, 10);
      const Grid grlog(1.0e-6, 100.0, pts, GridType::logarithmic, 0);

      std::vector<double> vll, vlog;
      for (const auto &r : grll.r)
        vll.push_back(func(r));
      for (const auto &r : grlog.r)
        vlog.push_back(func(r));

      // numerical integration (on grid):
      const auto intll = NumCalc::integrate(grll.du, 0, 0, vll, grll.drdu);
      const auto intlog = NumCalc::integrate(grlog.du, 0, 0, vlog, grlog.drdu);

      // Account for possibility r0, rmax slightly different:
      const auto exactll = Intfunc(grll.r.back()) - Intfunc(grll.r.front());
      const auto exactlog = Intfunc(grlog.r.back()) - Intfunc(grlog.r.front());

      pass &= qip::check_value(
          &obuff, "Quad int (log-lin) b=10 N=" + std::to_string(pts),
          (intll - exactll) / exactll, 0.0, 1.0e-13);
      pass &= qip::check_value(
          &obuff, "Quad int (logarithmic)  N=" + std::to_string(pts),
          (intlog - exactlog) / exactlog, 0.0, 1.0e-4 * (500.0 / double(pts)));

      // Compare w/ Mathematica: Works pretty well
      // std::vector<double> yy;
      // for (const auto &r : grll.r)
      //   yy.push_back(func(r) * func(r));
      // std::cout << "r:"
      //           << " " << grll.r[100] << " " << grll.r[300] << " "
      //           << grll.r[600] << "\n";
      // printf("%.7f %.7f %.7f\n", grll.r[100], grll.r[300], grll.r[600]);
      // for (int k = 0; k < 16; ++k) {
      //   auto y = helper::yk_naive(yy, grll, k);
      //   std::cout << k << " " << y[100] << " " << y[300] << " " << y[600]
      //             << "\n";
      //   std::cin.get();
      // }
    }
  }

  {
    // Don't need dense grid:
    Wavefunction wf({1000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
    // wf.hartreeFockValence("7sp5d"); //"7sp5d4f"
    wf.formBasis({"8spdfghi", 30, 7, 0.0, 1.0e-7, 30.0, false});

    // Split basis into core/excited
    std::vector<DiracSpinor> core, excited;
    for (const auto &Fb : wf.basis) {
      if (wf.isInCore(Fb.n, Fb.k)) {
        core.push_back(Fb);
      } else {
        excited.push_back(Fb);
      }
    }

    const Coulomb::YkTable Yce(wf.rgrid, &core, &excited);
    const Coulomb::YkTable Yij(wf.rgrid, &wf.basis);

    // const auto maxtj = std::max_element(wf.basis.cbegin(), wf.basis.cend(),
    //                                     DiracSpinor::comp_j)
    //                        ->twoj();
    // const auto &Ck = Yij.Ck();
    // const Angular::SixJ sj(maxtj, maxtj);

    {
      // check the 'diff' orbitals case:
      double del1 = helper::check_ykab_Tab(core, excited, Yce);
      // check the 'same' orbitals case:
      double del2 = helper::check_ykab_Tab(wf.basis, wf.basis, Yij);
      pass &= qip::check_value(&obuff, "Yk_ab tables", std::max(del1, del2),
                               0.0, 1.0e-17);
    }

    {
      // Check Yk formula:
      // FAILS - note: Seems like THIS is what is causing all the issues!
      // (Particularly blows up when i included [large k])
      auto delks = helper::check_ykab(wf.core, wf.basis);
      int k = 0;
      for (const auto &dk : delks) {
        pass &= qip::check_value(&obuff, "Yk_ab value k=" + std::to_string(k++),
                                 dk, 0.0, 1.0e-10);
      }
    }

    // XXX Add tests for R, Q, P, W, Z etc?
  }

  return pass;
}

} // namespace UnitTest
