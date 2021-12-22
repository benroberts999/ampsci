#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
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
// This is a Naiive (slow, but simple + correct) implementation of yk
// This forms the "baseline" to compare against
inline std::vector<double> yk_naive(const DiracSpinor &Fa,
                                    const DiracSpinor &Fb, int k);

// This takes two sets of spinors, and calculates each y^k_ab two ways (using
// YkTable, and using Coulomb:: functions). Also checks if y_ab = y_ba.
// Compares these, returns the worst y1(r)-y2(r). This is a test of the YkTable
// and symmetry, but not the Coulomb:: formula
inline double check_ykab_Tab(const std::vector<DiracSpinor> &a,
                             const std::vector<DiracSpinor> &b,
                             const Coulomb::YkTable &Yab);

// This takes a sets of spinors, and calculates each y^k_ab two ways (using
// Coulomb:: functions, and the above 'Naiive' formula). Compares
// these, returns the worst y1(r)-y2(r).
// Calculates only those with n values that differ by less than max_del_n.
// This is a test of the Coulomb:: formulas
inline std::vector<double> check_ykab(const std::vector<DiracSpinor> &orbs,
                                      int max_del_n = 99);

// This takes in a set of orbitals, and calculates a subset of the R^k radial
// integrals a number of ways. Subset only, otherwise takes too long. Calculates
// only those with n values that differ by less than max_del_n.
// Calculates R^k several ways: Using YkTable, and a few different Coulomb::
// functions. Compares them all, and resturns worst comparison.
inline double check_Rkabcd(const std::vector<DiracSpinor> &orbs,
                           int max_del_n = 99);

} // namespace helper

//******************************************************************************
//******************************************************************************
//! Unit tests for Coulomb integrals (y^k_ab, R^k_abcd, lookup tables etc).
//! Also: tests quadrature integation method
bool Coulomb(std::ostream &obuff) {
  bool pass = true;

  { // First, test quadrature integration method:

    // Define a wavefunction-like function, func,
    auto func = [](double x) {
      return x * std::exp(-0.2 * x) * (1.0 + 2.0 * std::sin(x));
    };
    // that has an exact integral, Intfunc = \int_0^x(func):
    auto Intfunc = [](double x) {
      return -(5.0 / 169.0) * std::exp(-0.2 * x) *
             (169.0 * (5.0 + x) + 5.0 * (5.0 + 13.0 * x) * std::cos(x) +
              (13.0 * x - 60.0) * std::sin(x));
    };

    // Perform the numerical integrals using a different grids, compare to exact
    const auto pts_lst = std::vector<std::size_t>{1000};
    for (const auto pts : pts_lst) {
      const Grid grll(1.0e-6, 100.0, pts, GridType::loglinear, 10);
      const Grid grlog(1.0e-6, 100.0, pts, GridType::logarithmic, 0);

      std::vector<double> vll, vlog;
      for (const auto &r : grll.r())
        vll.push_back(func(r));
      for (const auto &r : grlog.r())
        vlog.push_back(func(r));

      // numerical integration (on grid):
      const auto intll = NumCalc::integrate(grll.du(), 0, 0, vll, grll.drdu());
      const auto intlog =
          NumCalc::integrate(grlog.du(), 0, 0, vlog, grlog.drdu());

      // Account for possibility r0, rmax slightly different:
      const auto exactll = Intfunc(grll.r().back()) - Intfunc(grll.r().front());
      const auto exactlog =
          Intfunc(grlog.r().back()) - Intfunc(grlog.r().front());

      pass &=
          qip::check_value(&obuff,
                           "Quad[" + std::to_string(NumCalc::Nquad) +
                               "] int (log-lin) b=10 N=" + std::to_string(pts),
                           (intll - exactll) / exactll, 0.0, 1.0e-13);
      pass &= qip::check_value(
          &obuff,
          "Quad[" + std::to_string(NumCalc::Nquad) +
              "] int (logarithmic)  N=" + std::to_string(pts),
          (intlog - exactlog) / exactlog, 0.0, 1.0e-4 * (500.0 / double(pts)));
    }
  }

  //----------------------------------------------------------------------------

  // Test the Coulomb formulas
  // Don't need dense grid, and use a local potential (Hartree)
  Wavefunction wf({1000, 1.0e-6, 100.0, 10.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Hartree", 0.0, "[Xe]");
  wf.formBasis({"8spdfghi", 30, 9, 1.0e-4, 1.0e-6, 40.0, false});

  // Split basis into core/excited
  std::vector<DiracSpinor> core, excited;
  for (const auto &Fb : wf.basis) {
    if (wf.isInCore(Fb.n, Fb.k)) {
      core.push_back(Fb);
    } else {
      excited.push_back(Fb);
    }
  }
  // Form the Coulomb lookup tables:
  const Coulomb::YkTable Yce(core, excited);
  const Coulomb::YkTable Yij(wf.basis);

  { // Check the Ykab lookup-tables
    // check the 'different' orbitals case:
    double del1 = helper::check_ykab_Tab(core, excited, Yce);
    // check the 'same' orbitals case:
    double del2 = helper::check_ykab_Tab(wf.basis, wf.basis, Yij);
    pass &= qip::check_value(&obuff, "Yk_ab tables", std::max(del1, del2), 0.0,
                             1.0e-17);
  }

  { // Testing the Hartree Y functions formula:
    const auto delk_core = helper::check_ykab(wf.core, 2);
    const auto delk_basis = helper::check_ykab(wf.basis, 1);
    int k = 0;
    for (const auto &dk : delk_core) {
      pass &= qip::check_value(&obuff,
                               "Yk_ab (core) value  k=" + std::to_string(k++),
                               dk, 0.0, 1.0e-14);
    }
    k = 0;
    for (const auto &dk : delk_basis) {
      pass &= qip::check_value(&obuff,
                               "Yk_ab (basis) value k=" + std::to_string(k++),
                               dk, 0.0, 1.0e-14);
    }
  }

  // test R^k_abcd:
  const double eps_R = helper::check_Rkabcd(wf.core, 2);
  const double eps_R2 = helper::check_Rkabcd(wf.basis, 1);
  pass &= qip::check_value(&obuff, "Rk_abcd (core) ", eps_R, 0.0, 1.0e-13);
  pass &= qip::check_value(&obuff, "Rk_abcd (basis) ", eps_R2, 0.0, 1.0e-13);

  //****************************************************************************
  // Test "other" Coulomb integrals: P, Q, W (these are defined in terms of sums
  // over R and angular coeficients):
  {
    // Contruct a vector of DiracSpinors, with just single spinor of each kappa:
    std::vector<DiracSpinor> torbs;
    for (int kappa_index = 0;; ++kappa_index) {
      auto k = Angular::kappaFromIndex(kappa_index);
      auto phi = std::find_if(cbegin(wf.basis), cend(wf.basis),
                              [k](auto x) { return x.k == k; });
      if (phi == cend(wf.basis))
        break;
      torbs.emplace_back(*phi);
    }

    // test Q
    const auto maxtj = std::max_element(wf.basis.cbegin(), wf.basis.cend(),
                                        DiracSpinor::comp_j)
                           ->twoj();
    const auto &Ck = Yij.Ck();
    const Angular::SixJTable sj(2 * maxtj);

    double worstQ = 0.0;
    double worstP = 0.0;
    double worstW = 0.0;
    for (const auto &Fa : torbs) {
      for (const auto &Fb : torbs) {
        for (const auto &Fc : torbs) {
          for (const auto &Fd : torbs) {
            for (int k = 0; k <= Ck.max_k(); ++k) {

              double Q1 = 0.0;
              if (Angular::Ck_kk_SR(k, Fa.k, Fc.k) &&
                  Angular::Ck_kk_SR(k, Fb.k, Fd.k)) {

                Q1 = Yij.Q(k, Fa, Fb, Fc, Fd);
                const auto Q2 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
                const auto Q3 = Fa * Yij.Qkv_bcd(Fa.k, Fb, Fc, Fd, k);

                const auto Q4 =
                    Angular::neg1pow_2(2 * k + Fa.twoj() + Fb.twoj() + 2) *
                    Ck(k, Fa.k, Fc.k) * Ck(k, Fb.k, Fd.k) *
                    Coulomb::Rk_abcd(Fa, Fb, Fc, Fd, k);

                // test the 'Qk' version, including k_minmax_Q
                const auto [kmin, kmax] = Coulomb::k_minmax_Q(Fa, Fb, Fc, Fd);
                // This k_min should have correct parity rule too
                const auto Q5 =
                    (k >= kmin && k <= kmax && (kmin % 2 == k % 2)) ?
                        Yij.Q(k, Fa, Fb, Fc, Fd) :
                        0.0;

                const auto delQ =
                    std::abs(qip::max_difference(Q1, Q2, Q3, Q4, Q5));
                if (delQ > worstQ)
                  worstQ = delQ;
              }

              // Calc P
              const auto P1 = Yij.P(k, Fa, Fb, Fc, Fd);
              const auto P2 = Coulomb::Pk_abcd(Fa, Fb, Fc, Fd, k);
              double P4 = 0.0;
              for (int l = 0; l <= Ck.max_k(); ++l) {
                if (Angular::Ck_kk_SR(l, Fa.k, Fd.k) &&
                    Angular::Ck_kk_SR(l, Fb.k, Fc.k)) {
                  P4 += (2 * k + 1) * sj.get(Fa, Fc, k, Fb, Fd, l) *
                        Yij.Q(l, Fa, Fb, Fd, Fc);
                }
              }

              // test the 'Pk' version, including k_minmax_P
              const auto [kminP, kmaxP] = Coulomb::k_minmax_P(Fa, Fb, Fc, Fd);
              // This k_min CANNOT contain correct parity rule
              const auto P5 =
                  (k >= kminP && k <= kmaxP) ? Yij.P(k, Fa, Fb, Fc, Fd) : 0.0;

              const auto delP = std::abs(qip::max_difference(P1, P2, P4, P5));
              if (delP > worstP)
                worstP = delP;

              // calc W
              const auto W1 = Q1 + P1;
              const auto W2 = Coulomb::Wk_abcd(Fa, Fb, Fc, Fd, k);

              // test the 'Wk' version, including k_minmax_P
              const auto [kmin, kmax] = Coulomb::k_minmax_W(Fa, Fb, Fc, Fd);
              // This k_min CANNOT contain correct parity rule
              const auto W5 =
                  (k >= kmin && k <= kmax) ? Yij.W(k, Fa, Fb, Fc, Fd) : 0.0;

              const auto delW = std::abs(qip::max_difference(W1, W2, W5));
              if (delW > worstW)
                worstW = delW;

            } // k
          }   // d
        }     // c
      }       // b
    }         // a
    pass &= qip::check_value(&obuff, "Qk_abcd ", worstQ, 0.0, 5.0e-13);
    pass &= qip::check_value(&obuff, "Pk_abcd ", worstP, 0.0, 5.0e-14);
    pass &= qip::check_value(&obuff, "Wk_abcd ", worstW, 0.0, 5.0e-14);

    // Test the "Magic" 6J functions
    {
      double worst = 0.0;
      for (const auto &Fa : torbs) {
        for (const auto &Fb : torbs) {
          for (const auto &Fc : torbs) {
            for (const auto &Fd : torbs) {
              for (int k = 0; k < 15; ++k) {
                for (int l = 0; l < 15; ++l) {
                  auto sj1 = Angular::sixj_2(Fa.twoj(), Fb.twoj(), Fc.twoj(),
                                             Fd.twoj(), 2 * k, 2 * l);
                  auto sj2 = Coulomb::sixjTriads(Fa, Fb, Fc, Fd, k, l) ?
                                 Coulomb::sixj(Fa, Fb, Fc, Fd, k, l) :
                                 0.0;
                  if (std::abs(sj1 - sj2) > worst)
                    worst = std::abs(sj1 - sj2);
                }
              }
            }
          }
        }
      }
      pass &= qip::check_value(&obuff, "Magic 6J ", worst, 0.0, 1.0e-13);
    }
  }

  return pass;
}

} // namespace UnitTest

//****************************************************************************
//****************************************************************************

//****************************************************************************
inline std::vector<double> UnitTest::helper::yk_naive(const DiracSpinor &Fa,
                                                      const DiracSpinor &Fb,
                                                      int k) {
  const auto &gr = *Fa.rgrid;
  std::vector<double> yk(gr.r().size());
#pragma omp parallel for
  for (auto i = 0ul; i < yk.size(); ++i) {
    auto r = gr.r(i);

    auto rtkr = [](double x, double y, int kk) {
      return x < y ? std::pow(x / y, kk) / y : std::pow(y / x, kk) / x;
    };

    std::vector<double> f;
    f.reserve(gr.r().size());
    for (auto j = 0ul; j < yk.size(); ++j) {
      f.push_back(rtkr(r, gr.r(j), k) *
                  (Fa.f(j) * Fb.f(j) + Fa.g(j) * Fb.g(j)));
    }

    const auto p0 = std::max(Fa.min_pt(), Fb.min_pt());
    const auto pi = std::min(Fa.max_pt(), Fb.max_pt());
    yk[i] = NumCalc::integrate(gr.du(), p0, pi, f, gr.drdu());
  }

  return yk;
}

//****************************************************************************
inline double
UnitTest::helper::check_ykab_Tab(const std::vector<DiracSpinor> &a,
                                 const std::vector<DiracSpinor> &b,
                                 const Coulomb::YkTable &Yab) {
  //
  double worst = 0.0;
  for (const auto &Fa : a) {
    for (const auto &Fb : b) {
      const auto [kmin, kmax] = Coulomb::k_minmax(Fa, Fb);
      for (int k = kmin; k <= kmax; ++k) {
        // Only check if Angular factor is non-zero (since Ykab only calc'd in
        // this case)
        if (!Angular::Ck_kk_SR(k, Fa.k, Fb.k))
          continue;
        const auto y1 = Yab.get(k, Fa, Fb);
        const auto y2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto y3 = Coulomb::yk_ab(Fb, Fa, k); // symmetric
        if (y1 == nullptr) {
          std::cout << k << " " << Fa.symbol() << " " << Fb.symbol() << "\n";
        }
        const auto del = std::abs(qip::compare(*y1, y2).first) +
                         std::abs(qip::compare(y2, y3).first);
        if (del > worst)
          worst = del;
      }
    }
  }

  return worst;
}

//****************************************************************************
inline std::vector<double>
UnitTest::helper::check_ykab(const std::vector<DiracSpinor> &orbs,
                             int max_del_n) {

  // Compared Yk as calculated by fast Coulomb::yk_ab routine (used in the code)
  // to helper::yk_naive, a very slow, but simple version. In theory, should be
  // exactly the same

  std::vector<double> worst; // = 0.0;
  // nb: slow, so only check sub-set
  for (auto ia = 0ul; ia < orbs.size(); ++ia) {
    const auto &Fa = orbs[ia];
    for (auto ib = ia; ib < orbs.size(); ++ib) {
      const auto &Fb = orbs[ib];
      const auto [kmin, kmax] = Coulomb::k_minmax(Fa, Fb);
      for (int k = kmin; k <= kmax; ++k) {
        if (!Angular::Ck_kk_SR(k, Fa.k, Fb.k))
          continue;
        if (std::abs(Fa.n - Fb.n) > max_del_n)
          continue;
        if (std::size_t(k + 1) > worst.size())
          worst.resize(std::size_t(k + 1));
        const auto y2 = Coulomb::yk_ab(Fa, Fb, k);
        const auto y4 = helper::yk_naive(Fa, Fb, k); // slow

        const auto max = *std::max_element(cbegin(y4), cend(y4), qip::comp_abs);
        const auto del = std::abs(qip::compare(y2, y4).first / max);

        if (del > worst[std::size_t(k)])
          worst[std::size_t(k)] = del;
      }
    }
  }

  return worst;
}

//****************************************************************************
inline double
UnitTest::helper::check_Rkabcd(const std::vector<DiracSpinor> &orbs,
                               int max_del_n) {
  double eps_R = 0.0;
  const Coulomb::YkTable Yab(orbs);
#pragma omp parallel for
  for (auto ia = 0ul; ia < orbs.size(); ia++) {
    const auto &Fa = orbs[ia];
    for (auto ib = 0ul; ib < orbs.size(); ib += 2) {
      const auto &Fb = orbs[ib];
      for (auto ic = ia; ic < orbs.size(); ic++) {
        const auto &Fc = orbs[ic];
        if (std::abs(Fa.n - Fc.n) > max_del_n)
          continue;
        if (std::abs(Fb.n - Fc.n) > max_del_n)
          continue;
        for (auto id = ib; id < orbs.size(); id += 2) {
          const auto &Fd = orbs[id];
          if (std::abs(Fb.n - Fd.n) > max_del_n)
            continue;
          const auto [kmin, kmax] = Coulomb::k_minmax(Fa, Fc);
          for (int k = kmin; k <= kmax; ++k) {
            if (!Angular::Ck_kk_SR(k, Fa.k, Fc.k))
              continue;
            if (!Angular::Ck_kk_SR(k, Fb.k, Fd.k))
              continue;

            //---------
            const auto yac = Yab.get(k, Fa, Fc);
            const auto ybd = Yab.get(k, Fb, Fd);

            const auto r1a = Coulomb::Rk_abcd(Fa, Fb, Fc, Fd, k);
            const auto r1b = Coulomb::Rk_abcd(Fb, Fa, Fd, Fc, k);
            const auto r1c = Coulomb::Rk_abcd(Fc, Fd, Fa, Fb, k);
            const auto r2a = Coulomb::Rk_abcd(Fa, Fc, *ybd);
            const auto r2b = Coulomb::Rk_abcd(Fb, Fd, *yac);
            const auto r2c = Coulomb::Rk_abcd(Fc, Fa, *ybd);
            const auto r3 = Fa * Coulomb::Rkv_bcd(Fa.k, Fb, Fc, Fd, k);
            const auto r4 = Fa * Coulomb::Rkv_bcd(Fa.k, Fc, *ybd);
            const auto eps = std::max({r1a, r1b, r1c, r2a, r2b, r2c, r3, r4}) -
                             std::min({r1a, r1b, r1c, r2a, r2b, r2c, r3, r4});
#pragma omp critical(compare_epsR)
            if (eps > eps_R) {
              eps_R = eps;
            }
          }
        }
      }
    }
  }
  return eps_R;
}
