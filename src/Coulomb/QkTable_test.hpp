#pragma once
#include "Angular/SixJTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include <random>

namespace UnitTest {

//******************************************************************************
//******************************************************************************
//! Unit tests for Coulomb integrals (y^k_ab, R^k_abcd, lookup tables etc).
//! Also: tests quadrature integation method
bool QkTable(std::ostream &obuff) {
  bool pass = true;

  using namespace Coulomb;

  Wavefunction wf({1000, 1.0e-5, 50.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]");
  wf.formBasis({"10spdfg", 30, 7, 1.0e-5, 1.0e-6, 30.0, false});

  // Form the Coulomb lookup tables:

  // YkTable stores Hartee Y-functions Y_ab(r)
  // These save much time when calculating Q^k coeficients
  const Coulomb::YkTable yk(wf.basis);

  Coulomb::QkTable qk_t;
  // Coulomb::WkTable qk;
  // Coulomb::NkTable qk;
  qk_t.fill(wf.basis, yk);

  const std::string fname = "tmp_delete_me.qk";

  {
    IO::ChronoTimer t("Write to disk");
    qk_t.write(fname);
  }

  // Read in to qk (test of read/write)
  Coulomb::QkTable qk;
  {
    IO::ChronoTimer t("Read from disk");
    auto ok = qk.read(fname);
    std::cout << (ok ? "yes" : "no") << "\n";
  }
  std::cout << "\n";

  // Compare the speed of using Qk lookup table vs. direct calculation
  double dir_time = 0.0;
  double tab_time = 0.0;
  {
    IO::ChronoTimer t("Direct calc");
    double sum1 = 0.0;
#pragma omp parallel for reduction(+ : sum1)
    for (auto ia = 0ul; ia < wf.basis.size(); ++ia) {
      auto &a = wf.basis[ia];
      for (const auto &b : wf.basis) {
        for (const auto &c : wf.basis) {
          for (const auto &d : wf.basis) {
            const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
            for (int k = kmin; k <= kmax; k += 2) {
              const auto yk_bd = yk.get(k, b, d);
              if (yk_bd == nullptr)
                continue;
              sum1 += yk.Q(k, a, b, c, d);
            }
          }
        }
      }
    }
    std::cout << "sum1=" << sum1 << "\n" << std::flush;
    dir_time = t.reading_ms();
  }
  std::cout << "\n";

  {
    IO::ChronoTimer t("Use table");
    double sum2 = 0.0;
#pragma omp parallel for reduction(+ : sum2)
    for (auto ia = 0ul; ia < wf.basis.size(); ++ia) {
      const auto &a = wf.basis[ia];
      for (const auto &b : wf.basis) {
        for (const auto &c : wf.basis) {
          for (const auto &d : wf.basis) {
            const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
            for (int k = kmin; k <= kmax; k += 2) {
              sum2 += qk.Q(k, a, b, c, d);
            }
          }
        }
      }
    }
    std::cout << "sum2=" << sum2 << " (should = sum1)\n" << std::flush;
    tab_time = t.reading_ms();
  }

  std::cout << dir_time << "/" << tab_time
            << ": speed-up = " << dir_time / tab_time << "x\n";

  pass &= qip::check(&obuff, "QkTable: timing", tab_time < dir_time, true);

  {

    const auto max_2k = 2 * DiracSpinor::max_tj(wf.basis);
    Angular::SixJTable sjt{max_2k};

    // Test number of random instances of Q,P,R,W against direct way:
    // (most are zero)

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> rindex(0, wf.basis.size() - 1);

    // Test table versions of Q,P,W,R vs. 'from scratch' ones.
    // (only for subset, since otherwise v. slow)
    double max_devR = 0.0, max_devQ = 0.0;
    double max_devP = 0.0, max_devW = 0.0;
    const int num_to_test = 20000;
    int non_zero_count = 0;
    int total_count = 0; // more than num_to_test, because \sum_k
    for (int tries = 0; tries < num_to_test; ++tries) {
      const auto &a = wf.basis[rindex(gen)];
      const auto &b = wf.basis[rindex(gen)];
      const auto &c = wf.basis[rindex(gen)];
      const auto &d = wf.basis[rindex(gen)];

      const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
      // nb: test every k (test of k_minmax_Q..)
      for (int k = 0; k < kmax + 3; ++k) {
        const auto q1 = qk.Q(k, a, b, c, d);
        const auto p1 = qk.P(k, a, b, c, d);
        const auto p2 = qk.P(k, a, b, c, d, &sjt);
        const auto w1 = qk.W(k, a, b, c, d);
        const auto r1 = qk.R(k, a, b, c, d);
        const auto q0 = Coulomb::Qk_abcd(a, b, c, d, k);
        const auto p0 = Coulomb::Pk_abcd(a, b, c, d, k);
        const auto w0 = Coulomb::Wk_abcd(a, b, c, d, k);
        const auto r0 = q0 != 0.0 ? Coulomb::Rk_abcd(a, b, c, d, k) : 0.0;
        // nb: QkTable only stores Rk if Qk is non-zero, but Coulomb::Rk_abcd
        // will calculate it anyway.
        const auto devQ = std::abs(q1 - q0);
        const auto devP = std::max(std::abs(p1 - p0), std::abs(p2 - p0));
        const auto devW = std::abs(w1 - w0);
        const auto devR = std::abs(r1 - r0);
        if (devR > max_devR)
          max_devR = devR;
        if (devQ > max_devQ)
          max_devQ = devQ;
        if (devP > max_devP)
          max_devP = devP;
        if (devW > max_devW)
          max_devW = devW;

        ++total_count;
        if (q1 != 0.0)
          ++non_zero_count;
      }
    }
    std::cout << "tested: " << non_zero_count << "/" << total_count
              << " non-zero\n";
    pass &= qip::check_value(&obuff, "QkTable: Q", max_devQ, 0.0, 1.0e-13);
    pass &= qip::check_value(&obuff, "QkTable: R", max_devR, 0.0, 1.0e-13);
    pass &= qip::check_value(&obuff, "QkTable: P", max_devP, 0.0, 1.0e-13);
    pass &= qip::check_value(&obuff, "QkTable: W", max_devW, 0.0, 1.0e-13);
  }

  return pass;
}

} // namespace UnitTest
