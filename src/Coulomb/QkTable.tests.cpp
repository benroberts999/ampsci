#include "Coulomb/QkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <random>

//==============================================================================
//==============================================================================
//! Unit tests for Coulomb integrals (y^k_ab, R^k_abcd, lookup tables etc).
//! Also: tests quadrature integation method

TEST_CASE("Coulomb: Qk Table", "[Coulomb][QkTable]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: Qk Table, [Coulomb][QkTable]\n";

  using namespace Coulomb;

  Wavefunction wf({1000, 1.0e-5, 50.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]");
  wf.formBasis({"10spdfg", 30, 7, 1.0e-5, 1.0e-6, 30.0, false});

  // Form the Coulomb lookup tables:

  // YkTable stores Hartee Y-functions Y_ab(r)
  // These save much time when calculating Q^k coeficients
  const Coulomb::YkTable yk(wf.basis());

  Coulomb::QkTable qk_t;
  // Coulomb::WkTable qk;
  // Coulomb::NkTable qk;
  qk_t.fill(wf.basis(), yk);

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
    for (auto ia = 0ul; ia < wf.basis().size(); ++ia) {
      auto &a = wf.basis()[ia];
      for (const auto &b : wf.basis()) {
        for (const auto &c : wf.basis()) {
          for (const auto &d : wf.basis()) {
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
    for (auto ia = 0ul; ia < wf.basis().size(); ++ia) {
      const auto &a = wf.basis()[ia];
      for (const auto &b : wf.basis()) {
        for (const auto &c : wf.basis()) {
          for (const auto &d : wf.basis()) {
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

  REQUIRE(tab_time < dir_time);
  INFO("QkTable timing. tab time: " << tab_time
                                    << ", direct time: " << dir_time);

  {
    const auto max_2k = 2 * DiracSpinor::max_tj(wf.basis());
    Angular::SixJTable sjt{max_2k};

    // Test number of random instances of Q,P,R,W against direct way:
    // (most are zero)

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> rindex(0, wf.basis().size() - 1);

    // Test table versions of Q,P,W,R vs. 'from scratch' ones.
    // (only for subset, since otherwise v. slow)
    double max_devR = 0.0, max_devQ = 0.0;
    double max_devP = 0.0, max_devW = 0.0;
    const int num_to_test = 50000;
    int non_zero_count = 0;
    int total_count = 0; // more than num_to_test, because \sum_k
    for (int tries = 0; tries < num_to_test; ++tries) {
      const auto &a = wf.basis()[rindex(gen)];
      const auto &b = wf.basis()[rindex(gen)];
      const auto &c = wf.basis()[rindex(gen)];
      const auto &d = wf.basis()[rindex(gen)];

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
    REQUIRE(std::abs(max_devQ) < 1.0e-13);
    REQUIRE(std::abs(max_devR) < 1.0e-13);
    REQUIRE(std::abs(max_devP) < 1.0e-13);
    REQUIRE(std::abs(max_devW) < 1.0e-13);
  }

  //============================================================================
  {
    // check Normal Ordering:
    bool Qk_NormalOrder_ok = true;
    for (auto &a : wf.basis()) {
      for (auto &b : wf.basis()) {
        for (auto &c : wf.basis()) {
          for (auto &d : wf.basis()) {
            // Qk: abcd = cbad = adcb = cdab = badc = bcda = dabc = dcba
            const auto i_abcd = qk.NormalOrder(a, b, c, d);
            const auto i_adcb = qk.NormalOrder(a, d, c, b);
            const auto i_badc = qk.NormalOrder(b, a, d, c);
            const auto i_bcda = qk.NormalOrder(b, c, d, a);
            const auto i_cbad = qk.NormalOrder(c, b, a, d);
            const auto i_cdab = qk.NormalOrder(c, d, a, b);
            const auto i_dabc = qk.NormalOrder(d, a, b, c);
            const auto i_dcba = qk.NormalOrder(d, c, b, a);
            // If all are the same, max = min
            const auto i_max = std::max({i_abcd, i_adcb, i_badc, i_bcda, i_cbad,
                                         i_cdab, i_dabc, i_dcba});
            const auto i_min = std::min({i_abcd, i_adcb, i_badc, i_bcda, i_cbad,
                                         i_cdab, i_dabc, i_dcba});
            if (i_max != i_min)
              Qk_NormalOrder_ok = false;
            if (i_max != i_min) {
              std::cout << "Qk NormalOrder failed for:\n";
              std::cout << a << b << c << d << ":" << i_abcd << "\n";
              std::cout << a << d << c << b << ":" << i_adcb << "\n";
              std::cout << b << a << d << c << ":" << i_badc << "\n";
              std::cout << b << c << d << a << ":" << i_bcda << "\n";
              std::cout << c << b << a << d << ":" << i_cbad << "\n";
              std::cout << c << d << a << b << ":" << i_cdab << "\n";
              std::cout << d << a << b << c << ":" << i_dabc << "\n";
              std::cout << d << c << b << a << ":" << i_dcba << "\n";
            }
          }
        }
      }
    }
    REQUIRE(Qk_NormalOrder_ok);

    Coulomb::WkTable wk;
    bool Wk_NormalOrder_ok = true;
    for (auto &a : wf.basis()) {
      for (auto &b : wf.basis()) {
        for (auto &c : wf.basis()) {
          for (auto &d : wf.basis()) {
            // Qk: abcd = badc = cdab = dcba
            const auto i_abcd = wk.NormalOrder(a, b, c, d);
            const auto i_badc = wk.NormalOrder(b, a, d, c);
            const auto i_cdab = wk.NormalOrder(c, d, a, b);
            const auto i_dcba = wk.NormalOrder(d, c, b, a);
            // If all are the same, max = min
            const auto i_max = std::max({i_abcd, i_badc, i_cdab, i_dcba});
            const auto i_min = std::min({i_abcd, i_badc, i_cdab, i_dcba});
            if (i_max != i_min)
              Wk_NormalOrder_ok = false;
            if (i_max != i_min) {
              std::cout << "Wk NormalOrder failed for:\n";
              std::cout << a << b << c << d << ":" << i_abcd << "\n";
              std::cout << b << a << d << c << ":" << i_badc << "\n";
              std::cout << c << d << a << b << ":" << i_cdab << "\n";
              std::cout << d << c << b << a << ":" << i_dcba << "\n";
            }
          }
        }
      }
    }
    REQUIRE(Wk_NormalOrder_ok);

    Coulomb::LkTable lk;
    bool Lk_NormalOrder_ok = true;
    for (auto &a : wf.basis()) {
      for (auto &b : wf.basis()) {
        for (auto &c : wf.basis()) {
          for (auto &d : wf.basis()) {
            // Qk: abcd = badc = cdab = dcba
            const auto i_abcd = wk.NormalOrder(a, b, c, d);
            const auto i_badc = wk.NormalOrder(b, a, d, c);
            if (i_abcd != i_badc)
              Lk_NormalOrder_ok = false;
            if (i_abcd != i_badc) {
              std::cout << "Lk NormalOrder failed for:\n";
              std::cout << a << b << c << d << ":" << i_abcd << "\n";
              std::cout << b << a << d << c << ":" << i_badc << "\n";
            }
          }
        }
      }
    }
    REQUIRE(Lk_NormalOrder_ok);
  }
}
