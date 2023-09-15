#include "Coulomb/QkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"
#include <random>

//==============================================================================

//==============================================================================
TEST_CASE("Coulomb: Qk Table", "[Coulomb][QkTable][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: Qk Table\n";

  const auto radial_grid = std::make_shared<const Grid>(
      GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;
  const int lmax = 6;

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::vector<DiracSpinor> orbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      orbs.push_back(DiracSpinor::exactHlike(n_min, l, radial_grid, zeff));
    }
    orbs.push_back(DiracSpinor::exactHlike(n_min, -l - 1, radial_grid, zeff));
  }

  Coulomb::QkTable qk;
  Coulomb::YkTable yk(orbs);
  qk.fill(orbs, yk);

  // Test the Q,P,W formulas
  for (const auto &Fa : orbs) {
    for (const auto &Fb : orbs) {
      for (const auto &Fc : orbs) {
        for (const auto &Fd : orbs) {
          // go through _every_ k (will incllude zeros!)
          for (int k = 0; k <= 2 * lmax; ++k) {
            const auto q1 = qk.Q(k, Fa, Fb, Fc, Fd);
            const auto q2 = yk.Q(k, Fa, Fb, Fc, Fd);
            const auto eps_q = q2 == 0.0 ? q1 - q2 : (q1 - q2) / q2;
            REQUIRE(std::abs(eps_q) < 1.0e-14);

            const auto w1 = qk.W(k, Fa, Fb, Fc, Fd);
            const auto w2 = yk.W(k, Fa, Fb, Fc, Fd);
            const auto eps_w = std::abs(w2) < 1.0e-6 ? w1 - w2 : (w1 - w2) / w2;
            REQUIRE(std::abs(eps_w) < 1.0e-8);

            const auto p1 = qk.P(k, Fa, Fb, Fc, Fd);
            const auto p2 = w2 - q2;
            const auto eps_p = std::abs(p2) < 1.0e-6 ? p1 - p2 : (p1 - p2) / p2;
            REQUIRE(std::abs(eps_p) < 1.0e-8);
          }
        }
      }
    }
  }
}

//==============================================================================
TEST_CASE("Coulomb: Q,W,N k Table", "[Coulomb][QkTable][unit]") {

  // just 4 dummy DiracSpinors to test set/get for CoulombTable
  const DiracSpinor a{0, -1, nullptr};
  const DiracSpinor b{0, 1, nullptr};
  const DiracSpinor c{0, -2, nullptr};
  const DiracSpinor d{0, 2, nullptr};

  double x = 16.5;
  double y = 3.5;
  int k = 3;
  const auto rand_str = qip::random_string(6);
  std::string file_name{"deleteme_" + rand_str + ".qk"};

  // Qk
  {
    Coulomb::QkTable q;
    REQUIRE(q.emptyQ());
    REQUIRE(q.count() == 0);
    REQUIRE(!q.contains(k, a, b, c, d));
    REQUIRE(q.Q(k, a, b, c, d) == 0.0);
    q.add(k, a, b, c, d, x);
    REQUIRE(!q.emptyQ());
    REQUIRE(q.count() == 1);
    REQUIRE(q.contains(k, a, b, c, d));
    REQUIRE(q.Q(k, a, b, c, d) == x);
    REQUIRE(q.Q(k, a, d, c, b) == x);
    REQUIRE(q.Q(k, b, c, d, a) == x);
    REQUIRE(q.Q(k, b, a, d, c) == x);
    REQUIRE(q.Q(k, c, d, a, b) == x);
    REQUIRE(q.Q(k, c, b, a, d) == x);
    REQUIRE(q.Q(k, d, a, b, c) == x);
    REQUIRE(q.Q(k, d, c, b, a) == x);
    q.update(k, a, b, c, d, y);
    REQUIRE(q.Q(k, a, b, c, d) == y);

    q.write(file_name);
    Coulomb::QkTable q2;
    q2.read(file_name);
    REQUIRE(!q2.emptyQ());
    REQUIRE(q2.count() == q.count());
    REQUIRE(q2.Q(k, a, b, c, d) == q.Q(k, a, b, c, d));
  }

  // Wk
  {
    Coulomb::WkTable w;
    REQUIRE(w.emptyQ());
    REQUIRE(w.count() == 0);
    REQUIRE(!w.contains(k, a, b, c, d));
    REQUIRE(w.Q(k, a, b, c, d) == 0.0);
    w.add(k, a, b, c, d, x);
    REQUIRE(!w.emptyQ());
    REQUIRE(w.count() == 1);
    REQUIRE(w.contains(k, a, b, c, d));
    REQUIRE(w.Q(k, a, b, c, d) == x);
    REQUIRE(w.Q(k, b, a, d, c) == x);
    REQUIRE(w.Q(k, c, d, a, b) == x);
    REQUIRE(w.Q(k, d, c, b, a) == x);
    REQUIRE(w.Q(k, a, d, c, b) != x);
    REQUIRE(w.Q(k, b, c, d, a) != x);
    REQUIRE(w.Q(k, c, b, a, d) != x);
    REQUIRE(w.Q(k, d, a, b, c) != x);
    w.update(k, a, b, c, d, y);
    REQUIRE(w.Q(k, a, b, c, d) == y);
  }

  // Lk
  {
    Coulomb::LkTable l;
    REQUIRE(l.emptyQ());
    REQUIRE(l.count() == 0);
    REQUIRE(!l.contains(k, a, b, c, d));
    REQUIRE(l.Q(k, a, b, c, d) == 0.0);
    l.add(k, a, b, c, d, x);
    REQUIRE(!l.emptyQ());
    REQUIRE(l.count() == 1);
    REQUIRE(l.contains(k, a, b, c, d));
    REQUIRE(l.Q(k, a, b, c, d) == x);
    REQUIRE(l.Q(k, b, a, d, c) == x);
    REQUIRE(l.Q(k, a, d, c, b) != x);
    REQUIRE(l.Q(k, b, c, d, a) != x);
    REQUIRE(l.Q(k, c, d, a, b) != x);
    REQUIRE(l.Q(k, c, b, a, d) != x);
    REQUIRE(l.Q(k, d, a, b, c) != x);
    REQUIRE(l.Q(k, d, c, b, a) != x);
    l.update(k, a, b, c, d, y);
    REQUIRE(l.Q(k, a, b, c, d) == y);
  }

  // Nk
  {
    Coulomb::NkTable n;
    REQUIRE(n.emptyQ());
    REQUIRE(n.count() == 0);
    REQUIRE(!n.contains(k, a, b, c, d));
    REQUIRE(n.Q(k, a, b, c, d) == 0.0);
    n.add(k, a, b, c, d, x);
    REQUIRE(!n.emptyQ());
    REQUIRE(n.count() == 1);
    REQUIRE(n.contains(k, a, b, c, d));
    REQUIRE(n.Q(k, a, b, c, d) == x);
    REQUIRE(n.Q(k, a, d, c, b) != x);
    REQUIRE(n.Q(k, b, c, d, a) != x);
    REQUIRE(n.Q(k, b, a, d, c) != x);
    REQUIRE(n.Q(k, c, d, a, b) != x);
    REQUIRE(n.Q(k, c, b, a, d) != x);
    REQUIRE(n.Q(k, d, a, b, c) != x);
    REQUIRE(n.Q(k, d, c, b, a) != x);
    n.update(k, a, b, c, d, y);
    REQUIRE(n.Q(k, a, b, c, d) == y);
  }
}

//==============================================================================
//==============================================================================
TEST_CASE("Coulomb: Qk Table - with WF", "[Coulomb][QkTable][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: Qk Table - with WF\n";

  using namespace Coulomb;

  Wavefunction wf({1000, 1.0e-5, 50.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]");
  wf.formBasis({"12spdfg", 30, 7, 1.0e-3, 1.0e-3, 30.0, false});

  // Form the Coulomb lookup tables:

  // YkTable stores Hartee Y-functions Y_ab(r)
  // These save much time when calculating Q^k coeficients
  const Coulomb::YkTable yk(wf.basis());

  Coulomb::QkTable qk_t;
  // Coulomb::WkTable qk;
  // Coulomb::NkTable qk;
  qk_t.fill(wf.basis(), yk);

  const auto rand_str = qip::random_string(6);
  std::string fname{"deleteme_" + rand_str + ".qk"};

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
        const auto q0 = Coulomb::Qk_abcd(k, a, b, c, d);
        const auto p0 = Coulomb::Pk_abcd(k, a, b, c, d);
        const auto w0 = Coulomb::Wk_abcd(k, a, b, c, d);
        const auto r0 = q0 != 0.0 ? Coulomb::Rk_abcd(k, a, b, c, d) : 0.0;
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

  //----------------------------------------------------------------------------
  {
    Coulomb::QkTable qk_ab;
    qk_ab.fill_ab(wf.basis(), yk);
    for (auto &a : wf.basis()) {
      for (auto &b : wf.basis()) {
        {
          auto [kmin, kmax] = k_minmax_Q(a, b, a, b);
          for (auto k = kmin; k <= kmax; ++k) {
            auto q1 = yk.Q(k, a, b, a, b);
            auto q2 = qk_ab.Q(k, a, b, a, b);
            REQUIRE(q1 == Approx(q2));
          }
        }
        {
          auto [kmin, kmax] = k_minmax_Q(a, b, b, a);
          for (auto k = kmin; k <= kmax; ++k) {
            auto q1 = yk.Q(k, a, b, b, a);
            auto q2 = qk_ab.Q(k, a, b, b, a);
            if (q1 != Approx(q2)) {
              std::cout << k << " " << a << b << b << a << "\n";
              std::cout << yk.Q(k, a, b, b, a) << " " << qk_ab.Q(k, a, b, b, a)
                        << "\n";
              std::cout << yk.Q(k, a, a, b, b) << " " << qk_ab.Q(k, a, a, b, b)
                        << "\n";
            }
            REQUIRE(q1 == Approx(q2));
          }
        }
        {
          auto [kmin, kmax] = k_minmax_W(a, b, a, b);
          for (auto k = kmin; k <= kmax; ++k) {
            auto w1 = yk.W(k, a, b, a, b);
            auto w2 = qk_ab.W(k, a, b, a, b);
            REQUIRE(w1 == Approx(w2));
          }
        }
      }
    }
  }

  //----------------------------------------------------------------------------
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

//==============================================================================
TEST_CASE("Coulomb: Qk Table - performance",
          "[Coulomb][QkTable][performance][!mayfail]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: Qk Table - performance\n";

  // Use more realistic orbitals (have rmax, which makes 'Direct' calc faster)
  Wavefunction wf({1600, 1.0e-6, 120.0, 40.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", 0.0, "[Xe]");
  wf.formBasis({"10spdfghi", 30, 7, 1.0e-3, 1.0e-3, 30.0, false});
  const auto &basis = wf.basis();
  // wf.printCore();
  // wf.printBasis(basis);

  Coulomb::QkTable qk;
  Coulomb::YkTable yk(basis);
  qk.fill(basis, yk);
  const int num_runs = 3;

  // Compare the speed of using Qk lookup table vs. direct calculation
  double dir_time = 0.0;
  double tab_time = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  {
    IO::ChronoTimer t("Direct calc");
    for (int i = 0; i < num_runs; ++i) {
#pragma omp parallel for reduction(+ : sum1)
      for (auto ia = 0ul; ia < basis.size(); ++ia) {
        auto &a = basis[ia];
        for (const auto &b : basis) {
          for (const auto &c : basis) {
            for (const auto &d : basis) {
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
    }
    std::cout << "sum1=" << sum1 << "\n" << std::flush;
    dir_time = t.reading_ms();
  }
  std::cout << "\n";

  {
    IO::ChronoTimer t("Use table");
    for (int i = 0; i < num_runs; ++i) {
#pragma omp parallel for reduction(+ : sum2)
      for (auto ia = 0ul; ia < basis.size(); ++ia) {
        const auto &a = basis[ia];
        for (const auto &b : basis) {
          for (const auto &c : basis) {
            for (const auto &d : basis) {
              const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
              for (int k = kmin; k <= kmax; k += 2) {
                sum2 += qk.Q(k, a, b, c, d);
              }
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

  REQUIRE(std::abs(sum1 - sum2) < 1.0e-6);
  REQUIRE(tab_time < dir_time);
  INFO("QkTable timing. tab time: " << tab_time
                                    << ", direct time: " << dir_time);
}
