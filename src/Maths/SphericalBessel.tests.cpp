#include "SphericalBessel.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <iostream>
#include <vector>

TEST_CASE("Maths::SphericalBessel", "[jL][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "SphericalBessel\n";

  const auto r_list = qip::loglinear_range(1.0e-5, 10.0, 1.0, 100);

  for (int L = 0; L <= 15; ++L) {
    for (auto r : r_list) {
      REQUIRE(SphericalBessel::JL(L, r) ==
              Approx(SphericalBessel::exactGSL_JL(L, r)).margin(1.0e-15));

      REQUIRE(SphericalBessel::JL(L, r) ==
              Approx(SphericalBessel::exactGSL_JL(L, r)).epsilon(1.0e-3));
    }

    if (L < 4) {
      const auto jl = SphericalBessel::fillBesselVec(L, r_list);
      const auto jlkr = SphericalBessel::fillBesselVec_kr(L, 0.15, r_list);
      for (auto i = 0ul; i < r_list.size(); ++i) {
        const auto r = r_list.at(i);
        REQUIRE(jl.at(i) == Approx(SphericalBessel::JL(L, r)));
        REQUIRE(jlkr.at(i) == Approx(SphericalBessel::JL(L, 0.15 * r)));
      }
    }
  }
}

//==============================================================================
TEST_CASE("Maths::JL_table", "[JL_table][jL][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "SphericalBessel: JL_table\n";

  const std::vector<double> r = qip::logarithmic_range(1.0e-5, 100.0, 20);
  // dense q grid, for interp etc.
  const std::vector<double> q = qip::logarithmic_range(1.0e-2, 1000.0, 5000);

  const std::vector<std::size_t> iq_exact = {
    0ul,   1ul,   2ul,          15ul,         279ul,
    556ul, 999ul, q.size() - 3, q.size() - 2, q.size() - 1};

  const int max_L = 4;
  SphericalBessel::JL_table table(max_L, q, r);

  SECTION("at() returns correct values") {
    for (int L = 0; L <= max_L; ++L) {
      for (std::size_t iq = 0; iq < q.size(); iq += 150) {
        const auto &v = table.at(std::size_t(L), iq);
        REQUIRE(v.size() == r.size());
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::JL(L, q[iq] * r[ir])).epsilon(1e-12));
        }
      }
    }
  }

  SECTION("jL() returns first q_i >= q") {
    for (int L = 0; L <= max_L; ++L) {
      for (double qq :
           {0.99e-2, 1.0e-2, 0.1, 25.0, 164.7, 999.0, 1000.0, 1001.0}) {
        const auto &v = table.jL(L, qq);

        for (std::size_t ir = 0; ir < r.size(); ir++) {
          // This one may be quite rough, particularly for small jL values
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).margin(1e-2));
        }
      }

      // test "exact" cases
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto &v = table.jL(L, qq);
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      // test "exact" cases, with small error
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto eps = 1.0e-5;
        const auto &v1 = table.jL(L, qq * (1.0 + eps));
        const auto &v2 = table.jL(L, qq * (1.0 - eps));
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          // less exact: went to far!
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).margin(1e-2));
          // very exact
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::JL(
                           L, q[std::min(iq + 1, q.size() - 1)] * r[ir]))
                    .epsilon(1e-12));
          // very exact
          REQUIRE(v2[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      //test edge cases:
      const auto &v0 = table.jL(L, 0.1 * q.front());
      const auto &v1 = table.jL(L, 2.0 * q.back());
      for (std::size_t ir = 0; ir < r.size(); ir++) {
        REQUIRE(
          v0[ir] ==
          Approx(SphericalBessel::JL(L, q.front() * r[ir])).epsilon(1e-12));
        REQUIRE(
          v1[ir] ==
          Approx(SphericalBessel::JL(L, q.back() * r[ir])).epsilon(1e-12));
      }
    }
  }

  SECTION("jL_nearest() returns for closest q") {
    for (int L = 0; L <= max_L; ++L) {
      for (double qq :
           {0.99e-2, 1.0e-2, 0.1, 25.0, 164.7, 999.0, 1000.0, 1001.0}) {
        const auto &v = table.jL_nearest(L, qq);

        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).margin(1e-2));
        }
      }

      // test "exact" cases
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto &v = table.jL_nearest(L, qq);
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      // test "exact" cases, with small error
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto eps = 1.0e-5;
        const auto &v1 = table.jL_nearest(L, qq * (1.0 + eps));
        const auto &v2 = table.jL_nearest(L, qq * (1.0 - eps));
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-12));
          REQUIRE(v2[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      //test edge cases:
      const auto &v0 = table.jL_nearest(L, 0.1 * q.front());
      const auto &v1 = table.jL_nearest(L, 2.0 * q.back());
      for (std::size_t ir = 0; ir < r.size(); ir++) {
        REQUIRE(
          v0[ir] ==
          Approx(SphericalBessel::JL(L, q.front() * r[ir])).epsilon(1e-12));
        REQUIRE(
          v1[ir] ==
          Approx(SphericalBessel::JL(L, q.back() * r[ir])).epsilon(1e-12));
      }
    }
  }

  SECTION("jL_interp() interpolates correctly") {
    for (int L = 0; L <= max_L; ++L) {
      for (double qq : {1.0e-2, 0.1, 25.0, 164.7, 999.0, 1000.0}) {
        const auto &v = table.jL_interp(L, qq);
        REQUIRE(v.size() == r.size());
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          if (qq * r[ir] < 1.0) {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-4));
          } else if (qq * r[ir] < 10.0) {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-3));
          } else if (qq * r[ir] < 50.0) {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::JL(L, qq * r[ir])).margin(1e-3));
          } else {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::JL(L, qq * r[ir])).margin(1e-2));
          }
        }
      }

      // test "exact" cases
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto &v = table.jL_interp(L, qq);
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      // test "exact" cases, with small error
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto eps = 1.0e-8;
        const auto &v1 = table.jL_interp(L, qq * (1.0 + eps));
        const auto &v2 = table.jL_interp(L, qq * (1.0 - eps));
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          // interp still not expected to be exact
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-3));
          REQUIRE(v2[ir] ==
                  Approx(SphericalBessel::JL(L, qq * r[ir])).epsilon(1e-3));
        }
      }

      //test edge cases:
      const auto &v0 = table.jL_interp(L, 0.1 * q.front());
      const auto &v1 = table.jL_interp(L, 2.0 * q.back());
      for (std::size_t ir = 0; ir < r.size(); ir++) {
        REQUIRE(
          v0[ir] ==
          Approx(SphericalBessel::JL(L, q.front() * r[ir])).epsilon(1e-12));
        REQUIRE(
          v1[ir] ==
          Approx(SphericalBessel::JL(L, q.back() * r[ir])).epsilon(1e-12));
      }
    }
  }
}