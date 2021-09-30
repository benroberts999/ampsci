#pragma once
#include "Angular/SixJTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "qip/Check.hpp"
#include <string>

namespace UnitTest {

//******************************************************************************
//! Unit tests for angular functions/classes (threeJ symbols, lookup tables etc)
bool SixJTable(std::ostream &obuff) {
  bool pass = true;

  Angular::SixJTable sjt;
  // int max_2k = 15;

  auto fmt = [](auto i) {
    return i % 2 == 0 ? std::to_string(i / 2) : std::to_string(i) + "/2";
  };

  for (auto max_2k : {15, 20, 25}) {

    {
      IO::ChronoTimer t("Fill " + std::to_string(max_2k));
      sjt.fill(max_2k); // nb: each loop, "extends" the table
    }

    std::stringstream ss;
    auto max_del = 0.0;
    // compare values from table to direct calculation, no short-circuit
    for (int a = 0; a <= max_2k; ++a) {
      for (int b = 0; b <= max_2k; ++b) {
        for (int c = 0; c <= max_2k; ++c) {
          for (int d = 0; d <= max_2k; ++d) {
            for (int e = 0; e <= max_2k; ++e) {
              for (int f = 0; f <= max_2k; ++f) {
                const auto sj = gsl_sf_coupling_6j(a, b, c, d, e, f);
                const auto sj2 = sjt.get(a, b, c, d, e, f);
                const auto del = std::abs(sj - sj2);
                if (del > max_del) {
                  max_del = del;
                  ss.str("");
                  ss << "{" << fmt(a) << "," << fmt(b) << "," << fmt(c) << " | "
                     << fmt(d) << "," << fmt(e) << "," << fmt(f) << "}";
                }
              }
            }
          }
        }
      }
    }
    pass &= qip::check_value(&obuff,
                             "6J " + std::to_string(max_2k) + ":" + ss.str(),
                             max_del, 0.0, 1.0e-15);
  }

  // the following is performance test.
  // This shows ~2x speed up. However, in real-life application (ladder
  // diagrams), I observed a 20x speedup compared to manually calculating 6j!?

  auto max_2k = sjt.max_2jk();

  double t1 = 0.0, t2 = 0.0;

  {
    IO::ChronoTimer t("Direct calc");
    double sum = 0.0;
    for (int k1 = 0; k1 <= max_2k; k1 += 2) {
      for (int k2 = 0; k2 <= max_2k; k2 += 2) {
        for (int k3 = 0; k3 <= max_2k; k3 += 2) {
          for (int j1 = 1; j1 <= max_2k; j1 += 2) {
            for (int j2 = 1; j2 <= max_2k; j2 += 2) {
              for (int j3 = 1; j3 <= max_2k; j3 += 2) {
                auto sj = Angular::sixj_2(k1, k2, k3, j1, j2, j3);
                sum += sj;
              }
            }
          }
        }
      }
    }
    for (int j1 = 1; j1 <= max_2k; j1 += 2) {
      for (int j2 = 1; j2 <= max_2k; j2 += 2) {
        for (int k = 0; k <= max_2k; k += 2) {
          for (int j3 = 1; j3 <= max_2k; j3 += 2) {
            for (int j4 = 1; j4 <= max_2k; j4 += 2) {
              for (int l = 0; l <= max_2k; l += 2) {
                auto sj = Angular::sixj_2(j1, j2, k, j3, j4, l);
                sum += sj;
              }
            }
          }
        }
      }
    }
    std::cout << sum << "\n";
    t1 = t.lap_reading_ms();
  }

  {
    IO::ChronoTimer t("Table");
    double sum = 0.0;
    for (int k1 = 0; k1 <= max_2k; k1 += 2) {
      for (int k2 = 0; k2 <= max_2k; k2 += 2) {
        for (int k3 = 0; k3 <= max_2k; k3 += 2) {
          for (int j1 = 1; j1 <= max_2k; j1 += 2) {
            for (int j2 = 1; j2 <= max_2k; j2 += 2) {
              for (int j3 = 1; j3 <= max_2k; j3 += 2) {
                auto sj = sjt(k1, k2, k3, j1, j2, j3);
                sum += sj;
              }
            }
          }
        }
      }
    }
    for (int j1 = 1; j1 <= max_2k; j1 += 2) {
      for (int j2 = 1; j2 <= max_2k; j2 += 2) {
        for (int k = 0; k <= max_2k; k += 2) {
          for (int j3 = 1; j3 <= max_2k; j3 += 2) {
            for (int j4 = 1; j4 <= max_2k; j4 += 2) {
              for (int l = 0; l <= max_2k; l += 2) {
                auto sj = sjt(j1, j2, k, j3, j4, l);
                sum += sj;
              }
            }
          }
        }
      }
    }
    std::cout << sum << "\n";
    t2 = t.lap_reading_ms();
  }

  pass &= qip::check(&obuff, "6Jtab: >2x speed", t1 > 2 * t2, true);

  return pass;
}

} // namespace UnitTest
