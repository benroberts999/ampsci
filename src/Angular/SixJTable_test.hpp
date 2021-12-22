#pragma once
#include "Angular/SixJTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Check.hpp"
#include <string>

namespace UnitTest {

//******************************************************************************
//! Unit tests for angular functions/classes (threeJ symbols, lookup tables etc)
bool SixJTable(std::ostream &obuff) {
  bool pass = true;

  {

    Angular::SixJTable sjt;
    // int max_2k = 15;

    auto fmt = [](auto i) {
      return i % 2 == 0 ? std::to_string(i / 2) : std::to_string(i) + "/2";
    };

    for (auto max_k : {7, 10, 12}) {

      {
        IO::ChronoTimer t("Fill " + std::to_string(max_k));
        sjt.fill(max_k); // nb: each loop, "extends" the table
      }

      std::stringstream ss;
      auto max_del = 0.0;
      auto max_2k = 2 * max_k;
      // compare values from table to direct calculation, no short-circuit
      for (int a = 0; a <= max_2k; ++a) {
        for (int b = 0; b <= max_2k; ++b) {
          for (int c = 0; c <= max_2k; ++c) {
            for (int d = 0; d <= max_2k; ++d) {
              for (int e = 0; e <= max_2k; ++e) {
                for (int f = 0; f <= max_2k; ++f) {
                  const auto sj = gsl_sf_coupling_6j(a, b, c, d, e, f);
                  const auto sj2 = sjt.get_2(a, b, c, d, e, f);
                  const auto del = std::abs(sj - sj2);
                  if (del > max_del) {
                    max_del = del;
                    ss.str("");
                    ss << "{" << fmt(a) << "," << fmt(b) << "," << fmt(c)
                       << " | " << fmt(d) << "," << fmt(e) << "," << fmt(f)
                       << "}";
                  }
                }
              }
            }
          }
        }
      }
      pass &= qip::check_value(&obuff,
                               "6J " + std::to_string(max_k) + ":" + ss.str(),
                               max_del, 0.0, 1.0e-15);
    }
  }

  // Test DiracSpinor
  {
    const int l_max = 6;

    // Construct a dummy set of DiracOrbitals
    std::vector<DiracSpinor> basis;
    for (auto l = 0; l <= l_max; ++l) {
      // 2j = 2l+1, 2l-1
      if (l != 0) {
        const auto kappa = Angular::kappa_twojl(2 * l - 1, l);
        basis.emplace_back(0, kappa, nullptr);
      }
      const auto kappa = Angular::kappa_twojl(2 * l + 1, l);
      basis.emplace_back(0, kappa, nullptr);
    }

    // Fill 6J table as usually done:
    const auto max_k = DiracSpinor::max_tj(basis); // max_k = 2*max_j
    const Angular::SixJTable sjt(max_k);

    // loop through all kinds of 6J symbols
    std::stringstream ss;
    auto max_del = 0.0;
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          for (const auto &d : basis) {
            // add 5, ensure our table is not missing any non-zero
            const auto tjmax =
                std::max({a.twoj(), b.twoj(), c.twoj(), d.twoj()}) + 5;
            for (int k = 0; k <= tjmax; ++k) {
              for (int l = 0; l <= tjmax; ++l) {
                const auto sj1 = gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k,
                                                    c.twoj(), d.twoj(), 2 * l);
                const auto sj2 = sjt.get(a, b, k, c, d, l);

                // Compare, find worst offender:
                const auto del = std::abs(sj1 - sj2);
                if (del > max_del) {
                  max_del = del;
                  ss.str("");
                  ss << "{" << a.shortSymbol() << "," << b.shortSymbol() << ","
                     << k << " | " << c.shortSymbol() << "," << d.shortSymbol()
                     << "," << l << "}";
                }
              }
            }
          }
        }
      }
    }
    pass &=
        qip::check_value(&obuff, "6J " + std::to_string(max_k) + ":" + ss.str(),
                         max_del, 0.0, 1.0e-15);

    max_del = 0.0;
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          const auto tjmax = 2 * std::max({a.twoj(), b.twoj(), c.twoj()});
          for (int k = 0; k <= tjmax; ++k) {
            for (int l = 0; l <= tjmax; ++l) {
              for (int u = 0; u <= tjmax; ++u) {
                const auto sj1 = gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k,
                                                    2 * l, 2 * u, c.twoj());
                const auto sj2 = sjt.get(a, b, k, l, u, c);

                // Compare, find worst offender:
                const auto del = std::abs(sj1 - sj2);
                if (del > max_del) {
                  max_del = del;
                  ss.str("");
                  ss << "{" << a.shortSymbol() << "," << b.shortSymbol() << ","
                     << k << " | " << l << "," << u << "," << c.shortSymbol()
                     << "}";
                }
              }
            }
          }
        }
      }
    }
    pass &=
        qip::check_value(&obuff, "6J " + std::to_string(max_k) + ":" + ss.str(),
                         max_del, 0.0, 1.0e-15);

    {
      // performance test
      double t1 = 0.0, t2 = 0.0;
      {
        IO::ChronoTimer t("New");
        double sum1 = 0.0;
        for (const auto &a : basis) {
          for (const auto &b : basis) {
            for (const auto &c : basis) {
              for (const auto &d : basis) {
                const auto tjmax =
                    std::max({a.twoj(), b.twoj(), c.twoj(), d.twoj()});
                for (int k = 0; k <= tjmax; ++k) {
                  for (int l = 0; l <= tjmax; ++l) {
                    sum1 += sjt.get(a, b, k, c, d, l);
                  }
                }
              }
              const auto tjmax = std::max({a.twoj(), b.twoj(), c.twoj()});
              for (int k = 0; k <= tjmax; ++k) {
                for (int l = 0; l <= tjmax; ++l) {
                  for (int u = 0; u <= tjmax; ++u) {
                    sum1 += sjt.get(a, b, k, l, u, c);
                  }
                }
              }
            }
          }
        }
        t1 = t.lap_reading_ms();
      }
      {
        IO::ChronoTimer t("Old");
        double sum2 = 0.0;
        for (const auto &a : basis) {
          for (const auto &b : basis) {
            for (const auto &c : basis) {
              for (const auto &d : basis) {
                const auto tjmax =
                    std::max({a.twoj(), b.twoj(), c.twoj(), d.twoj()});
                for (int k = 0; k <= tjmax; ++k) {
                  for (int l = 0; l <= tjmax; ++l) {
                    sum2 += gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k,
                                               c.twoj(), d.twoj(), 2 * l);
                  }
                }
              }
              const auto tjmax = std::max({a.twoj(), b.twoj(), c.twoj()});
              for (int k = 0; k <= tjmax; ++k) {
                for (int l = 0; l <= tjmax; ++l) {
                  for (int u = 0; u <= tjmax; ++u) {
                    sum2 += gsl_sf_coupling_6j(a.twoj(), b.twoj(), 2 * k, 2 * l,
                                               2 * u, c.twoj());
                  }
                }
              }
            }
          }
        }
        t2 = t.lap_reading_ms();
      }

      pass &= qip::check(&obuff, "6Jtab: >2x speed", t2 > 1.5 * t1, true);
    }
  }

  return pass;
}

} // namespace UnitTest
