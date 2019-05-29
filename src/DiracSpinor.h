#pragma once
#include "ATI_atomInfo.h"
#include "Grid.h"
#include <string>
#include <vector>

class DiracSpinor {

public: // Data
  DiracSpinor(int in_n, int in_k, const Grid &rgrid)
      : n(in_n), k(in_k), pinf(rgrid.ngp - 1), p_rgrid(&rgrid) {
    f.resize(rgrid.ngp, 0);
    g.resize(rgrid.ngp, 0);
  }

  std::vector<double> f;
  std::vector<double> g;
  double en = 0;

  const int n;
  const int k;

  std::size_t pinf;
  const Grid *const p_rgrid;

  int its = -1;
  double eps = -1;
  double occ_frac = -1;

public: // Methods
  int l() const { return ATI::l_k(k); }
  double j() const { return ATI::j_k(k); }
  int twoj() const { return ATI::twoj_k(k); }
  int parity() const { return ATI::parity_k(k); }

  std::string symbol(bool gnuplot = false) const {
    // Readable symbol (s_1/2, p_{3/2} etc.).
    // gnuplot-firndly '{}' braces optional.
    std::string ostring1 = std::to_string(n) + ATI::l_symbol(l());
    std::string ostring2 = gnuplot ? "_{" + std::to_string(twoj()) + "/2}"
                                   : "_" + std::to_string(twoj()) + "/2";
    return ostring1 + ostring2;
  }

public: // comparitor overloads
  // Define the custom comparitors (based on n and k)
  // Is there a way to do this all together?
  bool operator==(const DiracSpinor &other) const {
    return n == other.n && k == other.k;
  }
  bool operator!=(const DiracSpinor &other) const { return !(*this == other); }

  // Note: these are a little slow
  bool operator<(const DiracSpinor &other) const {
    if (n == other.n)
      return ATI::indexFromKappa(k) < ATI::indexFromKappa(other.k);
    return n < other.n;
  }
  bool operator>=(const DiracSpinor &other) const { return !(*this < other); }
  bool operator>(const DiracSpinor &other) const {
    if (n == other.n)
      return ATI::indexFromKappa(k) > ATI::indexFromKappa(other.k);
    return n > other.n;
  }
  bool operator<=(const DiracSpinor &other) const { return !(*this > other); }
};
