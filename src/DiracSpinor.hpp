#pragma once
#include "ATI_atomInfo.hpp"
#include "Grid.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <string>
#include <vector>

class DiracSpinor {

public: // Data
  DiracSpinor(int in_n, int in_k, const Grid &rgrid, bool imaginary_g = true)
      : n(in_n), k(in_k), pinf(rgrid.ngp - 1), p_rgrid(&rgrid),
        imaginary_g(imaginary_g) {
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

  // determines relative sign in radial integral
  // Make private?
  const bool imaginary_g; // true by default. If false, means upper comp is i

  int its = -1;
  double eps = -1;
  double occ_frac = -1;

public: // Methods
  int l() const { return ATI::l_k(k); }
  double j() const { return ATI::j_k(k); }
  int twoj() const { return ATI::twoj_k(k); }
  int parity() const { return ATI::parity_k(k); }
  int k_index() const { return ATI::indexFromKappa(k); }

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

  double operator*(const DiracSpinor &other) const {
    int pinf = 0; // XXX goes to ngp...ok?

    // Change the relative sign based in Complex f or g component
    // (includes complex conjugation of lhs)
    int ffs = ((!this->imaginary_g) && (other.imaginary_g)) ? -1 : 1;
    int ggs = ((this->imaginary_g) && (!other.imaginary_g)) ? -1 : 1;

    auto ff =
        NumCalc::integrate(this->f, other.f, this->p_rgrid->drdu, 1, 0, pinf);
    auto gg =
        NumCalc::integrate(this->g, other.g, this->p_rgrid->drdu, 1, 0, pinf);
    return (ffs * ff + ggs * gg) * this->p_rgrid->du;
  }
  DiracSpinor operator+(const DiracSpinor &other) const {
    DiracSpinor sum(other.n, other.k, *other.p_rgrid, other.imaginary_g);
    sum.pinf = other.pinf; //?
    sum.en = other.en;     //?
    if (other.imaginary_g != this->imaginary_g) {
      std::cerr << "\nFAIL92 in DiracSpinor. Trying to add (re,im)+(im,re)!\n";
      std::cout << other.imaginary_g << " " << this->imaginary_g << "\n";
      std::abort();
    }
    for (std::size_t i = 0; i < other.p_rgrid->ngp; i++) {
      sum.f[i] = this->f[i] + other.f[i];
      sum.g[i] = this->g[i] + other.g[i];
    }
    return sum;
  }
};
