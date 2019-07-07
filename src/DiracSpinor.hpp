#pragma once
#include "AtomInfo.hpp"
#include "Grid.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <string>
#include <vector>

//******************************************************************************
class DiracSpinor {

public: // Data
  DiracSpinor(int in_n, int in_k, const Grid &rgrid, bool imaginary_g = true)
      : p_rgrid(&rgrid),                      //
        n(in_n), k(in_k), en(0.0),            //
        f(std::vector<double>(rgrid.ngp, 0)), //
        g(f),                                 //
        pinf(rgrid.ngp - 1),                  //
        imaginary_g(imaginary_g),             //
        its(-1), eps(-1), occ_frac(0),        //
        m_twoj(AtomInfo::twoj_k(in_k)),       //
        m_l(AtomInfo::l_k(in_k)),             //
        m_parity(AtomInfo::parity_k(in_k)),   //
        m_k_index(AtomInfo::indexFromKappa(in_k)) {}
  const Grid *const p_rgrid;

  const int n;
  const int k;

  double en = 0;
  std::vector<double> f;
  std::vector<double> g;
  std::size_t pinf;

  // determines relative sign in radial integral:
  // true by default. If false, means upper comp is i
  const bool imaginary_g;

  int its;
  double eps;
  double occ_frac;

private:
  const int m_twoj;
  const int m_l;
  const int m_parity;
  const int m_k_index;

public: // Methods
  int l() const { return m_l; }
  double j() const { return double(m_twoj) / 2; }
  int twoj() const { return m_twoj; }
  int parity() const { return m_parity; }
  int k_index() const { return m_k_index; }

  std::string symbol(bool gnuplot = false) const {
    // Readable symbol (s_1/2, p_{3/2} etc.).
    // gnuplot-firndly '{}' braces optional.
    std::string ostring1 = std::to_string(n) + AtomInfo::l_symbol(l());
    std::string ostring2 = gnuplot ? "_{" + std::to_string(twoj()) + "/2}"
                                   : "_" + std::to_string(twoj()) + "/2";
    return ostring1 + ostring2;
  }

  void normalise() {
    double norm = 1. / sqrt((*this) * (*this));
    for (auto &fa_r : f)
      fa_r *= norm;
    for (auto &ga_r : g)
      ga_r *= norm;
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
      return AtomInfo::indexFromKappa(k) < AtomInfo::indexFromKappa(other.k);
    return n < other.n;
  }
  bool operator>=(const DiracSpinor &other) const { return !(*this < other); }
  bool operator>(const DiracSpinor &other) const {
    if (n == other.n)
      return AtomInfo::indexFromKappa(k) > AtomInfo::indexFromKappa(other.k);
    return n > other.n;
  }
  bool operator<=(const DiracSpinor &other) const { return !(*this > other); }

public: // Operator overloads
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
