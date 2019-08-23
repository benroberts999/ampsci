#pragma once
#include <string>
#include <utility>
#include <vector>
class Grid;

//******************************************************************************
class DiracSpinor {

public: // Data
  DiracSpinor(int in_n, int in_k, const Grid &rgrid, bool in_imag_g = true);

  // Would be better if some of this were private... getters/setters.....
  const Grid *const p_rgrid;
  const int n;
  const int k;
  double en = 0;
  std::vector<double> f;
  std::vector<double> g;
  std::size_t pinf;

  // determines relative sign in radial integral:
  // true by default. If false, means upper comp is i
  // XXX Kill this! ?
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
  // double j() const { return 0.5 * double(m_twoj); }
  double jjp1() const { return 0.25 * double(m_twoj * (m_twoj + 2)); }
  int twoj() const { return m_twoj; }
  int twojp1() const { return m_twoj + 1; }
  int parity() const { return m_parity; }
  int k_index() const { return m_k_index; }

  std::string symbol(bool gnuplot = false) const;
  std::string shortSymbol() const;

  double norm() const;
  void scale(const double factor);
  void normalise(double norm_to = 1.0);
  std::pair<double, double> r0pinfratio() const;

public:
  // Operator overloads
  friend double operator*(const DiracSpinor &lhs, const DiracSpinor &rhs);

  DiracSpinor &operator+=(const DiracSpinor &rhs);
  DiracSpinor &operator-=(const DiracSpinor &rhs);
  friend DiracSpinor operator+(DiracSpinor lhs, const DiracSpinor &rhs);
  friend DiracSpinor operator-(DiracSpinor lhs, const DiracSpinor &rhs);

  DiracSpinor &operator*=(const double x);
  friend DiracSpinor operator*(DiracSpinor lhs, const double x);
  friend DiracSpinor operator*(const double x, DiracSpinor rhs);
  friend DiracSpinor operator*(const std::vector<double> &v, DiracSpinor rhs);

  DiracSpinor &operator=(const DiracSpinor &other);

  // comparitor overloads:
  friend bool operator==(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator!=(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator<(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator>(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator<=(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator>=(const DiracSpinor &lhs, const DiracSpinor &rhs);
};
