#pragma once
#include <memory>
#include <string>
#include <utility>
#include <vector>
class Grid;

//******************************************************************************
/*!
@brief Stores radial Dirac spinor: F_nk = (f, g)
@details
\f[
\psi_{n\kappa m} = \frac{1}{r}
\begin{pmatrix}
  f_{n\kappa}(r)\,\Omega_{\kappa m}\\
  g_{n\kappa}(r)\,\Omega_{-\kappa m}
\end{pmatrix},
\quad
F_{n\kappa} =
\begin{pmatrix}
  f_{n\kappa}(r)\\
  g_{n\kappa}(r)
\end{pmatrix}
\f]

\par  Construction.
Takes in constant n and k=kappa values + grid
  - A pointer to the Grid is stored (to avoid many copies of Grid)
  - This means Grid must survive [ensured to be true, see WF class]

\par Operator Overloads
  - Intuative operator overloads are provided (Fa, Fb are DiracSpinors):
  - v * Fa, where v is a double or vector works the obvious way
  - Fa * Fb = <Fa|Fb>
  - Fa == Fb returns true if {na,ka}=={nb,kb}
  - Fa > Fb : first compares n, and then kappa (via kappa_index)
  - Fa +/- Fb : Adds/subtracts the two spinors (and updates p0/pinf)
  - You can make copies: auto Fnew = Fa
  - And you can re-asign: Fb = Fa (provided Fa and Fb have same n and kappa!)
*/
class DiracSpinor {

public:
  // [[deprecated("Pass in a shared_ptr to avoid a copy")]] DiracSpinor(
  //     int in_n, int in_k, const Grid &rgrid);
  DiracSpinor(int in_n, int in_k, std::shared_ptr<const Grid> in_rgrid);

  //! Radial Grid; links F[i] to F(r)
  std::shared_ptr<const Grid> rgrid;
  //! Principal quantum number
  const int n;
  //! Dirac quantum number, kappa
  const int k;
  //! Single-particle energy, not including rest energy
  double en = 0.0;
  //! Upper (large) radial component
  std::vector<double> f;
  //! Lower (small) radial component
  std::vector<double> g;
  //! `practical zero': p0 is first non-zero point for f(r) [usually p0=0]
  std::size_t p0 = 0;
  //! `practical infinity': pinf is last non-zero point for f(r)
  std::size_t pinf;

  //! Number of iterations until energy convergence (for latest routine only)
  int its = -1;
  //! Fractional energy convergence: eps = |(en'-en)/en|
  double eps = -1.0;
  //! Occupation fraction. =1 for closed shells. =1/(2j+1) for valence
  double occ_frac = 0.0;

private:
  const int m_twoj;
  const int m_l;
  const int m_parity;
  const int m_k_index;

public:
  //! Orbital angular momentum Q number
  int l() const { return m_l; }
  //! j(j+1)
  double jjp1() const { return 0.25 * double(m_twoj * (m_twoj + 2)); }
  int twoj() const { return m_twoj; }
  //! 2j+1
  int twojp1() const { return m_twoj + 1; }
  //! (-1)^l, returns +/- 1
  int parity() const { return m_parity; }
  //! kappa index (see Angular)
  int k_index() const { return m_k_index; }

  //! Single-electron term symbol (e.g., 6s_1/2). Gnuplot=true => 6s_{1/2}
  std::string symbol(bool gnuplot = false) const;
  //! e.g., 6p_1/2 => 6p-, 6p_3/2 => 6p+
  std::string shortSymbol() const;

  //! norm = Sqrt[<a|a>]
  double norm() const;
  const DiracSpinor &scale(const double factor);
  const DiracSpinor &scale(const std::vector<double> &v);
  //! By default normalises to 1, but can normalise to other number.
  void normalise(double norm_to = 1.0);
  //! Returns [f[p0]/f_max , f[pinf]/f_max] - for tests
  std::pair<double, double> r0pinfratio() const;

  //! rho(r) = sum_m |Psi^2|(r) = (2j+1) * x_occ * |Psi^2|(r)
  std::vector<double> rho() const;

  //! Number of occupied electrons: (2j+1)*occ_frac
  int num_electrons() const;

  //! r0 = r[p0] (in atomic units)
  double r0() const;
  //! rinf = r[pinf]
  double rinf() const;

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

  DiracSpinor &operator*=(const std::vector<double> &v);
  friend DiracSpinor operator*(const std::vector<double> &v, DiracSpinor rhs);

  // comparitor overloads:
  friend bool operator==(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator!=(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator<(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator>(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator<=(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator>=(const DiracSpinor &lhs, const DiracSpinor &rhs);

  // Custom comparitors
  static bool comp_l(const DiracSpinor &lhs, const DiracSpinor &rhs) {
    return lhs.m_l < rhs.m_l;
  }
  static bool comp_j(const DiracSpinor &lhs, const DiracSpinor &rhs) {
    return lhs.m_twoj < rhs.m_twoj;
  }
  static bool comp_ki(const DiracSpinor &lhs, const DiracSpinor &rhs) {
    return lhs.m_k_index < rhs.m_k_index;
  }
  static bool comp_en(const DiracSpinor &lhs, const DiracSpinor &rhs) {
    return lhs.en < rhs.en;
  }

  //! Returns worst |<a|b>| (or |<a|b>-1| for a=b) {val, state_names}
  static std::pair<double, std::string>
  check_ortho(const std::vector<DiracSpinor> &a,
              const std::vector<DiracSpinor> &b);

  //! Returns formatted states string (e.g., '7sp5d') given list of orbs
  static std::string state_config(const std::vector<DiracSpinor> &orbs);

  DiracSpinor &operator=(const DiracSpinor &);
  DiracSpinor(const DiracSpinor &) = default;
  ~DiracSpinor() = default;
};
