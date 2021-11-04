#pragma once
#include <memory>
#include <string>
#include <utility>
#include <vector>
class Grid;

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
  - A shared pointer to the Grid is stored (to avoid many copies of Grid)

\par Operator Overloads
  - Intuative operator overloads are provided (Fa, Fb are DiracSpinors):
  - v * Fa, where v is a double or vector works the obvious way
  - Fa * Fb = <Fa|Fb>
  - Fa == Fb returns true if {na,ka}=={nb,kb}
  - Fa > Fb : first compares n, and then kappa (via kappa_index)
  - Fa +/- Fb : Adds/subtracts the two spinors (and updates p0/pinf)
  - You can make copies: auto Fnew = Fa
  - And you can re-asign: Fb = Fa (provided Fa and Fb have same n and kappa!)
  - n and kappa are constant, cannot be changed. Avoids angular errors.
  - all 'set_' functions return mutable references to variables
*/
class DiracSpinor {

public:
  DiracSpinor(int in_n, int in_k, std::shared_ptr<const Grid> in_rgrid);

  //! Radial Grid; links F[i] to F(r)
  const std::shared_ptr<const Grid> rgrid;
  //! Principal quantum number
  const int n;
  //! Dirac quantum number, kappa
  const int k;

private:
  // Single-particle energy, not including rest energy
  double m_en = 0.0;
  // Upper (large) radial component
  std::vector<double> m_f;
  // Lower (small) radial component
  std::vector<double> m_g;
  // `practical zero': p0 is first non-zero point for f(r) [usually p0=0]
  std::size_t m_p0 = 0;
  // `practical infinity': pinf is last non-zero point for f(r)
  std::size_t m_pinf;
  // Number of iterations until energy convergence (for latest routine only)
  int m_its = -1;
  // Fractional energy convergence: eps = |(en'-en)/en|
  double m_eps = -1.0;
  // Occupation fraction. =1 for closed shells. =1/(2j+1) for valence
  double m_occ_frac = 0.0;

  // 2j, l, pi, kappa_index (for convenience)
  const int m_twoj;
  const int m_l;
  const int m_parity;
  const int m_k_index;
  using Index = uint16_t;
  const Index m_nk_index;

public:
  //! Single-particle energy, not including rest energy
  auto en() const { return m_en; }
  auto &set_en() { return m_en; }

  //! Upper (large) radial component, f(r)
  const auto &f() const { return m_f; }
  auto &set_f() { return m_f; }
  const auto &f(std::size_t i) const { return m_f.at(i); }
  auto &set_f(std::size_t i) { return m_f.at(i); }

  //! Lower (small) radial component, g(r)
  const auto &g() const { return m_g; }
  auto &set_g() { return m_g; }
  const auto &g(std::size_t i) const { return m_g.at(i); }
  auto &set_g(std::size_t i) { return m_g.at(i); }

  //! First non-zero point (index for f[i])
  auto min_pt() const { return m_p0; }
  auto &set_min_pt() { return m_p0; }

  //! Last non-zero point (index for f[i])
  auto max_pt() const { return m_pinf; }
  auto &set_max_pt() { return m_pinf; }

  //! r0 = r[min_pt] (in atomic units)
  double r0() const;
  //! rinf = r[max_pt]
  double rinf() const;

  //! Number of iterations until energy convergence (for latest routine only)
  auto its() const { return m_its; }
  auto &set_its() { return m_its; }

  //! Occupation fraction. =1 for closed shells. =1/(2j+1) for valence
  auto occ_frac() const { return m_occ_frac; }
  auto &set_occ_frac() { return m_occ_frac; }

  //! Fractional energy convergence: eps = |(en'-en)/en|
  auto eps() const { return m_eps; }
  auto &set_eps() { return m_eps; }

  //! Orbital angular momentum Q number
  int l() const { return m_l; }
  //! j(j+1)
  double jjp1() const { return 0.25 * double(m_twoj * (m_twoj + 2)); }
  int twoj() const { return m_twoj; }
  //! 2j+1
  int twojp1() const { return m_twoj + 1; }
  //! (-1)^l, returns +/- 1
  int parity() const { return m_parity; }
  //! kappa index (see AtomData)
  int k_index() const { return m_k_index; }
  //! (n,kappa) index (see AtomData)
  Index nk_index() const { return m_nk_index; }

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

  //! Forces f(r) and g(r) to be zero outside of [p0,pinf)
  void zero_boundaries();

  //! Returns [f[p0]/f_max , f[pinf]/f_max] - for tests
  std::pair<double, double> r0pinfratio() const;

  //! rho(r) = sum_m |Psi^2|(r) = (2j+1) * x_occ * |F^2|(r)
  std::vector<double> rho() const;

  //! Number of occupied electrons: (2j+1)*occ_frac
  int num_electrons() const;

public:
  // Operator overloads

  //! Returns radial integral (Fa,Fb) = Int(fa*fb + ga*gb)
  friend double operator*(const DiracSpinor &Fa, const DiracSpinor &Fb);

  DiracSpinor &operator+=(const DiracSpinor &rhs);
  DiracSpinor &operator-=(const DiracSpinor &rhs);
  friend DiracSpinor operator+(DiracSpinor lhs, const DiracSpinor &rhs);
  friend DiracSpinor operator-(DiracSpinor lhs, const DiracSpinor &rhs);

  //! Scalar multiplication
  DiracSpinor &operator*=(const double x);
  friend DiracSpinor operator*(DiracSpinor Fa, const double x);
  friend DiracSpinor operator*(const double x, DiracSpinor Fa);

  //! Multiplication by array (function)
  DiracSpinor &operator*=(const std::vector<double> &v);
  friend DiracSpinor operator*(const std::vector<double> &v, DiracSpinor Fa);

  //! Comparitor overloads (compares n, then kappa):
  friend bool operator==(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator!=(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator<(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator>(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator<=(const DiracSpinor &lhs, const DiracSpinor &rhs);
  friend bool operator>=(const DiracSpinor &lhs, const DiracSpinor &rhs);

  //! Custom comparitors (for sorting): l, j, kappa_index, energy
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
    return lhs.en() < rhs.en();
  }

  // Static (helper) functions:

  //! Returns worst |<a|b>| (or |<a|b>-1| for a=b) {val, state_names}
  static std::pair<double, std::string>
  check_ortho(const std::vector<DiracSpinor> &a,
              const std::vector<DiracSpinor> &b);

  //! Returns formatted states string (e.g., '7sp5d') given list of orbs
  static std::string state_config(const std::vector<DiracSpinor> &orbs);

  //! Constructs H-like (pointlike) DiracSpinor - mainly for testing
  static DiracSpinor exactHlike(int n, int k, std::shared_ptr<const Grid> rgrid,
                                double zeff, double alpha = 0.0);

  //! Returns maximum (2j) found in {orbs}
  static int max_tj(const std::vector<DiracSpinor> &orbs);
  //! Returns maximum l found in {orbs}
  static int max_l(const std::vector<DiracSpinor> &orbs);
  //! Returns maximum kappa_index found in {orbs}
  static int max_kindex(const std::vector<DiracSpinor> &orbs);

  DiracSpinor &operator=(const DiracSpinor &);
  DiracSpinor(const DiracSpinor &) = default;
  ~DiracSpinor() = default;
};
