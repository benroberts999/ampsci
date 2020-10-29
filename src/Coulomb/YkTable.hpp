#pragma once
#include "Angular/Angular_tables.hpp"
#include <memory>
#include <utility>
#include <vector>
class Grid;
class DiracSpinor;

namespace Coulomb {

//! @brief Calculates + stores Hartree Y functions + Angular (w/ look-up)
/*! @details

Also stores an Angular::Angular::Ck_ab table (publically).

Note: Look-up must be done such that: Fa MUST be member of a_orbitals, Fb of
b_orbitals. Cheking is done with assers, so only in 'dev' or 'debug' mode -
ranges not checked in release mode.

Definitions:

\f[y^k_{ij}(r) = \int_0^\infty \frac{r_<^k}{r_>^{k+1}}\rho_{ij}(r')\,{\rm
d}r'\f]

\f[\rho(r) = f_i(r)f_j(r) + g_i(r)g_j(r)\f]

with \f$r_< = min(r,r')\f$
*/
class YkTable {
public:
  YkTable(std::shared_ptr<const Grid> in_grid,
          const std::vector<DiracSpinor> *const in_a_orbs,
          const std::vector<DiracSpinor> *const in_b_orbs = nullptr);
  YkTable &operator=(const YkTable &) = delete;
  YkTable(const YkTable &) = delete;
  ~YkTable() = default;

private:
  const std::vector<DiracSpinor> *const m_a_orbs;
  const std::vector<DiracSpinor> *const m_b_orbs;
  std::shared_ptr<const Grid> m_grid;
  const bool m_aisb;

private:
  std::size_t a_size = 0;
  std::size_t b_size = 0;
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_abkr = {};
  Angular::Ck_ab m_Ck = Angular::Ck_ab();
  Angular::SixJ m_6j = Angular::SixJ();

public:
  //! Re-calculates all y^k integrals [happens automatically on construct]
  void update_y_ints();
  //! Re-calculates y^k integrals involving single orbital Fn
  void update_y_ints(const DiracSpinor &Fn);

  const std::vector<DiracSpinor> &get_a() const { return *m_a_orbs; }
  const std::vector<DiracSpinor> &get_b() const { return *m_b_orbs; }

  //! Returns single y^k_ab(r) [by const ref]
  //! @details NOTE: Fa MUST be member of a_orbitals, Fb of b_or
  const std::vector<double> &get_yk_ab(const int k, const DiracSpinor &Fa,
                                       const DiracSpinor &Fb) const;
  //! Returns a pointer to (const) y^k_ab(r), if k valid. Else, returns nullptr
  const std::vector<double> *ptr_yk_ab(const int k, const DiracSpinor &Fa,
                                       const DiracSpinor &Fb) const {
    const auto [min, max] = k_minmax(Fa, Fb);
    if (k >= min && k <= max) {
      return &get_yk_ab(k, Fa, Fb);
    }
    return nullptr;
  }

  //! overload for get_yk_ab()
  const std::vector<double> &operator()(const int k, const DiracSpinor &Fa,
                                        const DiracSpinor &Fb) const {
    return get_yk_ab(k, Fa, Fb);
  }
  //! @brief Returns vector of y^k_ab(r) for each k [by const ref]
  //! @details Note: vector starts at first non-zero k (not zero)! use k_minmax
  //! to get k_min.
  //! NOTE: Fa MUST be member of a_orbitals, Fb of b_or
  const std::vector<std::vector<double>> &get_y_ab(const DiracSpinor &Fa,
                                                   const DiracSpinor &Fb) const;
  const std::vector<std::vector<double>> &
  operator()(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    return get_y_ab(Fa, Fb);
  }

  //! Access C^k angular tables (have public access to const methods only)
  const Angular::Ck_ab &Ck() const { return m_Ck; }
  void extend_Ck(const int max_twoj) { m_Ck.fill(max_twoj); }

  //! Returns maximum value of 2j in {a} and {b} orbitals
  int max_tj() const;

private:
  // Re-sizes the m_y_abkr [and a_size, b_size]
  void resize_y();

public:
  //! Returns min and max k (multipolarity) allowed for C^k_ab
  static std::pair<int, int> k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b);

  //! Returns Qk, Pk, Wk - uses existing yk table half the Coulomb integral.
  //! @details
  //! These are only safe to call if Y^k_ab {a}={b} i.e., constructed with a
  //! single set of states! NOTE: There is no safety check, undefined behaviour
  //! otherwise.
  double Qk(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
            const DiracSpinor &Fc, const DiracSpinor &Fd) const;
  double Pk(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
            const DiracSpinor &Fc, const DiracSpinor &Fd) const;
  double Wk(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
            const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Returns min and max k (multipolarity) allowed for Q^k_abcd, Wk, Pk. For Q
  //! (and Q only), parity rule is included, so you may safely call k+=2
  static std::pair<int, int> k_minmax_Q(const DiracSpinor &a,
                                        const DiracSpinor &b,
                                        const DiracSpinor &c,
                                        const DiracSpinor &d);
  //! DOES NOT contain parity rules (6j only) - so NOT safe to call k+=2
  static std::pair<int, int> k_minmax_P(const DiracSpinor &a,
                                        const DiracSpinor &b,
                                        const DiracSpinor &c,
                                        const DiracSpinor &d);
  //! DOES NOT contain parity rules (6j only) - so NOT safe to call k+=2
  static std::pair<int, int> k_minmax_W(const DiracSpinor &a,
                                        const DiracSpinor &b,
                                        const DiracSpinor &c,
                                        const DiracSpinor &d);

  //! "Magic" sixj symbol calculator. Integers (for k) and Spinors
  // do not multiply by 2!
  template <class A, class B, class C, class D, class E, class F>
  static double sixj(A a, B b, C c, D d, E e, F f) {
    return Angular::sixj_2(twojk(a), twojk(b), twojk(c), twojk(d), twojk(e),
                           twojk(f));
  }

private:
  // returns a.twoj() if a is a spinor, or 2*a if a is an integer
  // used for sixj symbol, just to simplify call site!
  template <class A> static int twojk(A a) {
    if constexpr (std::is_same_v<A, DiracSpinor>) {
      return a.twoj();
    } else {
      static_assert(std::is_same_v<A, int>);
      return 2 * a;
    }
  }
};

} // namespace Coulomb
