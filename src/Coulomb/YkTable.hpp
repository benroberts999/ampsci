#pragma once
#include "Angular/Angular_tables.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <utility>
#include <vector>

namespace Coulomb {

constexpr bool check_bounds = true;

//! @brief Calculates + stores Hartree Y functions + Angular (w/ look-up)
/*! @details

Also stores an Angular::Angular::Ck_ab table (publically).

Note: Look-up must be done such that: Fa MUST be member of a_orbitals, Fb of
b_orbitals. Turn the compile-time option "check_bounds" to true (in YkTable.hpp)
to turn on run-time bounds checking.

Definitions:

\f[y^k_{ij}(r) = \int_0^\infty \frac{r_<^k}{r_>^{k+1}}\rho_{ij}(r')\,{\rm
d}r'\f]

\f[\rho(r) = f_i(r)f_j(r) + g_i(r)g_j(r)\f]

with \f$r_< = min(r,r')\f$
*/
class YkTable {
public:
  YkTable(const Grid *const in_grid,
          const std::vector<DiracSpinor> *const in_a_orbs,
          const std::vector<DiracSpinor> *const in_b_orbs = nullptr);
  YkTable &operator=(const YkTable &) = delete;
  YkTable(const YkTable &) = delete;
  ~YkTable() = default;

private:
  const std::vector<DiracSpinor> *const m_a_orbs;
  const std::vector<DiracSpinor> *const m_b_orbs;
  const Grid *const m_grid;
  const bool m_aisb;

private:
  std::size_t a_size = 0;
  std::size_t b_size = 0;
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_abkr = {};
  Angular::Ck_ab m_Ck = Angular::Ck_ab();

public:
  //! Re-calculates all y^k integrals [happens automatically on construct]
  void update_y_ints();
  //! Re-calculates y^k integrals involving single orbital Fn
  void update_y_ints(const DiracSpinor &Fn);

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

  //! Returns maximum value of 2j in {a} and {b} orbitals
  int max_tj() const;

private:
  // Re-sizes the m_y_abkr [and a_size, b_size]
  void resize_y();

public:
  //! Returns min and max k (multipolarity) allowed for C^k_ab
  static std::pair<int, int> k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b);
};

} // namespace Coulomb
