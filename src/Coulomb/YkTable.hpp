#pragma once
#include "Angular/Angular_tables.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <utility>
#include <vector>

namespace Coulomb {

constexpr static bool check_bounds = true;

//! @brief Calculates + stores Hartree Y functions + Angular (w/ look-up)
/*! @details

Stores an Angular::Angular::Ck_ab table (publically)

Definitions:
  - y^k_ij(r)   := Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
  - Lambda^k_ij := 3js((ji,jj,k),(-1/2,1/2,0))^2 * parity(li+lj+k)
*/
class YkTable {
public:
  YkTable(const Grid *const in_grid,
          const std::vector<DiracSpinor> *const in_a_orbs,
          const std::vector<DiracSpinor> *const in_b_orbs = nullptr);
  YkTable &operator=(const YkTable &) = default;
  YkTable(const YkTable &) = default;
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
  //! Re-calculates all y^k integrals
  void update_y_ints();
  //! Re-calculates y^k integrals involving single orbital Fn
  void update_y_ints(const DiracSpinor &Fn);

  //! Returns single y^k_ab(r) [by const ref]
  const std::vector<double> &get_yk_ab(const int k, const DiracSpinor &Fa,
                                       const DiracSpinor &Fb) const;
  //! @brief Returns vector of y^k_ab(r) for each k [by const ref]
  //! @details Note: vector starts at first non-zero k (not zero)! use k_minmax
  //! to get k_min
  const std::vector<std::vector<double>> &get_y_ab(const DiracSpinor &Fa,
                                                   const DiracSpinor &Fb) const;

  //! Access C^k angular tables (have public access to const methods only)
  const Angular::Ck_ab &Ck() const { return m_Ck; }

private:
  // Returns maximum value of 2j in {a} and {b} orbitals
  int max_tj() const;
  // Re-sizes the m_y_abkr [and a_size, b_size]
  void resize_y();

public:
  //! Returns min and max k (multipolarity) allowed for C^k_ab
  static std::pair<int, int> k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b);
};

} // namespace Coulomb
