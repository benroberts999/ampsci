#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <unordered_map>
#include <vector>

namespace Coulomb {

//! @brief Calculates + stores Hartree Y functions + Angular (w/ look-up),
//! taking advantage of symmetry
/*! @details

Also stores a Ck and 6J table

Definitions:

\f[
y^k_{ij}(r) = \int_0^\infty \frac{r_<^k}{r_>^{k+1}}\rho_{ij}(r')\,{\rm d}r'
\f]

\f[\rho(r) = f_i(r)f_j(r) + g_i(r)g_j(r)\f]

with \f$r_< = min(r,r')\f$
*/
class YkTable {

private:
  std::vector<std::unordered_map<uint32_t, std::vector<double>>> m_Y{};
  Angular::CkTable m_Ck{};
  Angular::SixJTable m_6j{};

public:
  YkTable(const std::vector<DiracSpinor> &a_orbs,
          const std::vector<DiracSpinor> &b_orbs) {
    calculate(a_orbs, b_orbs);
  }
  YkTable(const std::vector<DiracSpinor> &a_orbs) { calculate(a_orbs); }
  YkTable() {}

  //----------------------------------------------------------------------------
  //! Re-calculates all y_ab functions (will over-ride existing ones); only
  //! calculates for a in as, b in bs
  //!@details Note: will re-calculate all, so don't use this to 'extend' the
  //! table
  void calculate(const std::vector<DiracSpinor> &as,
                 const std::vector<DiracSpinor> &bs);
  //! Re-calculates all y_ij functions (will over-ride existing ones) [i and j
  //! in as]
  void calculate(const std::vector<DiracSpinor> &as) { calculate(as, as); }

  //! Extends the Ck and 6J tables up to new maximum 2*j
  void extend_angular(int new_max_2j);

  //! Returns a (const ref) to Ck table [see Angular::CkTable]
  const Angular::CkTable &Ck() const { return m_Ck; }
  //! Returns a (const ref) to SixJ table [see Angular::SixJTable]
  const Angular::SixJTable &SixJ() const { return m_6j; }

  //! Returns a pointer to constant vector y^k_ab. If that integral is not
  //! stores, returns nullptr
  const std::vector<double> *get(const int k, const DiracSpinor &Fa,
                                 const DiracSpinor &Fb) const;

  //! Calculates Qk using the existing yk integrals. Note: Yk and Ck tables
  //! *must* include all required values, or behaviour not defined.
  [[nodiscard]] double Q(const int k, const DiracSpinor &Fa,
                         const DiracSpinor &Fb, const DiracSpinor &Fc,
                         const DiracSpinor &Fd) const;

  //! Calculates Pk using the existing yk integrals. Note: Yk and Ck tables
  //! *must* include all required values, or behaviour not defined.
  [[nodiscard]] double P(const int k, const DiracSpinor &Fa,
                         const DiracSpinor &Fb, const DiracSpinor &Fc,
                         const DiracSpinor &Fd) const;

  //! Calculates Wk=Qk+Pk using the existing yk integrals. Note: Yk and Ck
  //! tables *must* include all required values, or behaviour not defined.
  [[nodiscard]] double W(const int k, const DiracSpinor &Fa,
                         const DiracSpinor &Fb, const DiracSpinor &Fc,
                         const DiracSpinor &Fd) const;

  //! Calculates Q^K(v)_bcd using existing yk integrals
  [[nodiscard]] DiracSpinor Qkv_bcd(int kappa, const DiracSpinor &Fb,
                                    const DiracSpinor &Fc,
                                    const DiracSpinor &Fd, const int k) const;

  //! Calculates P^K(v)_bcd using existing yk integrals, including (optional)
  //! screening factors
  [[nodiscard]] DiracSpinor Pkv_bcd(int kappa, const DiracSpinor &Fb,
                                    const DiracSpinor &Fc,
                                    const DiracSpinor &Fd, const int k,
                                    const std::vector<double> &f2k = {}) const;

private:
  // Allocates space for the Yk table, but does not calculate Yk. This is
  // because allocation cannot be done in parallel, but once allocation is done,
  // calculation can be done in //
  void allocate_space(const std::vector<DiracSpinor> &a_orbs,
                      const std::vector<DiracSpinor> &b_orbs);

  // Returns key used for look-up table (unordered_map), considers symmetry
  uint32_t ab_key(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  // Returns reference to yk vector at 'key'; if no such vector exists, creates
  // it first.
  std::vector<double> &get_or_insert(std::size_t k, uint32_t key);

  // Returns reference to yk vector at 'key'. This vector *must* exist already
  std::vector<double> &get_ref(const int k, const DiracSpinor &Fa,
                               const DiracSpinor &Fb);
};

} // namespace Coulomb
