#pragma once
#include "Coulomb/meTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <cassert>
#include <utility>
#include <vector>

namespace HF {

// namespace hidden {
// struct Breit_Bk_ba; // forward decl
// }

//==============================================================================
namespace Breit_gb {

struct single_k_mop {
  // Class to hold the Breit-Coulomb integrals, for single k
public:
  single_k_mop() {}
  single_k_mop(const DiracSpinor &Fi, const DiracSpinor &Fj, int k) {
    calculate(Fi, Fj, k);
  }

  void calculate(const DiracSpinor &Fi, const DiracSpinor &Fj, int k);
  std::vector<double> b0_minus{}, bi_minus{};
  std::vector<double> b0_plus{}, bi_plus{};
  std::vector<double> g0_minus{}, gi_minus{};
  std::vector<double> g0_plus{}, gi_plus{};
};

struct single_k_n {
  // Class to hold the Breit-Coulomb integrals, for single k
public:
  single_k_n() {}
  single_k_n(const DiracSpinor &Fi, const DiracSpinor &Fj, int k) {
    calculate(Fi, Fj, k);
  }

  void calculate(const DiracSpinor &Fi, const DiracSpinor &Fj, int k);
  std::vector<double> g{};

private:
  std::vector<double> gi{};
};

} // namespace Breit_gb

//==============================================================================
//! Breit (Hartree-Fock Breit) interaction potential
class Breit {
  double m_scale;

  // Scaling factors for each term (mainly for tests). {M,N}=Gaunt; {O,P} rtrd
  double m_M{1.0}, m_N{1.0}, m_O{1.0}, m_P{1.0};

  // For speedy lookup, when only full integral is required
  std::vector<Coulomb::meTable<Breit_gb::single_k_mop>> m_gb{};
  std::vector<Coulomb::meTable<Breit_gb::single_k_n>> m_gb_N{};

public:
  //! Constructs Breit operator; scale is overall scaling factor, 1.0 typical
  Breit(double in_scale = 1.0) : m_scale(in_scale) {}

  //! Allows one to update scaling factor(s)
  void update_scale(double t_scale = 1.0, double t_M = 1.0, double t_N = 1.0,
                    double t_O = 1.0, double t_P = 1.0) {
    m_scale = t_scale;
    m_M = t_M;
    m_N = t_N;
    m_O = t_O;
    m_P = t_P;
  }

  static bool Bk_SR(int k, const DiracSpinor &v, const DiracSpinor &w,
                    const DiracSpinor &x, const DiracSpinor &y) {

    const auto have_mop = Angular::Ck_kk_SR(k, v.kappa(), x.kappa()) &&
                          Angular::Ck_kk_SR(k, w.kappa(), y.kappa()) &&
                          v != x && w != y;

    const auto have_n = Angular::Ck_kk_SR(k, -v.kappa(), x.kappa()) &&
                        Angular::Ck_kk_SR(k, -w.kappa(), y.kappa()) &&
                        (v.kappa() + x.kappa() != 0) &&
                        (w.kappa() + y.kappa() != 0);

    return (have_mop || have_n);
  }

  static std::pair<int, int> k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b,
                                      const DiracSpinor &c,
                                      const DiracSpinor &d) {
    const auto [k1, k2] = Coulomb::k_minmax_tj(a.twoj(), c.twoj());
    const auto [k3, k4] = Coulomb::k_minmax_tj(b.twoj(), d.twoj());
    return {std::max({1, k1, k3}), std::min(k2, k4)};
  }

  static std::pair<int, int> k_minmax_tj(int tja, int tjb, int tjc, int tjd) {
    const auto [k1, k2] = Coulomb::k_minmax_tj(tja, tjc);
    const auto [k3, k4] = Coulomb::k_minmax_tj(tjb, tjd);
    return {std::max({1, k1, k3}), std::min(k2, k4)};
  }

  void fill_gb(const std::vector<DiracSpinor> &basis, int t_max_k = 99);

  //! Resturns scaling factor
  double scale_factor() const { return m_scale; };

  //! Calculates V_br*Fa [Breit part of HF-Breit pot.]
  DiracSpinor VbrFa(const DiracSpinor &Fa,
                    const std::vector<DiracSpinor> &core) const;

  //! dV_b*Fa, dV_b is the RPA correction arising due to Fb -> Fb + dFb
  //! @details
  //! K is multipolarity of RPA operator, Fb is core state, with Xbeta and Ybeta
  //! perturbations. "reduced rhs"
  DiracSpinor dV_Br(int kappa, int K, const DiracSpinor &Fa,
                    const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                    const DiracSpinor &Ybeta) const;

  //! Reduced Breit integral (analogue of Coulomb Q^k).
  double Bk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Reduced Breit integral (analogue of Coulomb Q^k), faster version.
  //! Can only use if fill_gb() has been called
  double Bk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                   const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Reduced exchange Breit integral (analogue of Coulomb P^k).
  double BPk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                  const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Reduced exchange Breit integral (analogue of Coulomb P^k), faster version.
  //! Can only use if fill_gb() has been called
  double BPk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Reduced anti-symmetrised Breit integral (analogue of Coulomb W^k).
  //! BWk_abcd = Bk_abcd + BPk_abcd
  double BWk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                  const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Reduced anti-symmetrised Breit integral (analogue of Coulomb W^k), faster.
  //! Can only use if fill_gb() has been called.
  double BWk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Bkv_bcd defined: B^k_abcd = <a|Bkv_bcd> (direct part)
  DiracSpinor Bkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                      const DiracSpinor &Fc, const DiracSpinor &Fd) const;
  //! BPkv_bcd defined: P(B)^k_abcd = <a|BPkv_bcd> (exchange part)
  DiracSpinor BPkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                       const DiracSpinor &Fc, const DiracSpinor &Fd) const;
  //! BWkv_bcd defined: W(B)^k_abcd = B^k_abcd + P(B)^k_abcd= <a|BWkv_bcd>
  DiracSpinor BWkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                       const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  //! Hartree-Fock (one-particle) part of second-order Breit correction.
  //! With Breit-Coulomb Hartree-Fock, is included automatically.
  double de2_HF(const DiracSpinor &Fv, const std::vector<DiracSpinor> &holes,
                const std::vector<DiracSpinor> &excited) const;

  //! Sigma (two-particle) part of second-order Breit correction.
  double de2(const DiracSpinor &Fv, const std::vector<DiracSpinor> &holes,
             const std::vector<DiracSpinor> &excited) const;
};

} // namespace HF
