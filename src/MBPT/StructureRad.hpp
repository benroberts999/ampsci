#pragma once
#include "Coulomb/include.hpp"
#include "Coulomb/meTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace ExternalField {
class CorePolarisation;
}

//! Many-body perturbation theory
namespace MBPT {

/*!
@brief Calculates Structure Radiation + Normalisation of states, using diagram
method.
@details
Three functions: srTB(), srC(), norm().
Structure radiation is sum of srTB()+srC(); norm() gives normalisation of
states.

 - Typically, one should use splines for the 'legs' (outer states) of diagrams.
 However, code is written such that {w,v} (the valence states) always appear as
the _first_ state in the Q integrals; thus those states are always directly
integrated (other states may use the existing y^k integrals, which were
calculated using the spline states) * unless using QkTable (see below).

 - For using spline states as legs: This means you must typically ensure basis
is large enough to make relevant valence states "physical"

 - The user should check if the zeroth order <w|t|v> matches closely between
valence/spline states. If not, basis/cavity is too small and SR+N probably
meaningless

 - Option to use QkTable to calculate SR+N. Otherwise, calculates Coulomb
integrals on-the-fly. If QkTable is used, leads to ~10x speed-up, but uses large
amount of memory. Also, using QkTable means spline states are used for "legs" of
diagrams.


 - Formulas are presented in (for example):
   * W. R. Johnson, Z. W. Liu, and J. Sapirstein, At. Data Nucl. Data Tables 64,
279 (1996). doi:10.1006/adnd.1996.0024

*/
class StructureRad {

public:
  /*! @details
  en_core:  is defined such that states with e < en_core are in the core,
  while states with e > en_core are not. Typcially:
  en_core = max(e_core)-min(e_valence) =  wf.FermiLevel().
  nminmax is a pair{min, max}: we only used core states with n>=min, and only
  uses excited states with n<=nmax in the summations.
  Qk_fname is filename for QkTable - if given it will read/write QkTable to this
  file, and use QkTable to calculate SR. This ldeads to ~10x speedup in
  calculation, at cost of using much more memory
  */
  StructureRad(const std::vector<DiracSpinor> &basis, double en_core,
               std::pair<int, int> nminmax = {0, 999},
               const std::string &Qk_fname = "",
               const std::vector<double> &fk = {},
               const std::vector<double> &etak = {});

private:
  bool m_use_Qk;
  std::vector<double> m_root_fk; // sqrt - add to both sides!
  std::vector<double> m_etak;
  Coulomb::YkTable mY{};
  std::optional<Coulomb::QkTable> mQ{std::nullopt}; //?
  // nb: it seems conter-intuative, but this copy makes it FASTER!
  std::vector<DiracSpinor> mBasis{}, mCore{}, mExcited{};
  Coulomb::meTable<double> m_tab{};

public:
  const std::vector<DiracSpinor> &core() const { return mCore; }
  const std::vector<DiracSpinor> &excited() const { return mExcited; }

  //! Returns reference to Yk table. NOTE: may not be initialised!
  const Coulomb::YkTable &Yk() const { return mY; }
  const std::optional<Coulomb::QkTable> &Qk() const { return mQ; }

  //! Fill a meTable inside the SR class. Must be performed prior to using
  void fill_table(const DiracOperator::TensorOperator *const h,
                  const ExternalField::CorePolarisation *const dV = nullptr,
                  double omega = 0.0);

  //! Copy an existing meTable into SR class
  void fill_table(const Coulomb::meTable<double> &hab) { m_tab = hab; };

  //! Move existing meTable into SR class
  void fill_table(Coulomb::meTable<double> &&hab) { m_tab = std::move(hab); };

  const Coulomb::meTable<double> &hab_table() const { return m_tab; }

  //! Returns sum of Top+Bottom (SR) diagrams, reduced ME: <w||T+B||v>. Returns
  //! a pair: {TB, TB+dV}: second includes RPA (if dV given)
  std::pair<double, double>
  srTB(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
       const DiracSpinor &v, double omega = 0.0,
       const ExternalField::CorePolarisation *const dV = nullptr) const;

  //! Returns Centre (SR) diagrams, reduced ME: <w||C||v>. Returns
  //! a pair: {C, C+dV}: second includes RPA (if dV given)
  std::pair<double, double>
  srC(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
      const DiracSpinor &v,
      const ExternalField::CorePolarisation *const dV = nullptr) const;

  //! Returns Normalisation of states, reduced ME: <w||h||v>_norm. Returns
  //! a pair: {N, N+dV}: second includes RPA (if dV given).
  std::pair<double, double>
  norm(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
       const DiracSpinor &v,
       const ExternalField::CorePolarisation *const dV = nullptr) const;

  //! Returns sum of SR+Norm diagrams, reduced ME: <w||T+B+C+N||v>. Returns
  //! a pair: {SRN, SRN+dV}: second includes RPA (if dV given)
  std::pair<double, double>
  srn(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
      const DiracSpinor &v, double omega = 0.0,
      const ExternalField::CorePolarisation *const dV = nullptr) const {
    const auto [tb, dvtb] = srTB(h, w, v, omega, dV);
    const auto [c, dvc] = srC(h, w, v, dV);
    const auto [n, dvn] = norm(h, w, v, dV);
    return {tb + c + n, dvtb + dvc + dvn};
  }

  //! Effective screening factor for Coulomb lines
  double f_root_scr(int k) const {
    const auto sk = std::size_t(k);
    return sk < m_root_fk.size() ? m_root_fk[sk] : 1.0;
  }

  //! Effective hole-particle factor for polarisation loops
  double eta_hp(int k) const {
    const auto sk = std::size_t(k);
    return sk < m_etak.size() ? m_etak[sk] : 1.0;
  }

  //! Returns Brueckner orbital contribution to the, reduced ME: <w||h||v>_norm.
  //! Returns a pair: {N, N+dV}: second includes RPA (if dV given).
  //! @details from 10.1006/adnd.1996.0024
  std::pair<double, double>
  BO(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
     const DiracSpinor &v,
     const ExternalField::CorePolarisation *const dV = nullptr, double fw = 1.0,
     double fv = 1.0) const;

  //! Calculates <v|Sigma|w>
  double Sigma_vw(const DiracSpinor &v, const DiracSpinor &w) const;

  //! constructs an me table of {srn, srn+dv} for each pair or {a,b}
  Coulomb::meTable<std::pair<double, double>>
  srn_table(const DiracOperator::TensorOperator *const h,
            const std::vector<DiracSpinor> &as,
            const std::vector<DiracSpinor> &tbs = {}, double omega = 0.0,
            const ExternalField::CorePolarisation *const dV = nullptr) const;

private:
  // "Top" diagrams
  double t1234(int k, const DiracSpinor &w, const DiracSpinor &r,
               const DiracSpinor &v, const DiracSpinor &c) const;
  double t1(int k, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t2(int k, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t3(int k, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t4(int k, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;

  // "Bottom" diagrams
  double b1234(int k, const DiracSpinor &w, const DiracSpinor &c,
               const DiracSpinor &v, const DiracSpinor &r) const;

  // "Centre" diagrams (most important)
  double c1(int k, const DiracSpinor &w, const DiracSpinor &a,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double c2(int k, const DiracSpinor &w, const DiracSpinor &a,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double d1(int k, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &m) const;
  double d2(int k, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &m) const;

  // "Normalisation" terms (most important)
  double n1(const DiracSpinor &v) const;
  double n2(const DiracSpinor &v) const;
  double dSigma_dE(const DiracSpinor &v, const DiracSpinor &i,
                   const DiracSpinor &j, const DiracSpinor &k) const;

  std::pair<double, double>
  z_bo(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
       const DiracSpinor &v, bool transpose = false,
       const ExternalField::CorePolarisation *const dV = nullptr) const;

  //
  double Q(int k, const DiracSpinor &a, const DiracSpinor &b,
           const DiracSpinor &c, const DiracSpinor &d) const {
    // inlcudes effective Coulomb screening
    const auto fk = f_root_scr(k);
    return m_use_Qk ? fk * mQ->Q(k, a, b, c, d) : fk * mY.Q(k, a, b, c, d);
  }

  double P(int k, const DiracSpinor &a, const DiracSpinor &b,
           const DiracSpinor &c, const DiracSpinor &d) const {
    // inlcudes effective Coulomb screening
    return m_use_Qk ? mQ->P2(k, a, b, c, d, mY.SixJ(), m_root_fk) :
                      mY.P2(k, a, b, c, d, m_root_fk);
  }

  double W(int k, const DiracSpinor &a, const DiracSpinor &b,
           const DiracSpinor &c, const DiracSpinor &d) const {
    // inlcudes effective Coulomb screening
    return Q(k, a, b, c, d) + P(k, a, b, c, d);
  }
};

} // namespace MBPT
