#pragma once
#include "Coulomb/include.hpp"
#include "Coulomb/meTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <optional>
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

/*! @brief Calculates Structure Radiation + Normalisation of states, using diagram
method.

  @details
  Two main functions: reducedME(), norm(), which retrun the reduced matrix elements
  for structure radiation and normalisation, respectively.

  @note
  For Structure Radiation, must call StructureRad::solve_core() first.

  Structure radiation is sum of srTB()+srC(); norm() gives normalisation of
  states.

  * Typically, one should use splines for the 'legs' (outer states) of diagrams.
  However, code is written such that {w,v} (the valence states) always appear as
  the _first_ state in the Q integrals; thus those states are always directly
  integrated (other states may use the existing y^k integrals, which were
  calculated using the spline states) * unless using QkTable (see below).

  * For using spline states as legs: This means you must typically ensure basis
  is large enough to make relevant valence states "physical"

  * The user should check if the zeroth order <w|t|v> matches closely between
  valence/spline states. If not, basis/cavity is too small and SR+N probably
  meaningless

  * Option to use QkTable to calculate SR+N. Otherwise, calculates Coulomb
  integrals on-the-fly. If QkTable is used, leads to ~10x speed-up, but uses large
  amount of memory. Also, using QkTable means spline states are used for "legs" of
  diagrams.

  * Norm() works differently; It directly taken in the operator to calculate norm.
  This is because norm depends only on Coulomb operators .
  In some sense, would make sense to keep norm in a seperate class.
  However, they are typically calculated together, and use the same Coulomb integrals.

  * Formulas are presented in (for example):
    * W. R. Johnson, Z. W. Liu, and J. Sapirstein, [At. Data Nucl. Data Tables 64,
  279 (1996)](https://www.sciencedirect.com/science/article/pii/S0092640X96900248).

  * Specific form of formulas used here are presented in: 
    [Phys. Rev. A **107**, 052812 (2023)](https://arxiv.org/abs/2211.11134)

*/
class StructureRad {

public:
  //! Construct a Structeure Radiation (MBPT) object
  /*! @details
    en_core defines the core boundary: states with e < en_core are treated as
    core states, while states with e > en_core are treated as excited states.
    Typically
    en_core = max(e_core) - min(e_valence) = wf.FermiLevel().

    nminmax is the pair {nmin, nmax}. Only core states with n >= nmin are used,
    and only excited states with n <= nmax are included in the summations.

    Qk_fname is the filename for the QkTable. If provided, the QkTable will be
    read from / written to this file and then used to calculate SR. This leads
    to ~10× speedup in the calculation at the cost of significantly higher
    memory usage.

    fk and etak are effective screening factors for screening and hole-particle.
    These should only be used as a test.
  */
  StructureRad(const std::vector<DiracSpinor> &basis, double en_core,
               std::pair<int, int> nminmax = {0, 999},
               const std::string &Qk_fname = "",
               const std::vector<double> &fk = {},
               const std::vector<double> &etak = {}, bool verbose = true);

private:
  // Switch to use Qk (rather than Yk) table
  bool m_use_Qk;
  // effective screening; sqrt - add to both sides!
  std::vector<double> m_root_fk;
  // effective hole-particle
  std::vector<double> m_etak;

  Coulomb::YkTable mY{};
  std::optional<Coulomb::QkTable> mQ{std::nullopt};
  // Keeping a local copy is significantly faster
  std::vector<DiracSpinor> mCore{}, mExcited{}, mBasis{};

  // matrix element table; set by solve_core
  Coulomb::meTable<double> mTab{};
  // Rank of current operator; set by solve_core
  int m_K{-1};

public:
  //! const reference to employed _subset_ of core orbitals
  const std::vector<DiracSpinor> &core() const { return mCore; }
  //! const reference to employed _subset_ of excited orbitals
  const std::vector<DiracSpinor> &excited() const { return mExcited; }
  //! const reference to employed _subset_ of core+excited orbitals
  const std::vector<DiracSpinor> &basis() const { return mBasis; }
  //! const reference to underlying matrix element table. May be empty
  const Coulomb::meTable<double> &me_Table() const { return mTab; }

  //! const reference to Yk table. NOTE: may not be initialised!
  const Coulomb::YkTable &Yk() const { return mY; }
  //! const reference to underlying Qk table. Optional; may not be set
  const std::optional<Coulomb::QkTable> &Qk() const { return mQ; }

  //! Prepares
  /*! @details
    - Pre-computes the full set of <i||T||j> matrix elements
    - Does for all core-core, core-excited, excited-excited orbitals
    - Optionally includes core polarisation; dV may be null
    - h must not be null
    - For frequency-dependent operators, this assumes operator has already been updated
    - Assumes Core Polarisation has already been solves

    @note
    This must be called before structure radiation can be calculated
  */
  void solve_core(const DiracOperator::TensorOperator *const h,
                  const ExternalField::CorePolarisation *const dV = nullptr);

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

  //! Reduced matrix element of structure radiation **only**. Does not include norm?
  double reducedME(const DiracOperator::TensorOperator *const h,
                   const DiracSpinor &w, const DiracSpinor &v,
                   const ExternalField::CorePolarisation *const dV,
                   double omega) const {
    return srTB(h, w, v, omega, dV).first + srC(h, w, v, dV).first;
  }

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
