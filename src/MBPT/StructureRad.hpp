#pragma once
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace HF {
class ExternalField;
}

namespace MBPT {

/*!
Calculates Structure Radiation + Normalisation of states, using diagram method.

@details

 - Note: MUST use splines for 'legs' of diagrams; i.e., for w and v DiracSpinors
in srTB, srC, and norm.
 - This is because the states become "mixed" in the coulomb integrals (i.e.,
inner program does not distingush between {w,v} and the basis states that have
same quantum numbers!

 - Note: This means you must typically ensure basis is large enough to make
relevant valence states "physical"

 - The user should check if the zeroth order <w|t|v> matches closely between
valence/spline states. If not, basis/cavity is too small and SR+N probably
meaningless

 - Formulas are presented in (for example):
   * W. R. Johnson, Z. W. Liu, and J. Sapirstein, At. Data Nucl. Data Tables 64,
279 (1996). doi:10.1006/adnd.1996.0024

*/
class StructureRad {

  /*
  TODO:
  1. ALSO: add new Yk.Qk, Xk, Wk, Zk etc. tests to Coulomb Tests!
  2. Transform from X,Z -> Q,W/P
  3. Try to make 'v' and 'w' appear _first_ (so always uses actual states) (??)
  */

public:
  /*! @details
  en_core:  is defined such that states with e < en_core are in the core,
  while states with e > en_core are not. Typcially:
  en_core = max(e_core)-min(e_valence).
  nminmax is a pair{min, max}: we only used core states with n>=min, and only
  uses excited states with n>=nmax in the summations.
  */
  StructureRad(const std::vector<DiracSpinor> &basis, double en_core,
               std::pair<int, int> nminmax = {0, 999});

private:
  // Store local copy of basis (seems to make faster)
  std::vector<DiracSpinor> mBasis;
  Coulomb::YkTable mY;
  // nb: it seems conter-intuative, but this copy makes it FASTER!
  std::vector<DiracSpinor> mCore{}, mExcited{};

public:
  //! Returns sum of Top+Bottom (SR) diagrams, reduced ME: <w||T+B||v>. Returns
  //! a pair: {TB, TB+dV}: second includes RPA (if dV given)
  std::pair<double, double>
  srTB(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
       const DiracSpinor &v, double omega,
       const HF::ExternalField *const dV = nullptr) const;

  //! Returns Centre (SR) diagrams, reduced ME: <w||C||v>. Returns
  //! a pair: {C, C+dV}: second includes RPA (if dV given)
  std::pair<double, double>
  srC(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
      const DiracSpinor &v, const HF::ExternalField *const dV = nullptr) const;

  //! Returns Normalisation of states, reduced ME: <w||h||v>_norm. Returns
  //! a pair: {N, N+dV}: second includes RPA (if dV given).
  std::pair<double, double>
  norm(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
       const DiracSpinor &v, const HF::ExternalField *const dV = nullptr) const;

private:
  // "Top" diagrams
  double t1(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t2(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t3(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t4(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;

  // "Bottom" diagrams
  double b1(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;
  double b2(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;
  double b3(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;
  double b4(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;

  // "Centre" diagrams (most important)
  double c1(int K, const DiracSpinor &w, const DiracSpinor &a,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double c2(int K, const DiracSpinor &w, const DiracSpinor &a,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double c3(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &m) const;
  double c4(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &m) const;

  // "Normalisation" terms (most important)
  double n1(const DiracSpinor &v) const;
  double n2(const DiracSpinor &v) const;
};

} // namespace MBPT
