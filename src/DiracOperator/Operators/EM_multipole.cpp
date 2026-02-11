#include "DiracOperator/Operators/EM_multipole.hpp"

namespace DiracOperator {

//==============================================================================

DiracSpinor VEk_Len::radial_rhs(const int kappa_a,
                                const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();
  if (isZero(kappa_a, Fb.kappa()) || m_K == 0) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto c1 = double(kappa_a - Fb.kappa()) / (K + 1.0);
  const auto cx = std::sqrt((K + 1.0) / K);
  Rab_rhs(+1, *p_jK, &dF, Fb, cx);
  Pab_rhs(+1, *p_jKp1, &dF, Fb, -c1 * cx);
  Pab_rhs(-1, *p_jKp1, &dF, Fb, -cx);
  return dF;
}

} // namespace DiracOperator