#include "DiracOperator/Operators/EM_multipole.hpp"

namespace DiracOperator {

//==============================================================================
// VEk_Len::radial_rhs
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

//==============================================================================
// VEk::radial_rhs
//==============================================================================

DiracSpinor VEk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa()) || m_K == 0) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto dk = double(kappa_a - Fb.kappa());
  assert(m_K != 0); // should already be discounted!
  const auto cx = std::sqrt((K + 1.0) / K);
  if (dk != 0.0) {
    Pab_rhs(+1, *p_jK_on_qr, &dF, Fb, cx * dk);
    Pab_rhs(+1, *p_jKp1, &dF, Fb, -cx * dk / (K + 1.0));
  }
  Pab_rhs(-1, *p_jK_on_qr, &dF, Fb, -cx * K);

  return dF;
}

//==============================================================================
// VLk::radial_rhs
//==============================================================================

DiracSpinor VLk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa())) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto dk = double(kappa_a - Fb.kappa());

  if (dk != 0.0)
    Pab_rhs(+1, *p_jK_on_qr, &dF, Fb, -dk);
  if (K != 0.0)
    Pab_rhs(-1, *p_jK_on_qr, &dF, Fb, K);
  Pab_rhs(-1, *p_jKp1, &dF, Fb, -1.0);
  return dF;
}

//==============================================================================
// VMk::radial_rhs
//==============================================================================

DiracSpinor VMk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa()) || (kappa_a == -Fb.kappa()) || m_K == 0) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto sk = double(kappa_a + Fb.kappa());
  assert(m_K != 0); // should already be discounted!
  const auto ck = sk / std::sqrt(K * (K + 1.0));

  Pab_rhs(+1, *p_jK, &dF, Fb, -ck);
  return dF;
}

//==============================================================================
// Phik::radial_rhs
//==============================================================================

DiracSpinor Phik::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa())) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  Rab_rhs(+1, *p_jK, &dF, Fb);
  return dF;
}

//==============================================================================
// Sk::radial_rhs
//==============================================================================

DiracSpinor Sk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa())) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  Rab_rhs(-1, *p_jK, &dF, Fb);
  return dF;
}

//==============================================================================
// E5k_w::radial_rhs
//==============================================================================

DiracSpinor E5k_w::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa()) || m_K == 0) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto dk = double(kappa_a + Fb.kappa());
  assert(m_K != 0); // should already be discounted!
  const auto cx = std::sqrt((K + 1.0) / K);

  Rab_rhs(-1, jK_on_qr, &dF, Fb, cx * dk);
  Rab_rhs(-1, jKp1, &dF, Fb, -cx * dk / (K + 1.0));
  Rab_rhs(+1, jK_on_qr, &dF, Fb, -cx * K);

  return dF;
}

//==============================================================================
// L5k_w::radial_rhs
//==============================================================================

DiracSpinor L5k_w::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa())) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto dk = double(kappa_a + Fb.kappa());

  Rab_rhs(-1, jK_on_qr, &dF, Fb, -dk);
  Rab_rhs(+1, jK_on_qr, &dF, Fb, K);
  Rab_rhs(+1, jKp1, &dF, Fb, -1.0);
  return dF;
}

//==============================================================================
// M5k_w::radial_rhs
//==============================================================================

DiracSpinor M5k_w::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa()) || (kappa_a == -Fb.kappa()) || m_K == 0) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto sk = double(kappa_a - Fb.kappa());
  assert(m_K != 0); // should already be discounted!
  const auto ck = sk / std::sqrt(K * (K + 1.0));

  Rab_rhs(-1, jK, &dF, Fb, -ck);
  return dF;
}

//==============================================================================
// Phi5k_w::radial_rhs
//==============================================================================

DiracSpinor Phi5k_w::radial_rhs(const int kappa_a,
                                const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa())) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  Pab_rhs(-1, jk, &dF, Fb);
  return dF;
}

//==============================================================================
// S5k_w::radial_rhs
//==============================================================================

DiracSpinor S5k_w::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa())) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  Pab_rhs(+1, jk, &dF, Fb);
  return dF;
}

} // namespace DiracOperator