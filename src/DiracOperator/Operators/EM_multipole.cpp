#include "DiracOperator/Operators/EM_multipole.hpp"

namespace DiracOperator {

//==============================================================================
// VEk_Len
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

//------------------------------------------------------------------------------
double VEk_Len::radialIntegral(const DiracSpinor &Fa,
                               const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
    return 0.0;
  }

  const auto K = double(m_K);
  const auto cc = double(Fa.kappa() - Fb.kappa()) / (K + 1.0);
  const auto cx = std::sqrt((K + 1.0) / K);

  return cx * (Rab(+1, *p_jK, Fa, Fb) - (cc + 1.0) * Vab(*p_jKp1, Fa, Fb) +
               (1.0 - cc) * Wab(*p_jKp1, Fa, Fb));
}

//------------------------------------------------------------------------------
void VEk_Len::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    // nb: may not be exact! Ensure lookup table is dense enough!
    p_jK = &m_jl->jL_nearest(m_K, q);
    p_jKp1 = &m_jl->jL_nearest(m_K + 1, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &m_jKp1);
    p_jK = &m_jK;
    p_jKp1 = &m_jKp1;
  }
}

//==============================================================================
// VEk
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

//------------------------------------------------------------------------------
double VEk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
    return 0.0;
  }

  const auto K = double(m_K);
  assert(m_K != 0); // should already be discounted!
  const auto cx = std::sqrt((K + 1.0) / K);
  const auto dk = double(Fa.kappa() - Fb.kappa());

  const auto Pp1 = Pab(+1, *p_jK_on_qr, Fa, Fb);
  const auto Pp2 = Pab(+1, *p_jKp1, Fa, Fb);
  const auto Pm1 = Pab(-1, *p_jK_on_qr, Fa, Fb);

  return cx * (dk * (Pp1 - Pp2 / (K + 1)) - K * Pm1);
}

//------------------------------------------------------------------------------
void VEk::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    // nb: may not be exact! Ensure lookup table is dense enough!
    p_jK_on_qr = &m_jl->jL_on_qr_nearest(m_K, q);
    p_jKp1 = &m_jl->jL_nearest(m_K + 1, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &m_jKp1);
    using namespace qip::overloads;
    m_jK_on_qr /= (q * m_vec);
    p_jK_on_qr = &m_jK_on_qr;
    p_jKp1 = &m_jKp1;
  }
}

//==============================================================================
// VLk
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

//------------------------------------------------------------------------------
double VLk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  const auto K = double(m_K);
  const auto dk = double(Fa.kappa() - Fb.kappa());

  return -dk * Pab(+1, *p_jK_on_qr, Fa, Fb) + K * Pab(-1, *p_jK_on_qr, Fa, Fb) -
         Pab(-1, *p_jKp1, Fa, Fb);
}

//------------------------------------------------------------------------------
void VLk::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    // nb: may not be exact! Ensure lookup table is dense enough!
    p_jK_on_qr = &m_jl->jL_on_qr_nearest(m_K, q);
    p_jKp1 = &m_jl->jL_nearest(m_K + 1, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &m_jKp1);
    using namespace qip::overloads;
    m_jK_on_qr /= (q * m_vec);
    p_jK_on_qr = &m_jK_on_qr;
    p_jKp1 = &m_jKp1;
  }
}

//==============================================================================
// VMk
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

//------------------------------------------------------------------------------
double VMk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa()) ||
      m_K == 0) {
    return 0.0;
  }

  const auto K = double(m_K);
  const auto sk = double(Fa.kappa() + Fb.kappa());
  assert(m_K != 0); // should already be discounted!
  const auto ck = sk / std::sqrt(K * (K + 1.0));

  return -ck * Pab(+1, *p_jK, Fa, Fb);
}

//------------------------------------------------------------------------------
void VMk::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);
  if (m_jl) {
    // nb: may not be exact! Ensure lookup table is dense enough!
    p_jK = &m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
    p_jK = &m_jK;
  }
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

//------------------------------------------------------------------------------
double Phik::radialIntegral(const DiracSpinor &Fa,
                            const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  return Rab(+1, *p_jK, Fa, Fb);
}

//------------------------------------------------------------------------------
void Phik::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);
  if (m_jl) {
    // nb: may not be exact! Ensure lookup table is dense enough!
    p_jK = &m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
    p_jK = &m_jK;
  }
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

//------------------------------------------------------------------------------
double Sk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  return Rab(-1, *p_jK, Fa, Fb);
}

//------------------------------------------------------------------------------
void Sk::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);
  if (m_jl) {
    // nb: may not be exact! Ensure lookup table is dense enough!
    p_jK = &m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
    p_jK = &m_jK;
  }
}

//==============================================================================
// AEk::radial_rhs
//==============================================================================

DiracSpinor AEk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

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

//------------------------------------------------------------------------------
double AEk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
    return 0.0;
  }

  const auto K = double(m_K);
  assert(m_K != 0); // should already be discounted!
  const auto cx = std::sqrt((K + 1.0) / K);
  const auto dk = double(Fa.kappa() + Fb.kappa());

  const auto Pp1 = Rab(-1, jK_on_qr, Fa, Fb);
  const auto Pp2 = Rab(-1, jKp1, Fa, Fb);
  const auto Pm1 = Rab(+1, jK_on_qr, Fa, Fb);

  return cx * (dk * (Pp1 - Pp2 / (K + 1)) - K * Pm1);
}

//------------------------------------------------------------------------------
void AEk::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    jK_on_qr = m_jl->jL_on_qr_nearest(m_K, q);
    jKp1 = m_jl->jL_nearest(m_K + 1, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jK_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &jKp1);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      jK_on_qr[i] /= (q * m_vec[i]);
    }
  }
}

//------------------------------------------------------------------------------
DiracSpinor ALk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

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

//------------------------------------------------------------------------------
double ALk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  const auto K = double(m_K);
  const auto dk = double(Fa.kappa() + Fb.kappa());

  return -dk * Rab(-1, jK_on_qr, Fa, Fb) + K * Rab(+1, jK_on_qr, Fa, Fb) -
         Rab(+1, jKp1, Fa, Fb);
}

//------------------------------------------------------------------------------
void ALk::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    jK_on_qr = m_jl->jL_on_qr_nearest(m_K, q);
    jKp1 = m_jl->jL_nearest(m_K + 1, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jK_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &jKp1);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      jK_on_qr[i] /= (q * m_vec[i]);
    }
  }
}

//==============================================================================
// AMk::radial_rhs
//==============================================================================

DiracSpinor AMk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

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

//------------------------------------------------------------------------------
double AMk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa()) ||
      m_K == 0) {
    return 0.0;
  }

  const auto K = double(m_K);
  const auto sk = double(Fa.kappa() - Fb.kappa());
  assert(m_K != 0); // should already be discounted!
  const auto ck = sk / std::sqrt(K * (K + 1.0));

  return -ck * Rab(-1, jK, Fa, Fb);
}

//------------------------------------------------------------------------------
void AMk::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    jK = m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jK);
  }
}

//==============================================================================
// Phi5k::radial_rhs
//==============================================================================

DiracSpinor Phi5k::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

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

//------------------------------------------------------------------------------
double Phi5k::radialIntegral(const DiracSpinor &Fa,
                             const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  return Pab(-1, jk, Fa, Fb);
}

//------------------------------------------------------------------------------
void Phi5k::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    jk = m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }
}

//==============================================================================
// S5k::radial_rhs
//==============================================================================

DiracSpinor S5k::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

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

//------------------------------------------------------------------------------
double S5k::radialIntegral(const DiracSpinor &Fa,
                             const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  return Pab(+1, jk, Fa, Fb);
}

//------------------------------------------------------------------------------
void S5k::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    jk = m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }
}

} // namespace DiracOperator