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
  const auto dk_int = kappa_a + Fb.kappa();
  const auto dk = double(dk_int);
  assert(m_K != 0); // should already be discounted!
  const auto cx = std::sqrt((K + 1.0) / K);

  Rab_rhs(-1, *p_jKp1, &dF, Fb, -cx * dk / (K + 1.0));

  if (dk_int == m_K) {
    // FF terms cancel!
    // R^- - R^+ = -2G
    Gab_rhs(*p_jK_on_qr, &dF, Fb, -2.0 * cx * dk);
    return dF;
  } else {
    Rab_rhs(-1, *p_jK_on_qr, &dF, Fb, cx * dk);
    Rab_rhs(+1, *p_jK_on_qr, &dF, Fb, -cx * K);
    return dF;
  }
}

//------------------------------------------------------------------------------
double AEk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
    return 0.0;
  }

  const auto K = double(m_K);
  assert(m_K != 0); // should already be discounted!
  const auto cx = std::sqrt((K + 1.0) / K);

  const auto dk_int = Fa.kappa() + Fb.kappa();
  const auto dk = double(dk_int);

  const auto Rmp1 = Rab(-1, *p_jKp1, Fa, Fb) / (K + 1);

  if (dk_int == m_K) {
    // FF terms cancel!
    // R^- - R^+ = -2G
    // Always? Or only same kappa?
    const auto GG = Gab(*p_jK_on_qr, Fa, Fb);
    return cx * dk * (-2.0 * GG - Rmp1);
  }

  const auto Rm1 = Rab(-1, *p_jK_on_qr, Fa, Fb);
  const auto Rp1 = Rab(+1, *p_jK_on_qr, Fa, Fb);

  return cx * (dk * Rm1 - K * Rp1 - dk * Rmp1);
}

//------------------------------------------------------------------------------
void AEk::updateFrequency(const double omega) {
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

  Rab_rhs(-1, *p_jK_on_qr, &dF, Fb, -dk);
  Rab_rhs(+1, *p_jK_on_qr, &dF, Fb, K);
  Rab_rhs(+1, *p_jKp1, &dF, Fb, -1.0);
  return dF;
}

//------------------------------------------------------------------------------
double ALk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  const auto K = double(m_K);
  const auto dk = double(Fa.kappa() + Fb.kappa());

  return -dk * Rab(-1, *p_jK_on_qr, Fa, Fb) + K * Rab(+1, *p_jK_on_qr, Fa, Fb) -
         Rab(+1, *p_jKp1, Fa, Fb);
}

//------------------------------------------------------------------------------
void ALk::updateFrequency(const double omega) {
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
// AMk::radial_rhs
//==============================================================================

DiracSpinor AMk::radial_rhs(const int kappa_a, const DiracSpinor &Fb) const {

  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa()) || (kappa_a == Fb.kappa()) || m_K == 0) {
    dF.min_pt() = 0;
    dF.max_pt() = 0;
    return dF;
  }

  const auto K = double(m_K);
  const auto sk = double(kappa_a - Fb.kappa());
  assert(m_K != 0); // should already be discounted!
  const auto ck = sk / std::sqrt(K * (K + 1.0));

  Rab_rhs(-1, *p_jK, &dF, Fb, -ck);
  return dF;
}

//------------------------------------------------------------------------------
double AMk::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
    return 0.0;
  }

  const auto K = double(m_K);
  const auto sk = double(Fa.kappa() - Fb.kappa());
  if (sk == 0.0)
    return 0.0;
  assert(m_K != 0); // should already be discounted!
  const auto ck = sk / std::sqrt(K * (K + 1.0));

  return -ck * Rab(-1, *p_jK, Fa, Fb);
}

//------------------------------------------------------------------------------
void AMk::updateFrequency(const double omega) {
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

  Pab_rhs(-1, *p_jK, &dF, Fb);
  return dF;
}

//------------------------------------------------------------------------------
double Phi5k::radialIntegral(const DiracSpinor &Fa,
                             const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  return Pab(-1, *p_jK, Fa, Fb);
}

//------------------------------------------------------------------------------
void Phi5k::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    p_jK = &m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
    p_jK = &m_jK;
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

  Pab_rhs(+1, *p_jK, &dF, Fb);
  return dF;
}

//------------------------------------------------------------------------------
double S5k::radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  if (isZero(Fa.kappa(), Fb.kappa())) {
    return 0.0;
  }

  return Pab(+1, *p_jK, Fa, Fb);
}

//------------------------------------------------------------------------------
void S5k::updateFrequency(const double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);

  if (m_jl) {
    // nb: may not be exact! Ensure lookup table is dense enough!
    p_jK = &m_jl->jL_nearest(m_K, q);
  } else {
    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
    p_jK = &m_jK;
  }
}

//------------------------------------------------------------------------------
// Factory for relativistic multipole operators.
std::unique_ptr<DiracOperator::TensorOperator>
MultipoleOperator(const Grid &grid, int k, double omega, char type, char comp,
                  bool low_q, const SphericalBessel::JL_table *jl) {

  // normalise to upper case
  type = static_cast<char>(std::toupper(static_cast<unsigned char>(type)));
  comp = static_cast<char>(std::toupper(static_cast<unsigned char>(comp)));

  const bool Vector = (type == 'V');
  const bool AxialVector = (type == 'A');
  const bool Scalar = (type == 'S');
  const bool PseudoScalar = (type == 'P');

  const bool Electric = (comp == 'E');
  const bool Magnetic = (comp == 'M');
  const bool Longitudinal = (comp == 'L');
  // nb: don't confuse with transverse!!
  const bool Temporal = (comp == 'T' || comp == '0');

  // basic validation
  if (!(Vector || AxialVector || Scalar || PseudoScalar)) {
    std::cout << "make_tensor_op: invalid type '" << type
              << "' (expected V,A,S,P)\n";
    return nullptr;
  }

  // For V/A we require a valid component; for S/P we ignore comp but still validate if given.
  if ((Vector || AxialVector) &&
      !(Electric || Magnetic || Longitudinal || Temporal)) {
    std::cout << "make_tensor_op: invalid component '" << comp << "' for type "
              << type << " (expected T,E,M,L)\n";
    return nullptr;
  }
  if ((Scalar || PseudoScalar) &&
      !(Electric || Magnetic || Longitudinal || Temporal)) {
    // If user passed junk comp for S/P, complain (helps catch typos)
    std::cout << "make_tensor_op: invalid component '" << comp << "' for type "
              << type
              << " (S/P ignore component, but you gave an invalid one; "
                 "expected T,E,M,L)\n";
    return nullptr;
  }

  // -------- low-q branch --------
  if (low_q && Electric && Vector)
    return std::make_unique<VEk_lowq>(grid, k, omega);
  if (low_q && Electric && AxialVector)
    return std::make_unique<AEk_lowq>(grid, k, omega);

  if (low_q && Longitudinal && Vector)
    return std::make_unique<VLk_lowq>(grid, k, omega);
  if (low_q && Longitudinal && AxialVector)
    return std::make_unique<ALk_lowq>(grid, k, omega);

  if (low_q && Magnetic && Vector)
    return std::make_unique<VMk_lowq>(grid, k, omega);
  if (low_q && Magnetic && AxialVector)
    return std::make_unique<AMk_lowq>(grid, k, omega);

  if (low_q && Temporal && Vector)
    return std::make_unique<Phik_lowq>(grid, k, omega);
  if (low_q && Temporal && AxialVector)
    return std::make_unique<Phi5k_lowq>(grid, k, omega);

  if (low_q && Scalar)
    return std::make_unique<Sk_lowq>(grid, k, omega);
  if (low_q && PseudoScalar)
    return std::make_unique<S5k_lowq>(grid, k, omega);

  // -------- normal-q branch --------

  if (Electric && Vector)
    return std::make_unique<VEk>(grid, k, omega, jl);
  if (Electric && AxialVector)
    return std::make_unique<AEk>(grid, k, omega, jl);

  if (Longitudinal && Vector)
    return std::make_unique<VLk>(grid, k, omega, jl);
  if (Longitudinal && AxialVector)
    return std::make_unique<ALk>(grid, k, omega, jl);

  if (Magnetic && Vector)
    return std::make_unique<VMk>(grid, k, omega, jl);
  if (Magnetic && AxialVector)
    return std::make_unique<AMk>(grid, k, omega, jl);

  if (Temporal && Vector)
    return std::make_unique<Phik>(grid, k, omega, jl);
  if (Temporal && AxialVector)
    return std::make_unique<Phi5k>(grid, k, omega, jl);

  if (Scalar)
    return std::make_unique<Sk>(grid, k, omega, jl);
  if (PseudoScalar)
    return std::make_unique<S5k>(grid, k, omega, jl);

  // Should be unreachable if validation above is correct
  std::cout << "make_tensor_op: no operator matched (type=" << type
            << ", comp=" << comp << ", low_q=" << low_q << ")\n";
  return nullptr;
}

} // namespace DiracOperator