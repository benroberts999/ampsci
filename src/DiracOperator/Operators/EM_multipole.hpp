#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

// include the 'low qr' form here, simply for convenience
#include "EM_multipole_lowqr.hpp"

namespace DiracOperator {

//==============================================================================
//! @brief Vector electric multipole (transition) operator, Length-form: ($T^{(+1),{\rm Len}}_K$).
/*!
  @details
  - nb: q = alpha*omega, so "omega" = c*q
  - Electric multipole operator in the L-form (transition form
    convention): t^E_L. The radial dependence is expressed in terms of
    spherical Bessel functions j_L(q*r) with q = alpha * omega.
  - The operator supports using either on-the-fly Bessel evaluation or a
    precomputed `SphericalBessel::JL_table` supplied via the constructor
    (pointer `jl`). When `jl` is non-null the table is used to look up the
    closest j_L(q) vector for performance.
  - Note: The q value for jL(qr) will be the _nearest_ to the requested q 
    - you should ensure the lookup table is close enough, or this can lead to errors.
*/
class VEk_Len final : public TensorOperator {
public:
  VEk_Len(const Grid &gr, int K, double omega,
          const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::imaginary, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("V^E(Len)_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;
  //                         {

  //   DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  //   dF.min_pt() = Fb.min_pt();
  //   dF.max_pt() = Fb.max_pt();
  //   if (isZero(kappa_a, Fb.kappa()) || m_K == 0) {
  //     dF.min_pt() = 0;
  //     dF.max_pt() = 0;
  //     return dF;
  //   }

  //   const auto K = double(m_K);
  //   const auto c1 = double(kappa_a - Fb.kappa()) / (K + 1.0);
  //   const auto cx = std::sqrt((K + 1.0) / K);
  //   Rab_rhs(+1, *p_jK, &dF, Fb, cx);
  //   Pab_rhs(+1, *p_jKp1, &dF, Fb, -c1 * cx);
  //   Pab_rhs(-1, *p_jKp1, &dF, Fb, -cx);
  //   return dF;
  // }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto cc = double(Fa.kappa() - Fb.kappa()) / (K + 1.0);
    const auto cx = std::sqrt((K + 1.0) / K);

    return cx * (Rab(+1, *p_jK, Fa, Fb) - (cc + 1.0) * Vab(*p_jKp1, Fa, Fb) +
                 (1.0 - cc) * Wab(*p_jKp1, Fa, Fb));
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
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

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};

  // This is to avoid copy if not required:
  std::vector<double> m_jK{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  // Default shallow copy semantics
  VEk_Len(const VEk_Len &) = default;
  VEk_Len &operator=(const VEk_Len &) = default;
};

//==============================================================================
//! @brief Vector electric multipole (V-form) operator: $V^E_K = T^{(+1)}_K(q)$
/*!
  @details
  - nb: q = alpha*omega, so "omega" = c*q
  - Vector (spatial) electric multipole operator used in the V-form
    (transition convention). Radial dependence uses spherical
    Bessel functions j_L(q*r) and derived combinations like j_L(q*r)/(q*r)
    where needed; q = alpha * omega.
  - The constructor takes an optional `const SphericalBessel::JL_table *jl` to
    enable lookup from a precomputed table for improved performance.
  - Note: The q value for jL(qr) will be the _nearest_ to the requested q 
    - you should ensure the lookup table is close enough, or this can lead to errors.
*/
class VEk final : public TensorOperator {
public:
  VEk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::imaginary, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("V^E_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

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

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
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

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};

  // This is to avoid copy if not required:
  std::vector<double> m_jK_on_qr{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK_on_qr{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  VEk(const VEk &) = default;
  VEk &operator=(const VEk &) = default;
};

//==============================================================================
//! @brief Vector longitudinal multipole operator (V-form): $V^L_K = T^{(-1)}_K(q)$
/*!
  @details
  - nb: q = alpha*omega, so "omega" = c*q
  - Implements the longitudinal (scalar-like) component of the vector
    multipole operator. Radial dependence uses spherical Bessel functions
    j_L(q*r) and, when required, the combination j_L(q*r)/(q*r).
  - The constructor takes an optional `const SphericalBessel::JL_table *jl` to
    enable lookup from a precomputed table for improved performance.
  - Note: The q value for jL(qr) will be the _nearest_ to the requested q 
    - you should ensure the lookup table is close enough, or this can lead to errors.
*/
class VLk final : public TensorOperator {
public:
  VLk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::imaginary, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("V^L_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto dk = double(Fa.kappa() - Fb.kappa());

    return -dk * Pab(+1, *p_jK_on_qr, Fa, Fb) +
           K * Pab(-1, *p_jK_on_qr, Fa, Fb) - Pab(-1, *p_jKp1, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
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

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};

  // This is to avoid copy if not required:
  std::vector<double> m_jK_on_qr{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK_on_qr{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  VLk(const VLk &) = default;
  VLk &operator=(const VLk &) = default;
};

//==============================================================================
//! @brief Vector magnetic multipole operator: $V^M_K = T^{(0)}_K(q)$
/*!
  @details
  - nb: q = alpha*omega, so "omega" = c*q
  - Implements the magnetic multipole operator. The radial dependence is
    expressed via spherical Bessel functions j_L(q*r) with q = alpha * omega.
  - The constructor takes an optional `const SphericalBessel::JL_table *jl` to
    enable lookup from a precomputed table for improved performance.
  - Note: The q value for jL(qr) will be the _nearest_ to the requested q 
    - you should ensure the lookup table is close enough, or this can lead to errors.
*/
class VMk final : public TensorOperator {
public:
  VMk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::imaginary, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }

  std::string name() const override final {
    return std::string("V^M_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, -ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

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

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);
    if (m_jl) {
      // nb: may not be exact! Ensure lookup table is dense enough!
      p_jK = &m_jl->jL_nearest(m_K, q);
    } else {
      SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
      p_jK = &m_jK;
    }
  }

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};

  // This is to avoid copy if not required:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  VMk(const VMk &) = default;
  VMk &operator=(const VMk &) = default;
};

//==============================================================================
//! @brief Temporal component of the vector multipole operator: $\Phi_K = t^K(q)$
/*!
  @details
  - nb: q = alpha*omega, so "omega" = c*q
  - Implements the time-like (temporal) component of the vector multipole
    operator with explicit frequency dependence via j_L(q*r).
  - The constructor takes an optional `const SphericalBessel::JL_table *jl` to
    enable lookup from a precomputed table for improved performance.
  - Note: The q value for jL(qr) will be the _nearest_ to the requested q 
    - you should ensure the lookup table is close enough, or this can lead to errors.
*/
class Phik final : public TensorOperator {
public:
  Phik(const Grid &gr, int K, double omega,
       const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("Phi_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Rab(+1, *p_jK, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);
    if (m_jl) {
      // nb: may not be exact! Ensure lookup table is dense enough!
      p_jK = &m_jl->jL_nearest(m_K, q);
    } else {
      SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
      p_jK = &m_jK;
    }
  }

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};

  // This is to avoid copy if not required:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  Phik(const Phik &) = default;
  Phik &operator=(const Phik &) = default;
};

//==============================================================================
//! @brief Scalar multipole operator: $S_K = t^K(q)\gamma^0$
/*!
  @details
  - nb: q = alpha*omega, so "omega" = c*q
  - Implements the scalar multipole operator whose radial factor is
    e^{i q r} (represented via spherical Bessel functions j_L(q*r)). The
    operator includes the gamma^0 Dirac matrix structure.
  - The constructor takes an optional `const SphericalBessel::JL_table *jl` to
    enable lookup from a precomputed table for improved performance.
  - Note: The q value for jL(qr) will be the _nearest_ to the requested q 
    - you should ensure the lookup table is close enough, or this can lead to errors.
*/
class Sk final : public TensorOperator {
public:
  Sk(const Grid &gr, int K, double omega,
     const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("t^S_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Rab(-1, *p_jK, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);
    if (m_jl) {
      // nb: may not be exact! Ensure lookup table is dense enough!
      p_jK = &m_jl->jL_nearest(m_K, q);
    } else {
      SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &m_jK);
      p_jK = &m_jK;
    }
  }

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};

  // This is to avoid copy if not required:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  Sk(const Sk &) = default;
  Sk &operator=(const Sk &) = default;
};

//==============================================================================
//==============================================================================
// Gamma^5 versions!

//==============================================================================
//! @brief Electric multipole operator (gamma5 variant), V-form, frequency-dependent.
/*!
  @details
  - Gamma^5 variant of the electric multipole (vector) operator. Functions
    analogously to `VEk` but with the gamma^5 Dirac structure applied to the
    angular part.
  - Radial functions use spherical Bessel combinations j_L(q*r) or
    j_L(q*r)/(q*r) with q = alpha * omega.
  - Accepts an optional `const SphericalBessel::JL_table *jl` to use precomputed
    Bessel vectors; otherwise computes them on demand.
*/
class E5k_w final : public TensorOperator {
public:
  E5k_w(const Grid &gr, int K, double omega,
        const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^E5_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

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

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
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

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};
  std::vector<double> jK_on_qr{};
  std::vector<double> jKp1{};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  E5k_w(const E5k_w &) = default;
  E5k_w &operator=(const E5k_w &) = default;
};

//==============================================================================
//! @brief Longitudinal multipole operator (gamma5 variant), V-form.
/*!
  @details
  - Gamma^5 variant of the longitudinal multipole operator. Works like
    `VLk` but with the gamma^5 Dirac structure applied to the angular part.
  - Supports optional `const SphericalBessel::JL_table *jl` for lookup-table
    acceleration; otherwise uses on-the-fly Bessel evaluation.
*/
class L5k_w final : public TensorOperator {
public:
  L5k_w(const Grid &gr, int K, double omega,
        const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^5L_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto dk = double(Fa.kappa() + Fb.kappa());

    return -dk * Rab(-1, jK_on_qr, Fa, Fb) + K * Rab(+1, jK_on_qr, Fa, Fb) -
           Rab(+1, jKp1, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
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

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};
  std::vector<double> jK_on_qr{};
  std::vector<double> jKp1{};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  L5k_w(const L5k_w &) = default;
  L5k_w &operator=(const L5k_w &) = default;
};

//==============================================================================
//! @brief Magnetic multipole operator (gamma5 variant), frequency-dependent.
/*!
  @details
  - Gamma^5 variant of the magnetic multipole operator. Analogous to
    `VMk` but with the gamma^5 Dirac structure applied where appropriate.
  - Uses spherical Bessel functions j_L(q*r) for radial dependence and
    accepts an optional `const SphericalBessel::JL_table *jl`.
*/
class M5k_w final : public TensorOperator {
public:
  M5k_w(const Grid &gr, int K, double omega,
        const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }

  std::string name() const override final {
    return std::string("t^5M_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

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

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    if (m_jl) {
      jK = m_jl->jL_nearest(m_K, q);
    } else {
      SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jK);
    }
  }

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};
  std::vector<double> jK{};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  M5k_w(const M5k_w &) = default;
  M5k_w &operator=(const M5k_w &) = default;
};

//==============================================================================
//! @brief Temporal component of vector multipole operator (gamma5 variant).
/*!
  @details
  - Gamma^5 variant of the temporal (time-like) component of the vector
    multipole operator. Functions like `Phik` but with gamma^5 applied to
    the appropriate spin-angular structure.
  - Radial dependence uses j_L(q*r) and the constructor accepts an
    optional `const SphericalBessel::JL_table *jl`.
*/
class Phi5k_w final : public TensorOperator {
public:
  Phi5k_w(const Grid &gr, int K, double omega,
          const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("t^5_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Pab(-1, jk, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    if (m_jl) {
      jk = m_jl->jL_nearest(m_K, q);
    } else {
      SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
    }
  }

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};
  std::vector<double> jk{};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  Phi5k_w(const Phi5k_w &) = default;
  Phi5k_w &operator=(const Phi5k_w &) = default;
};

//==============================================================================
//! @brief Pseudoscalar multipole operator (gamma5 scalar), frequency-dependent.
/*!
  @details
  - Implements the pseudoscalar multipole operator ~ e^{i q r} gamma^0 gamma^5.
  - Radial dependence is provided via spherical Bessel functions j_L(q*r).
  - Supports an optional `const SphericalBessel::JL_table *jl` for precomputed
    Bessel lookup; otherwise computes on demand.
*/
class S5k_w final : public TensorOperator {
public:
  S5k_w(const Grid &gr, int K, double omega,
        const SphericalBessel::JL_table *jl = nullptr)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K),
        m_jl(jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^5S_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

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

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Pab(+1, jk, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    if (m_jl) {
      jk = m_jl->jL_nearest(m_K, q);
    } else {
      SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
    }
  }

private:
  int m_K;
  const SphericalBessel::JL_table *m_jl{nullptr};
  std::vector<double> jk{};

public:
  // Default shallow copy semantics (pointer member is shallow-copied)
  S5k_w(const S5k_w &) = default;
  S5k_w &operator=(const S5k_w &) = default;
};

//==============================================================================
//==============================================================================

//! Helper functions for the multipole operators
namespace multipole {

//! Convert from "transition form" to "moment form"
inline double moment_factor(int K, double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);
  // -1 for electric, +1 for magnetic?
  return -1.0 * qip::double_factorial(2 * K + 1) / qip::pow(q, K) *
         std::sqrt(K / (K + 1.0));
}

//! Scalar(pseudoscalar) multipole operator. Note: q/omega not set
inline std::unique_ptr<DiracOperator::TensorOperator>
S_K(const Grid &grid, int k, bool gamma5 = false,
    const SphericalBessel::JL_table *jl = nullptr) {
  if (gamma5)
    return std::make_unique<S5k_w>(grid, k, 1.0e-4, jl);
  return std::make_unique<Sk>(grid, k, 1.0e-4, jl);
}

//! Temporal part of vector(pseudovector) multipole operator. Note: q/omega not set
inline std::unique_ptr<DiracOperator::TensorOperator>
Phi_K(const Grid &grid, int k, bool gamma5 = false,
      const SphericalBessel::JL_table *jl = nullptr) {
  if (gamma5)
    return std::make_unique<Phi5k_w>(grid, k, 1.0e-4, jl);
  return std::make_unique<Phik>(grid, k, 1.0e-4, jl);
}

//! Spatial part of vector(pseudovector) multipole operator. Note: q/omega not set
inline std::unique_ptr<DiracOperator::TensorOperator>
V_sigma_K(const Grid &grid, int sigma, int k, bool gamma5 = false,
          const SphericalBessel::JL_table *jl = nullptr) {

  switch (sigma) {

    // "Electric"
  case +1:
    if (gamma5)
      return std::make_unique<E5k_w>(grid, k, 1.0e-4, jl);
    else
      return std::make_unique<VEk>(grid, k, 1.0e-4, jl);

    // "Longitudanal"
  case -1:
    if (gamma5)
      return std::make_unique<L5k_w>(grid, k, 1.0e-4, jl);
    else
      return std::make_unique<VLk>(grid, k, 1.0e-4, jl);

    // "Magnetic"
  case 0:
    if (gamma5)
      return std::make_unique<M5k_w>(grid, k, 1.0e-4, jl);
    else
      return std::make_unique<VMk>(grid, k, 1.0e-4, jl);
  }

  assert(false && "make_op: sigma must be -1, 0, or +1");
  return {}; // never reached when assertions are enabled
}

} // namespace multipole

//==============================================================================
//==============================================================================

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek_w_L(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  return std::make_unique<VEk_Len>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<E5k_w>(wf.grid(), k, omega);
  return std::make_unique<VEk>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Mk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for M1, =2 for M2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<M5k_w>(wf.grid(), k, omega);
  return std::make_unique<VMk>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Lk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<L5k_w>(wf.grid(), k, omega);
  return std::make_unique<VLk>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Vk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<Phi5k_w>(wf.grid(), k, omega);
  return std::make_unique<Phik>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Sk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<S5k_w>(wf.grid(), k, omega);
  return std::make_unique<Sk>(wf.grid(), k, omega);
}

} // namespace DiracOperator
