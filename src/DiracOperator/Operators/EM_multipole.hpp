#pragma once
#include "DiracOperator/Operators/EM_multipole_base.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

// include the 'low qr' form here, simply for convenience
// (EM_multipole_lowqr.hpp already includes EM_multipole_base.hpp)
#include "EM_multipole_lowqr.hpp"

namespace DiracOperator {

//==============================================================================
//! @brief Vector electric multipole (transition) operator, Length-form: ( \f$ T^{(+1),{\rm Len}}_K \f$ ).
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
class VEk_Len final : public EM_multipole {
public:
  VEk_Len(const Grid &gr, int K, double omega,
          const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                   gr.r(), Realness::imaginary, true, &gr, 'V', 'E', false, jl,
                   'L') {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  VEk_Len(const VEk_Len &other)
    : EM_multipole(other),
      m_jK(other.m_jK),
      m_jKp1(other.m_jKp1),
      p_jK(other.p_jK == &other.m_jK ? &m_jK : other.p_jK),
      p_jKp1(other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1) {}
  VEk_Len &operator=(const VEk_Len &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK = other.m_jK;
      m_jKp1 = other.m_jKp1;
      p_jK = other.p_jK == &other.m_jK ? &m_jK : other.p_jK;
      p_jKp1 = other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Vector electric multipole (V-form) operator: \f$ V^E_K = T^{(+1)}_K(q) \f$
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
class VEk final : public EM_multipole {
public:
  VEk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                   gr.r(), Realness::imaginary, true, &gr, 'V', 'E', false,
                   jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK_on_qr{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK_on_qr{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  VEk(const VEk &other)
    : EM_multipole(other),
      m_jK_on_qr(other.m_jK_on_qr),
      m_jKp1(other.m_jKp1),
      p_jK_on_qr(other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr :
                                                         other.p_jK_on_qr),
      p_jKp1(other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1) {}
  VEk &operator=(const VEk &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK_on_qr = other.m_jK_on_qr;
      m_jKp1 = other.m_jKp1;
      p_jK_on_qr =
        other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr : other.p_jK_on_qr;
      p_jKp1 = other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Vector longitudinal multipole operator (V-form): \f$ V^L_K = T^{(-1)}_K(q) \f$
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
class VLk final : public EM_multipole {
public:
  VLk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                   gr.r(), Realness::imaginary, true, &gr, 'V', 'L', false,
                   jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK_on_qr{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK_on_qr{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  VLk(const VLk &other)
    : EM_multipole(other),
      m_jK_on_qr(other.m_jK_on_qr),
      m_jKp1(other.m_jKp1),
      p_jK_on_qr(other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr :
                                                         other.p_jK_on_qr),
      p_jKp1(other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1) {}
  VLk &operator=(const VLk &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK_on_qr = other.m_jK_on_qr;
      m_jKp1 = other.m_jKp1;
      p_jK_on_qr =
        other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr : other.p_jK_on_qr;
      p_jKp1 = other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Vector magnetic multipole operator: \f$ V^M_K = T^{(0)}_K(q) \f$
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
class VMk final : public EM_multipole {
public:
  VMk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                   gr.r(), Realness::imaginary, true, &gr, 'V', 'M', false,
                   jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  VMk(const VMk &other)
    : EM_multipole(other),
      m_jK(other.m_jK),
      p_jK(other.p_jK == &other.m_jK ? &m_jK : other.p_jK) {}
  VMk &operator=(const VMk &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK = other.m_jK;
      p_jK = other.p_jK == &other.m_jK ? &m_jK : other.p_jK;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Temporal component of the vector multipole operator: \f$ \Phi_K = t^K(q) \f$
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
class Phik final : public EM_multipole {
public:
  Phik(const Grid &gr, int K, double omega,
       const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                   gr.r(), Realness::real, true, &gr, 'V', 'T', false, jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  Phik(const Phik &other)
    : EM_multipole(other),
      m_jK(other.m_jK),
      p_jK(other.p_jK == &other.m_jK ? &m_jK : other.p_jK) {}
  Phik &operator=(const Phik &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK = other.m_jK;
      p_jK = other.p_jK == &other.m_jK ? &m_jK : other.p_jK;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Scalar multipole operator: \f$ S_K = t^K(q)\gamma^0 \f$
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
class Sk final : public EM_multipole {
public:
  Sk(const Grid &gr, int K, double omega,
     const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                   gr.r(), Realness::real, true, &gr, 'S', 'T', false, jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  Sk(const Sk &other)
    : EM_multipole(other),
      m_jK(other.m_jK),
      p_jK(other.p_jK == &other.m_jK ? &m_jK : other.p_jK) {}
  Sk &operator=(const Sk &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK = other.m_jK;
      p_jK = other.p_jK == &other.m_jK ? &m_jK : other.p_jK;
    }
    return *this;
  }
};

//==============================================================================
//==============================================================================
// Gamma^5 versions!

//==============================================================================
//! @brief Axial electric multipole operator: \f$ A^E_K = T^{(+1)}_K(q)\gamma^5 \f$
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
class AEk final : public EM_multipole {
public:
  AEk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                   gr.r(), Realness::real, false, &gr, 'A', 'E', false, jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK_on_qr{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK_on_qr{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  AEk(const AEk &other)
    : EM_multipole(other),
      m_jK_on_qr(other.m_jK_on_qr),
      m_jKp1(other.m_jKp1),
      p_jK_on_qr(other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr :
                                                         other.p_jK_on_qr),
      p_jKp1(other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1) {}
  AEk &operator=(const AEk &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK_on_qr = other.m_jK_on_qr;
      m_jKp1 = other.m_jKp1;
      p_jK_on_qr =
        other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr : other.p_jK_on_qr;
      p_jKp1 = other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Axial longitudinal multipole operator: \f$ A^L_K = T^{(-1)}_K(q)\gamma^5 \f$
/*!
  @details
  - Gamma^5 variant of the longitudinal multipole operator. Works like
    `VLk` but with the gamma^5 Dirac structure applied to the angular part.
  - Supports optional `const SphericalBessel::JL_table *jl` for lookup-table
    acceleration; otherwise uses on-the-fly Bessel evaluation.
*/
class ALk final : public EM_multipole {
public:
  ALk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                   gr.r(), Realness::real, true, &gr, 'A', 'L', false, jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK_on_qr{};
  std::vector<double> m_jKp1{};
  const std::vector<double> *p_jK_on_qr{nullptr};
  const std::vector<double> *p_jKp1{nullptr};

public:
  ALk(const ALk &other)
    : EM_multipole(other),
      m_jK_on_qr(other.m_jK_on_qr),
      m_jKp1(other.m_jKp1),
      p_jK_on_qr(other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr :
                                                         other.p_jK_on_qr),
      p_jKp1(other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1) {}
  ALk &operator=(const ALk &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK_on_qr = other.m_jK_on_qr;
      m_jKp1 = other.m_jKp1;
      p_jK_on_qr =
        other.p_jK_on_qr == &other.m_jK_on_qr ? &m_jK_on_qr : other.p_jK_on_qr;
      p_jKp1 = other.p_jKp1 == &other.m_jKp1 ? &m_jKp1 : other.p_jKp1;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Axial magnetic multipole operator: \f$ A^M_K = T^{(0)}_K(q)\gamma^5 \f$
/*!
  @details
  - Gamma^5 variant of the magnetic multipole operator. Analogous to
    `VMk` but with the gamma^5 Dirac structure applied where appropriate.
  - Uses spherical Bessel functions j_L(q*r) for radial dependence and
    accepts an optional `const SphericalBessel::JL_table *jl`.
*/
class AMk final : public EM_multipole {
public:
  AMk(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                   gr.r(), Realness::real, true, &gr, 'A', 'M', false, jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  AMk(const AMk &other)
    : EM_multipole(other),
      m_jK(other.m_jK),
      p_jK(other.p_jK == &other.m_jK ? &m_jK : other.p_jK) {}
  AMk &operator=(const AMk &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK = other.m_jK;
      p_jK = other.p_jK == &other.m_jK ? &m_jK : other.p_jK;
    }
    return *this;
  }
};

//==============================================================================
//! @brief Temporal component of the axial vector multipole operator
/*!
  @details
  \f$ \Theta_K = \Phi^5_K = t^K(q)\gamma^5 \f$

  - Gamma^5 variant of the temporal (time-like) component of the vector
    multipole operator. Functions like `Phik` but with gamma^5 applied to
    the appropriate spin-angular structure.
  - Radial dependence uses j_L(q*r) and the constructor accepts an
    optional `const SphericalBessel::JL_table *jl`.
*/
class Phi5k final : public EM_multipole {
public:
  Phi5k(const Grid &gr, int K, double omega,
        const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                   gr.r(), Realness::real, true, &gr, 'A', 'T', false, jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  Phi5k(const Phi5k &other)
    : EM_multipole(other),
      m_jK(other.m_jK),
      p_jK(other.p_jK == &other.m_jK ? &m_jK : other.p_jK) {}
  Phi5k &operator=(const Phi5k &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK = other.m_jK;
      p_jK = other.p_jK == &other.m_jK ? &m_jK : other.p_jK;
    }
    return *this;
  }
};

//==============================================================================
/*!
  @brief Pseudoscalar multipole operator: t^k (i g^0 g^5)
  @details

  \f[ P_K = S^5_K = t^K(q)(i\gamma^0\gamma^5) \f]

  - Implements the pseudoscalar multipole operator ~ \f$ e^{i q r} i \gamma^0 \gamma^5. \f$
  - Radial dependence is provided via spherical Bessel functions \f$ j_L(q*r). \f$
  - Supports an optional `const SphericalBessel::JL_table *jl` for precomputed
    Bessel lookup; otherwise computes on demand.
*/
class S5k final : public EM_multipole {
public:
  S5k(const Grid &gr, int K, double omega,
      const SphericalBessel::JL_table *jl = nullptr)
    : EM_multipole(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                   gr.r(), Realness::real, true, &gr, 'P', 'T', false, jl) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final;

private:
  std::vector<double> m_jK{};
  const std::vector<double> *p_jK{nullptr};

public:
  S5k(const S5k &other)
    : EM_multipole(other),
      m_jK(other.m_jK),
      p_jK(other.p_jK == &other.m_jK ? &m_jK : other.p_jK) {}
  S5k &operator=(const S5k &other) {
    if (this != &other) {
      EM_multipole::operator=(other);
      m_jK = other.m_jK;
      p_jK = other.p_jK == &other.m_jK ? &m_jK : other.p_jK;
    }
    return *this;
  }
};

//==============================================================================
//==============================================================================

//! Helper functions for the multipole operators
namespace multipole {

//! Convert from "transition form" to "moment form"
inline double moment_factor(int K, double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);
  return qip::double_factorial(2 * K + 1) / qip::pow(q, K) *
         std::sqrt(K / (K + 1.0));
}

} // namespace multipole

//==============================================================================
//==============================================================================

inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Multipole(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({
    {"", "Note: This function cannot use the Spherical Bessel looup table. If "
         "require efficiency for large number of q values, construct directly"},
    {"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
    {"omega", "Frequency: nb: q := alpha*omega [1.0e-4]"},
    {"type", "V,A,S,P (Vector, Axial, Scalar, Pseudoscalar) [V]"},
    {"low_q", "bool. Use low-q formulas (K=0 and 1 only, no L-form) [false]"},
    {"component", "E,M,L,T (electric, magnetic, longitudanel, temporal). "
                  "Temporal is forced if type = S or P. [E]"},
    {"form", "L,V (Length, Velocity); only for electric vector [L]"},
  });
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);

  const auto low_q = input.get("low_q", false);

  using namespace std::string_literals;
  const auto type = input.get("type", "V"s);
  const auto component = input.get("component", "E"s);
  const auto form = input.get("form", "V"s);

  const bool Vector = qip::ci_wc_compare(type, "V*");
  const bool AxialVector = qip::ci_wc_compare(type, "A*");
  const bool Scalar = qip::ci_wc_compare(type, "S*");
  const bool PseudoScalar = qip::ci_wc_compare(type, "P*");

  const bool Electric = qip::ci_wc_compare(component, "E*");
  const bool Magnetic = qip::ci_wc_compare(component, "M*");
  const bool Longitudinal = qip::ci_wc_compare(component, "L*");
  const bool Temporal = qip::ci_wc_compare(component, "T*");

  const bool LengthForm = qip::ci_wc_compare(form, "L*");

  if (LengthForm && !(Electric && Vector)) {
    std::cout << "Fail; Length form only valid for Electric Vector\n";
  }

  if (low_q) {
    if (Electric && Vector)
      return std::make_unique<VEk_lowq>(wf.grid(), k, omega);
    if (Electric && AxialVector)
      return std::make_unique<AEk_lowq>(wf.grid(), k, omega);

    // Longitudinal
    if (Longitudinal && Vector)
      return std::make_unique<VLk_lowq>(wf.grid(), k, omega);
    if (Longitudinal && AxialVector)
      return std::make_unique<ALk_lowq>(wf.grid(), k, omega);

    // Magnetic
    if (Magnetic && Vector)
      return std::make_unique<VMk_lowq>(wf.grid(), k, omega);
    if (Magnetic && AxialVector)
      return std::make_unique<AMk_lowq>(wf.grid(), k, omega);

    // Temporal
    if (Temporal && Vector)
      return std::make_unique<Phik_lowq>(wf.grid(), k, omega);
    if (Temporal && AxialVector)
      return std::make_unique<Phi5k_lowq>(wf.grid(), k, omega);

    if (Scalar)
      return std::make_unique<Sk_lowq>(wf.grid(), k, omega);
    if (PseudoScalar)
      return std::make_unique<S5k_lowq>(wf.grid(), k, omega);
  }

  // Electric:
  if (Electric && LengthForm && Vector)
    return std::make_unique<VEk_Len>(wf.grid(), k, omega);
  if (Electric && Vector)
    return std::make_unique<VEk>(wf.grid(), k, omega);
  if (Electric && AxialVector)
    return std::make_unique<AEk>(wf.grid(), k, omega);

  // Longitudinal
  if (Longitudinal && Vector)
    return std::make_unique<VLk>(wf.grid(), k, omega);
  if (Longitudinal && AxialVector)
    return std::make_unique<ALk>(wf.grid(), k, omega);

  // Magnetic
  if (Magnetic && Vector)
    return std::make_unique<VMk>(wf.grid(), k, omega);
  if (Magnetic && AxialVector)
    return std::make_unique<AMk>(wf.grid(), k, omega);

  // Temporal
  if (Temporal && Vector)
    return std::make_unique<Phik>(wf.grid(), k, omega);
  if (Temporal && AxialVector)
    return std::make_unique<Phi5k>(wf.grid(), k, omega);

  if (Scalar)
    return std::make_unique<Sk>(wf.grid(), k, omega);
  if (PseudoScalar)
    return std::make_unique<S5k>(wf.grid(), k, omega);

  std::cout << "Fail; Invalid Combination\n";
  return std::make_unique<NullOperator>();
}

//------------------------------------------------------------------------------
/*!
  @brief Factory for relativistic multipole operators.
  @details 
  Constructs and returns a specific multipole operator derived from
  DiracOperator::TensorOperator, based on the requested Lorentz structure,
  multipole component, and momentum-transfer regime.

  These are the \f$ t^K_Q \tilde\gamma \f$, \f$ T^{(\sigma)}_{KQ} \tilde\gamma \f$
  operators from the vector expansion:

  \f[
  \begin{align}
    e^{i\vec{q}\cdot\vec{r}} 
      &= \sqrt{4\pi}\sum_{KQ}\sqrt{[K]} \, 
        i^K \, {Y^*_{KQ}}{(\hat q)} \, t^K_Q(q,r),\\
    \vec{\alpha} \, e^{i\vec{q}\cdot\vec{r}}
      & = \sqrt{4\pi} \sum_{KQ\sigma} \sqrt{[K]} \, i^{K-\sigma} \, 
            \vec{Y}_{KQ}^{(\sigma)*}(\hat{{q}}) \, 
            T^{(\sigma)}_{KQ}.
  \end{align}  
  \f]

  The operator corresponds to a spherical multipole of rank @p k with
  frequency/energy transfer @p omega. The type and component determine
  the Lorentz structure and spatial character of the interaction.

  These are "frequency"-dependent operators, via momentum transfer, q.
  For electromagnetic interactions, \f$ \omega = q c\f$.
  The "updateFrequency()" function expects these units, so even when momentum
  transfer is not equal to energy exchange, we should pass \f$ qc \f$ to this function.

  @param grid   Radial grid on which the operator acts.
  @param k      Multipole rank (total angular momentum of the operator).
  @param omega  Energy (frequency) transfer.
  @param type   Interaction type:
                - 'V' : Vector
                - 'A' : Axial-vector
                - 'S' : Scalar
                - 'P' : Pseudoscalar
  @param comp   Multipole component (for V/A types):
                - 'E' : Electric
                - 'M' : Magnetic
                - 'L' : Longitudinal
                - 'T' : Temporal [caution - NOT transverse!]
                (Ignored for scalar and pseudoscalar operators.)
  @param low_q  If true, construct the low-momentum (long-wavelength)
                approximation of the operator.
  @param jl     Optional pointer to a precomputed spherical Bessel table.
                If provided, radial Bessel functions are taken from this
                table to avoid recomputation. If nullptr, they are generated
                internally as needed.

  @return A std::unique_ptr to the requested TensorOperator.

  @note Length-form electric operators (if implemented separately) are
        only valid for vector-electric ('V','E') combinations.

  @warning Invalid combinations of @p type and @p comp may result in
            a nullptr being returned.
*/
std::unique_ptr<DiracOperator::TensorOperator>
MultipoleOperator(const Grid &grid, int k, double omega, char type, char comp,
                  bool low_q, const SphericalBessel::JL_table *jl = nullptr);

} // namespace DiracOperator
