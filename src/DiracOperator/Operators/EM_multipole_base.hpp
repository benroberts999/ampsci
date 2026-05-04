#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include <memory>

namespace DiracOperator {

/*!
  @brief Intermediate abstract base class for all EM relativistic multipole
  operators.

  @details

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
            T^{(\sigma)}_{KQ},
  \end{align}  
  \f]

  with \f$ \tilde\gamma = 1,\f$ \f$ \gamma^5, \f$ \f$ \gamma^0, \f$ \f$ i\gamma^0\gamma^5 \f$
  for vector (V), axial(pseudo)-vector (A), scalar (S), or pseudoscalar (P).

  Stores the construction parameters needed to reconstruct the concrete
  operator via the MultipoleOperator factory, enabling type-erased cloning
  without the caller knowing the derived type.

  Specifically, the stored parameters are: a pointer to the radial grid, the
  current frequency @p omega (updated by each derived class's
  updateFrequency()), the Lorentz type character @p m_type, the component
  character @p m_comp, the low-q flag @p m_low_q, an optional Bessel table
  pointer @p m_jl (shallow-owned; shared between original and clone), and a
  form character @p m_form (used to distinguish the length-form electric
  operator VEk_Len from the velocity-form VEk).

  All concrete EM multipole operators -- VEk, VMk, VLk, Phik, Sk, AEk, ALk,
  AMk, Phi5k, S5k, and their low-q counterparts -- derive from this class.

  @note This class does not implement angularF() and so remains abstract.
  Operators must be constructed via the concrete derived types or the
  MultipoleOperator factory.
*/
class EM_multipole : public TensorOperator {
protected:
  /*!
    @brief Initialise the EM_multipole base layer.
    @details Called exclusively from derived-class constructors.

    @param rank_k   Tensor rank of the operator.
    @param pi       Parity.
    @param constant Overall multiplicative constant.
    @param vec      Radial vector (usually gr.r(), sometimes empty).
    @param RorI     Real or imaginary matrix elements.
    @param freq_dep True if the operator is frequency-dependent.
    @param grid     Pointer to the radial Grid (stored for clone()).
    @param type     Lorentz type: 'V', 'A', 'S', or 'P'.
    @param comp     Component: 'E', 'M', 'L', or 'T'.
    @param low_q    True for the low-momentum (long-wavelength) approximation.
    @param jl       Optional pointer to a precomputed Bessel table.
    @param form     'V' (velocity/V-form) or 'L' (length-form, VEk_Len only).
  */
  EM_multipole(int rank_k, Parity pi, double constant,
               const std::vector<double> &vec, Realness RorI, bool freq_dep,
               const Grid *grid, char type, char comp, bool low_q,
               const SphericalBessel::JL_table *jl = nullptr, char form = 'V')
    : TensorOperator(rank_k, pi, constant, vec, 0, RorI, freq_dep),
      m_grid(grid),
      m_type(type),
      m_comp(comp),
      m_low_q(low_q),
      m_form(form),
      m_jl(jl) {}

public:
  //! Returns the precomputed Bessel table pointer (may be nullptr).
  const SphericalBessel::JL_table *jl() const { return m_jl; }

  //! Returns a human-readable label, e.g. "T^E_1", "T^M5_2", "t_1", "P_1".
  std::string name() const override {
    const bool axial = (m_type == 'A');
    const std::string g5 = axial ? "5" : "";
    std::string base;
    if (m_type == 'V' || m_type == 'A') {
      switch (m_comp) {
      case 'E':
        base = std::string("T^E") + g5;
        break;
      case 'M':
        base = std::string("T^M") + g5;
        break;
      case 'L':
        base = std::string("T^L") + g5;
        break;
      case 'T':
        base = std::string("t") + g5;
        break;
      default:
        base = "?";
        break;
      }
    } else if (m_type == 'S') {
      base = "S";
    } else if (m_type == 'P') {
      base = "P";
    } else {
      base = "?";
    }
    if (m_form == 'L')
      base += "(Len)";
    if (m_low_q)
      base += "(lq)";
    return base + "_" + std::to_string(m_rank);
  }

  /*!
    @brief Angular factor linking the radial integral to the reduced matrix
    element: \f$ \langle a \| h \| b \rangle = F(\kappa_a,\kappa_b) \cdot
    \int \! dr \f$.
    @details
    Returns \f$ C^K(\kappa_a, \kappa_b) \f$ or
    \f$ C^K(\kappa_a, -\kappa_b) \f$ depending on the Lorentz structure.
    The flipped form applies to: vector-magnetic (VM), all axial components
    except axial-magnetic (AE, AL, AT), and pseudoscalar (P).

    @note \f$ C^K(\kappa_a,-\kappa_b) = C^K(-\kappa_a,\kappa_b) \f$.
    @note VMk_lowq overrides this with an extra \f$(κ_a+κ_b)\f$ prefactor.
  */
  double angularF(const int ka, const int kb) const override {
    const bool flip = (m_type == 'V' && m_comp == 'M') ||
                      (m_type == 'A' && m_comp != 'M') || (m_type == 'P');
    return flip ? Angular::Ck_kk(m_rank, ka, -kb) :
                  Angular::Ck_kk(m_rank, ka, kb);
  }

  /*!
    @brief Updates the tensor rank and adjusts parity accordingly.
    @details
    Parity follows the tensor rank: even K gives even parity for most
    operators, but is flipped (even K → odd parity) for vector-magnetic,
    all axial components except axial-magnetic, and pseudoscalar operators.
    The same flip rule as angularF() applies.

    @note Must be followed by a call to updateFrequency() to recompute
    the radial Bessel vectors for the new rank.
  */
  void updateRank(int new_K) override {
    m_rank = new_K;
    const bool flip = (m_type == 'V' && m_comp == 'M') ||
                      (m_type == 'A' && m_comp != 'M') || (m_type == 'P');
    m_parity = (Angular::evenQ(new_K) != flip) ? Parity::even : Parity::odd;
  }

  /*!
    @brief Creates a fully independent copy of this operator at its current
    (rank, frequency) state via the MultipoleOperator factory.
    @details
    The Bessel table pointer @p m_jl is shallow-copied (the clone shares the
    same external table), but all computed radial vectors are regenerated
    independently on the first updateFrequency() call.

    @return A std::unique_ptr<TensorOperator> to the cloned operator.
  */
  std::unique_ptr<TensorOperator> clone() const override;

protected:
  const Grid *m_grid{nullptr};
  double m_omega{
    0.0}; //!< Current frequency; cached by each derived updateFrequency().
  char m_type{};
  char m_comp{};
  bool m_low_q{};
  char m_form{};
  const SphericalBessel::JL_table *m_jl{nullptr};

public:
  EM_multipole(const EM_multipole &) = default;
  EM_multipole &operator=(const EM_multipole &) = default;
  EM_multipole(EM_multipole &&) = default;
  EM_multipole &operator=(EM_multipole &&) = default;
};

} // namespace DiracOperator
