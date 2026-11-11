#pragma once
#include "CorePolarisation.hpp"
#include "TDHF.hpp"
#include <string>
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace MBPT {
class CorrelationPotential;
}

namespace HF {
class HartreeFock;
class Breit;
} // namespace HF

namespace ExternalField {

/*!
  @brief Double core polarisation (DCP) by second-order TDHF: correction to
  a second-order amplitude of two external fields.

  @details
  For two fields, \f$ t \f$ (rank \f$ k_t \f$, frequency \f$ \omega_t \f$)
  and \f$ s \f$ (\f$ k_s \f$, \f$ \omega_s \f$), each core orbital
  \f$ \phi_a \f$ has the first-order corrections \f$ \varphi^a_{t\pm} \f$,
  \f$ \varphi^a_{s\pm} \f$ (see @ref TDHF), and second-order corrections at
  \f$ \pm(\omega_t+\omega_s) \f$:
  \f[
    (h_{\rm HF} - \en_a \mp \omega_t \mp \omega_s)\varphi^a_{\pm\pm}
      = -(\tilde t_\pm - \delta\en^a_{t\pm})\varphi^a_{s\pm}
        -(\tilde s_\pm - \delta\en^a_{s\pm})\varphi^a_{t\pm}
        -(\delta^2V_{\pm\pm} - \delta^2\en^a_{\pm\pm})\phi_a ,
  \f]
  with \f$ \tilde t_\pm = t_\pm + \delta V^t_\pm \f$. The correction to the
  Hartree-Fock potential, \f$ \delta^2V_{\pm\pm} \f$, is built from pairs of
  the first-order corrections and from the second-order corrections
  themselves, so the equations are solved self-consistently.

  The DCP correction to the amplitude is
  \f$ \matel{w}{\delta^2V_{++}}{v} \f$; coupled to total rank \f$ K \f$,
  this is \f$ \delta A^K \f$, returned by dV().

  The second-order corrections are stored in channels \f$ \beta \f$ of
  definite angular momentum: \f$ [\varphi^a_{++}]^K_\beta \f$ as dX, and
  \f$ [\varphi^a_{--}]^K_\beta \f$ as dY.

  \par Construction
  Requires the two first-order @ref TDHF objects (for \f$ t \f$ and
  \f$ s \f$), the total rank \f$ K \f$, and a const @ref HF::HartreeFock
  object. The first-order equations must be solved (TDHF::solve_core) before
  the second-order ones.
*/
class TDHF_DCP : public CorePolarisation {

protected:
  // Second-order corrections: dX = [phi^a_{++}]^K_beta (m_X), and
  // dY = [phi^a_{--}]^K_beta (m_Y). For each core orbital a, one for each
  // channel beta:
  // j(beta) = j(a)-K, ..., j(a)+K.  And: pi(beta) = pi(a)*pi(t)*pi(s)
  std::vector<std::vector<DiracSpinor>> m_X{};
  std::vector<std::vector<DiracSpinor>> m_Y{};
  std::vector<std::vector<DiracSpinor>> m_hFcore{};
  // can just write these to disk! Read them in, continue as per normal

  const HF::HartreeFock *const p_hf;
  const std::vector<DiracSpinor> m_core;
  // const std::vector<double> m_Hmag;
  const double m_alpha;
  const HF::Breit *const p_VBr;

  // First-order TDHF for each of the two fields: hold the operators t, s,
  // and the first-order corrections (Xt, Yt), (Xs, Ys)
  const ExternalField::TDHF *m_T;
  const ExternalField::TDHF *m_S;

public:
  /*!
    @brief Constructs second-order TDHF for the pair of fields t, s, coupled
    to total rank K.

    @details
    The second-order corrections have rank @p K, parity
    \f$ \pi_t\pi_s \f$, and matrix elements that are imaginary if those of
    exactly one of \f$ t \f$, \f$ s \f$ are.

    @param T   First-order TDHF for the field \f$ t \f$.
    @param S   First-order TDHF for the field \f$ s \f$.
    @param K   Total rank: \f$ |k_t-k_s| \le K \le k_t+k_s \f$.
    @param hf  @ref HF::HartreeFock object defining the core; the same as
               used for @p T and @p S.

    @note Stores the pointers @p T and @p S: those objects must remain
    alive for the lifetime of this one.
  */
  TDHF_DCP(const ExternalField::TDHF *const T,
           const ExternalField::TDHF *const S, int K,
           const HF::HartreeFock *const hf);

  //! Solves TDHF equations self-consistantly for core electrons at frequency
  //! omega.
  /*! @details Will iterate up to a maximum of max_its. Set max_its=1
     to get first-order correction [note: no damping is used for first
     itteration]. If print=true, will write progress to screen.
  */
  virtual void solve_core(const double omega, int max_its = 100,
                          const bool print = true) override;

  //! Returns RPA method
  virtual Method method() const override { return Method::TDHF; }

  //! Clears the dPsi orbitals (sets to zero)
  virtual void clear() override final;

  /*!
    @brief Reduced DCP amplitude between Fa and Fb, of total rank K.

    @details
    \f$ \delta A^K_{ab} \f$: the matrix element of \f$ \delta^2V_{++} \f$
    (at \f$ \omega_t+\omega_s \f$), coupled to rank \f$ K \f$. If @p conj,
    that of the \f$ \delta^2V_{--} \f$ component instead.

    @param Fa    Final state.
    @param Fb    Initial state.
    @param conj  False: the \f$ ++ \f$ component; true: the \f$ -- \f$ one.
  */
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb,
            bool conj = false) const override;

  //! Returns "reduced partial matrix element RHS": dV||Fb}.
  //! Note: Fa * dV_rhs(..) equiv to dV(..)
  DiracSpinor dV_rhs(const int kappa_n, const DiracSpinor &Fm,
                     bool conj = false) const override;

  //! Returns const ref to the second-order corrections (all channels) of
  //! core orbital Fc: dX (XorY = X) or dY (XorY = Y).
  const std::vector<DiracSpinor> &get_dPsis(const DiracSpinor &Fc,
                                            dPsiType XorY) const;
  //! Returns const ref to the second-order correction of core orbital Fc
  //! in the channel kappa_x.
  const DiracSpinor &get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                const int kappa_x) const;

  //! Forms \delta Psi_v for valence state Fv (including core pol.) - 1 kappa
  /*! @details 
   Solves
   \f[ (H + \Sigma - \epsilon - \omega)X = -(h + \delta V
   - \delta\epsilon)\psi \f]
   or
   \f[ (H + \Sigma - \epsilon + \omega)Y = -(h^\dagger + \delta V^\dagger
     - \delta\epsilon)Psi\f]
   Returns \f$ \chi_\beta \f$ for given kappa_beta, where
   \f[ X_{j,m} = (-1)^{j_\beta-m}tjs(j,k,j;-m,0,m)\chi_j \f]
   XorY takes values: dPsiType::X or dPsiType::Y.
   st takes values: StateType::ket or StateType::bra
  */
  DiracSpinor solve_dPsi(const DiracSpinor &Fa, const double omega,
                         dPsiType XorY, const int kappa_beta) const;
  //! Forms \delta Psi_v for valence state Fv for all kappas (see solve_dPsi)
  std::vector<DiracSpinor> solve_dPsis(const DiracSpinor &Fv,
                                       const double omega, dPsiType XorY) const;

  // //! Writes dPsi (f-component) to textfile
  // void print(const std::string &ofname = "dPsi.txt") const;

private:
  void initialise_dPsi();

  // Single iteration of TDHF equations
  std::pair<double, std::string> tdhf_core_it(double omega, double eta_damp);
  // Forms the set of sources of the second-order equations, for all core
  // orbitals and all channels: of the dX equations (XorY = X), or the dY
  // equations (XorY = Y). Second-order analogue of h*Fc.
  std::vector<std::vector<DiracSpinor>> form_hFcore(dPsiType XorY) const;
  // Source of the second-order equation for core orbital Fa, in the channel
  // kappa_beta: the first-order operators acting on the first-order
  // corrections
  DiracSpinor source_rhs(const DiracSpinor &Fa, dPsiType XorY,
                         const int kappa_beta) const;

public:
  TDHF_DCP &operator=(const TDHF_DCP &) = delete;
  TDHF_DCP(const TDHF_DCP &) = default;
  ~TDHF_DCP() = default;
};

} // namespace ExternalField
