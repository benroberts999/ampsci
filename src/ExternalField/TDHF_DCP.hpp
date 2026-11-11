#pragma once
#include "CorePolarisation.hpp"
#include "TDHF.hpp"
#include <string>
#include <utility>
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

  /*!
    @brief Returns [d2V Fv]_kappa: the DCP correction to the potential
    acting on Fv, projected onto kappa_w.

    @details
    \f$ [\delta^2V_{++}\phi_v]^K_w \f$, such that Fw * dV_rhs(kappa_w, Fv)
    is \f$ \delta A^K_{wv} \f$. Sum of d2V_11_rhs(), d2V_Norm_rhs(), and
    d2V_2_rhs():
    \f[
      \delta A^K = \delta A^K(C1_{ts}) + \delta A^K(C1_{st})
        + \delta A^K(C2_{ts}) + \delta A^K(C2_{st}) + \delta A^K_{RL} .
    \f]

    @param kappa_w  Kappa of the projection (final state).
    @param Fv       State acted on (initial state).
    @param conj     False: the \f$ ++ \f$ component; true: the \f$ -- \f$
                    one.

    @note The \f$ -- \f$ component is the Hermitian conjugate:
    \f[
      \redmatel{w}{\delta^2V_{--}}{v}
        = \iota_t\iota_s\,(-1)^{j_w-j_v}\,\redmatel{v}{\delta^2V_{++}}{w} .
    \f]
    With @p conj = true, the result carries the factor
    \f$ \iota_t\iota_s \f$, as dY does.
  */
  DiracSpinor dV_rhs(const int kappa_w, const DiracSpinor &Fv,
                     bool conj = false) const override;

  /*!
    @brief Part of dV_rhs() from pairs of first-order corrections, one from
    each operator.

    @details
    The \f$ \delta^2V_{(1,1)} \f$ part; its matrix element is
    \f$ \delta A^K(C1_{ts}) + \delta A^K(C1_{st}) \f$:
    \f[
      \delta A^K(C1_{ts}) = \iota_s\sum_{a\alpha\beta}
        \frac{(-1)^{j_w+j_a}}{\sqrt{[K]}}
        \sixj{k_t}{k_s}{K}{j_\beta}{j_\alpha}{j_a}\,
        W^K_{w\beta_{s-}v\alpha_{t+}} ,
    \f]
    \f[
      \delta A^K(C1_{st}) = \iota_t\sum_{a\alpha\beta}
        \frac{(-1)^{K+k_s+k_t+j_w+j_a}}{\sqrt{[K]}}
        \sixj{k_t}{k_s}{K}{j_\alpha}{j_\beta}{j_a}\,
        W^K_{w\beta_{t-}v\alpha_{s+}} .
    \f]
    Does not depend on the second-order corrections.

    @param kappa_w  Kappa of the projection (final state).
    @param Fv       State acted on (initial state).
    @param conj     False: the \f$ ++ \f$ component; true: the \f$ -- \f$
                    one (the Hermitian conjugate: the first-order X and Y
                    swap, with the phase \f$ (-1)^{K+k_t+k_s} \f$).

    @note The Breit interaction is not included.
  */
  DiracSpinor d2V_11_rhs(const int kappa_w, const DiracSpinor &Fv,
                         bool conj = false) const;

  /*!
    @brief Part of dV_rhs() from the normalisation of the second-order
    corrections.

    @details
    The second-order corrections are solved orthogonal to their core
    orbital; the component along it is fixed by the normalisation,
    through the overlaps of the first-order corrections. The matrix
    element is \f$ \delta A^K(C2_{ts}) + \delta A^K(C2_{st}) \f$:
    \f[
      \delta A^K(C2_{ts}) = \iota_s\sum_{a\nu}
        \frac{(-1)^{K+k_s+k_t+j_w-j_\nu}}{\sqrt{[K]}}
        \sixj{k_t}{k_s}{K}{j_a}{j_a}{j_\nu}\,W^K_{wava}\,
        \braket{[\varphi^a_{s-}]_\nu}{[\varphi^a_{t+}]_\nu} ,
    \f]
    \f[
      \delta A^K(C2_{st}) = \iota_t\sum_{a\nu}
        \frac{(-1)^{j_w-j_\nu}}{\sqrt{[K]}}
        \sixj{k_t}{k_s}{K}{j_a}{j_a}{j_\nu}\,W^K_{wava}\,
        \braket{[\varphi^a_{t-}]_\nu}{[\varphi^a_{s+}]_\nu} .
    \f]
    Does not depend on the second-order corrections. The same for the
    \f$ ++ \f$ and \f$ -- \f$ components.

    @param kappa_w  Kappa of the projection (final state).
    @param Fv       State acted on (initial state).

    @note The Breit interaction is not included.
  */
  DiracSpinor d2V_Norm_rhs(const int kappa_w, const DiracSpinor &Fv) const;

  /*!
    @brief Part of dV_rhs() from the second-order corrections dX, dY.

    @details
    The \f$ \delta^2V_{(2)} \f$ part: the first-order formula for
    \f$ \delta V \f$, at rank \f$ K \f$, with the second-order corrections
    in place of the first-order ones. Its matrix element is
    \f[
      \delta A^K_{RL} = \sum_{a\beta}\frac{(-1)^{K+j_w-j_\beta}}{[K]}
        \Big(W^K_{wav\beta_{++}} + \iota_t\iota_s\,W^K_{w\beta_{--}va}\Big) .
    \f]
    This is the part that is iterated.

    @param kappa_w  Kappa of the projection (final state).
    @param Fv       State acted on (initial state).
    @param conj     False: the \f$ ++ \f$ component; true: the \f$ -- \f$
                    one (dX and dY swap).
  */
  DiracSpinor d2V_2_rhs(const int kappa_w, const DiracSpinor &Fv,
                        bool conj = false) const;

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
  /*!
    @brief Sets up the second-order corrections dX, dY (as zero): one for
    each core orbital and each channel allowed by the selection rules.
  */
  void initialise_dPsi();

  /*!
    @brief Single iteration of the TDHF equations.

    @param omega     Frequency.
    @param eta_damp  Damping factor.

    @return {eps, label of the worst channel}.
  */
  std::pair<double, std::string> tdhf_core_it(double omega, double eta_damp);

  /*!
    @brief Forms the sources of the second-order equations, for all core
    orbitals and all channels.

    @details
    Second-order analogue of the set of \f$ [t\phi_a]_\beta \f$ of TDHF:
    the part of the right-hand side that is fixed by the first-order
    solutions. See source_rhs().

    @param XorY  X: sources of the dX (\f$ ++ \f$) equations;
                 Y: of the dY (\f$ -- \f$) equations.

    @return Sources, indexed as dX and dY: [core orbital][channel].
  */
  std::vector<std::vector<DiracSpinor>> form_hFcore(dPsiType XorY) const;

  /*!
    @brief Source of the second-order equation for one core orbital and
    channel, from the first-order operators and corrections.

    @details
    For the \f$ ++ \f$ equation (@p XorY = X), returns
    \f$ [\sigma^a_{++}]^K_\beta - [\lambda^a_{++}]^K_\beta \f$, where
    \f[
      [\sigma^a_{++}]^K_\beta = \sum_\alpha\Big(
        c^1_{\beta\alpha a}\,[\tilde t_+[\varphi^a_{s+}]_\alpha]_\beta
        + c^2_{\beta\alpha a}\,[\tilde s_+[\varphi^a_{t+}]_\alpha]_\beta
      \Big) ,
    \f]
    \f[
      [\lambda^a_{++}]^K_\beta =
        c^2_{\beta aa}\,\redmatel{a}{\tilde t_+}{a}\,[\varphi^a_{s+}]_\beta
        + c^1_{\beta aa}\,\redmatel{a}{\tilde s_+}{a}\,[\varphi^a_{t+}]_\beta .
    \f]
    For the \f$ -- \f$ equation (@p XorY = Y), the same with the \f$ - \f$
    components of the operators and corrections, and
    \f$ c^{1,2}_{\beta\alpha a} \to c^{2,1}_{a\alpha\beta} \f$.
    The \f$ \delta^2V \f$ term is not included.

    @param Fa          Core orbital \f$ \phi_a \f$.
    @param XorY        X: the \f$ ++ \f$ equation; Y: the \f$ -- \f$ one.
    @param kappa_beta  Kappa of the channel \f$ \beta \f$.

    @return The source, a spinor of kappa @p kappa_beta.

    @note The first-order corrections Xt, Xs are taken to be the \f$ + \f$
    components of \f$ t \f$ and \f$ s \f$: \f$ \omega_t, \omega_s \ge 0 \f$.
  */
  DiracSpinor source_rhs(const DiracSpinor &Fa, dPsiType XorY,
                         const int kappa_beta) const;

  /*!
    @brief Angular coefficients of the terms of source_rhs() that pass
    through an intermediate state of angular momentum j_n.

    @details
    In terms of the coefficients
    \f[
      c^1_{wnv} = (-1)^{K+j_w+j_v}\sqrt{[K]}
        \sixj{K}{k_s}{k_t}{j_n}{j_w}{j_v} ,
      \qquad
      c^2_{wnv} = (-1)^{k_t+k_s+j_w+j_v}\sqrt{[K]}
        \sixj{K}{k_t}{k_s}{j_n}{j_w}{j_v} :
    \f]
    - \f$ ++ \f$ equation:
      \f$ \{c^1_{\beta na},\, c^2_{\beta na}\} \f$
    - \f$ -- \f$ equation:
      \f$ \{c^2_{an\beta},\, c^1_{an\beta}\} \f$

    @param twoj_a     2j of the core orbital.
    @param twoj_n     2j of the intermediate state.
    @param twoj_beta  2j of the channel.
    @param conj       False: the \f$ ++ \f$ equation; true: the \f$ -- \f$ one.

    @return {c_ts, c_st}: c_ts for t acting on an s correction; c_st for
    s acting on a t correction.
  */
  std::pair<double, double> source_coefs(int twoj_a, int twoj_n, int twoj_beta,
                                         bool conj) const;

  /*!
    @brief First-order operator, with its core polarisation, acting on F:
    reduced projection onto kappa.

    @details
    \f$ [\tilde h F]_\kappa \f$, with
    \f$ \tilde h = h + \delta V^h \f$, for \f$ h = t \f$ or \f$ s \f$.

    @param h_tilde  First-order TDHF of \f$ h \f$: holds the operator and
                    its \f$ \delta V^h \f$.
    @param kappa    Kappa of the projection.
    @param F        Spinor acted on (of definite kappa).
    @param conj     If true, the emission components \f$ h_- \f$,
                    \f$ \delta V^h_- \f$, in the convention of the Y
                    corrections.
  */
  DiracSpinor h_tilde_rhs(const TDHF &h_tilde, const int kappa,
                          const DiracSpinor &F, bool conj) const;

  /*!
    @brief First-order energy shift of a core orbital, <a||h + dV^h||a>.

    @param h_tilde  First-order TDHF of \f$ h \f$: holds the operator and
                    its \f$ \delta V^h \f$.
    @param Fa       Core orbital \f$ \phi_a \f$.
    @param conj     If true, that of the emission components; see
                    h_tilde_rhs().
  */
  double h_tilde_de(const TDHF &h_tilde, const DiracSpinor &Fa,
                    bool conj) const;

public:
  TDHF_DCP &operator=(const TDHF_DCP &) = delete;
  TDHF_DCP(const TDHF_DCP &) = default;
  ~TDHF_DCP() = default;
};

} // namespace ExternalField
