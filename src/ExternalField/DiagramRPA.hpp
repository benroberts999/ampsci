#pragma once
#include "CorePolarisation.hpp"
#include "Coulomb/QkTable.hpp"
#include "HF/Breit.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace HF {
class HartreeFock;
}

namespace ExternalField {

//! RPA correction to matrix elements, using Diagram technique
/*! @details
  The basic equation is:

  \f[
      \langle{w}|{\delta V_{\pm}}|{v}\rangle &= 
      \sum_{na}\frac{t_{na}\,\widetilde g_{wavn}}{\varepsilon_a - \varepsilon_n \pm \omega}
      +
      \sum_{na}\frac{t_{an}\,\widetilde g_{wnva}}{\varepsilon_a - \varepsilon_n \mp \omega},
  \f]

  These should be solved self-consistently for all core-excited matrix elements.
  Should be solved independently for each frequency, including negative frequencies. To maintain Hermicity:

  \f[
    [\delta V_{\pm}(\omega)]^\dag = \delta V_{\mp}(-\omega).
  \f]

  The reduced matrix elements are:

  \f[
    \langle{w}\|{\delta V_{\pm}}\|{v}\rangle = 
    \sum_{na}\frac{(-1)^{k+w-n}}{[k]}\Bigg(
      \frac{T_{na}\,W^k_{wavn}}{\varepsilon_a - \varepsilon_n \pm \omega}
      +
      \frac{(-1)^{a-n}T_{an}\,W^k_{wnva}}{\varepsilon_a - \varepsilon_n \mp \omega}
    \Bigg).
  \f]

  @note
  Stores a pointer to the external field operator.
  That operator must remain alive for life of DiagramRPA object.
  Also, for frequency-dependent operators, updating the frequency of the operator externally will affect 
  the results from this class (as it should).
  The DiagramRPA::solve_core() function should be re called if the external operator is updated.

  @note
  Stores required W^k integrals (only those which are required).
  Writes them to disk by default.
  This can use significant memory, particularly for large basis.
  At the moment, always uses full basis -- this may optionally change in the future.

  @note
  Calling DiagramRPA::dV() without calling DiagramRPA::solve_core() will return 0.
  This is different from previous behaviour, where it returned dV^1, but in line with others.
  Call with max_its = 0 will also return 0.
  Call with max_its = 1 will return lowest-order dV correction.
*/
class DiagramRPA : public CorePolarisation {

private:
  const HF::HartreeFock *p_hf;
  std::optional<HF::Breit> m_Br{};
  std::vector<DiracSpinor> m_holes{};
  std::vector<DiracSpinor> m_excited{};

  // t0 (no RPA). These stay constant for frequency-independent operators
  // But change for f-dependent. Set during solve_core
  std::vector<std::vector<double>> m_t0am{};
  std::vector<std::vector<double>> m_t0ma{};
  // t's updated each solve_core itteration
  std::vector<std::vector<double>> m_tam{};
  std::vector<std::vector<double>> m_tma{};

  // Note: W's depend on rank (also parity)! Can re-use!?
  // These are probably an excellent candidate for unordered_map?
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wanmb{};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wabmn{};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wmnab{};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wmban{};

public:
  //! @brief Constructs DiagramRPA from a basis and a Hartree-Fock object
  /*! @details
    Splits @p basis into core (hole) and excited states, then computes the
    Coulomb W-matrix integrals required by the diagram RPA equations.
    If @p atom is non-empty the W matrices are read from a binary cache file
    (named from the atom label, operator rank, parity, and basis string) to
    avoid recalculation; if no matching file exists the matrices are computed
    and the file is written for future use.
    The Breit interaction is included automatically when present in @p in_hf.
    @param h      Pointer to the external-field tensor operator; must not be null
    @param basis  Complete single-particle basis (holes and excited states)
    @param in_hf  Pointer to the Hartree-Fock object; must not be null
    @param atom   Atom label used to form the W-matrix cache filename;
                  pass an empty string to disable file I/O
    @param print  If true, prints progress information to stdout
  */
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const std::vector<DiracSpinor> &basis,
             const HF::HartreeFock *in_hf, const std::string &atom = "",
             bool print = true);

  //! @brief Constructs DiagramRPA by reusing W matrices from an existing instance
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const DiagramRPA *const drpa);

public:
  //! @brief Iterates the RPA equations to convergence for the core electrons
  /*! @details
    Performs a self-consistent iterative solution of the diagram-method RPA
    equations at external-field frequency @p omega.
    The zeroth-order matrix elements \f$ t^{(0)}_{am} \f$ are calculated (so long as max_its>0).
    Each iteration updates the RPA matrix elements \f$ t_{am} \f$ and
    \f$ t_{ma} \f$ using the precomputed Coulomb W-matrix integrals until the
    maximum relative change falls below eps_target(), or @p max_its iterations
    have been performed.
    @param omega    External-field frequency in atomic units
    @param max_its  Maximum number of iterations (default 200)
    @param print    If true, prints per-iteration progress to stdout
    
    @note
    Calling DiagramRPA::dV() without calling DiagramRPA::solve_core() will return 0.
    Call with max_its = 0 will mean that dV() also returns 0.
    Call with max_its = 1 will mean that dV() returns the lowest-order dV correction.
  */
  void solve_core(double omega, int max_its = 200,
                  bool print = true) override final;

  //! @brief Returns the RPA method identifier (diagram)
  Method method() const final { return Method::diagram; }

  //! @brief Returns the RPA correction to a reduced matrix element
  /*! @details
    Computes the many-body RPA correction to the reduced matrix element
    \f[ \langle a \| \delta V \| b \rangle \f]
    using the converged RPA matrix elements \f$ t_{am} \f$ from the most
    recent solve_core() call.  Should be called after solve_core() has
    converged.
    @param Fa  Bra state
    @param Fb  Ket state
    @return    Reduced matrix element of the RPA correction

    @note
    Not Hermitian - solve_core() must be called with negaitve omega
  */
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const override final;

  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb, bool conj) const;

  //! @brief Calculates reduced right-hand-side, projected onto kappa: [dV|phi_m]_kappa
  DiracSpinor dV_rhs(int kappa, const DiracSpinor &Fm,
                     bool conj = false) const final;

  //! @brief Resets the RPA matrix elements to their unperturbed (zeroth-order) values
  /*! @details
    Clears the \f$ t_{am} \f$ and \f$ t_{ma} \f$ matrices (which encode the
    hole-excited RPA corrections) back to the state prior to any
    solve_core() call.
    Use this to discard a failed or poorly converged solution before restarting,
    for example with a different frequency or damping parameter.
  */
  void clear() final;

  //! @brief Updates the zeroth-order matrix elements and resets the RPA solution
  /*! @details
    Recomputes the lowest-order (no-RPA) matrix elements \f$ t^{(0)}_{am} \f$
    for operator @p h . Does not update the rpa T_am values.
    If @p h is null the currently stored operator is used unchanged.
    @param h  Optional replacement tensor operator; pass nullptr to keep the
              current operator
  */
  void update_t0s(const DiracOperator::TensorOperator *const h = nullptr);

private:
  // Note: only writes W (depends on k/pi, and basis). Do not write t's, since
  // they depend on operator. This makes it very fast when making small changes
  // to operator (don't need to re-calc W)
  // Note: doesn't depend on grid!
  bool read_write(const std::string &fname, IO::FRW::RoW rw, bool print = true);

  // Calculates all required W^k integrals (uses rank, parity of h)
  void fill_W_matrix(const DiracOperator::TensorOperator *const h, bool print);
  // Calculates all t_am matrix elements (without RPA)
  void setup_ts(const DiracOperator::TensorOperator *const h);

public:
  DiagramRPA &operator=(const DiagramRPA &) = delete;
  DiagramRPA(const DiagramRPA &) = default;
  ~DiagramRPA() = default;
};

} // namespace ExternalField
