#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>
class Wavefunction;
class Grid;
namespace HF {
class HartreeFock;
}

/*!
  @brief Stores set of continuum orbitals, and solves using Hartree-Fock method.

  @details
  Solves in the Hartree-Fock potential (solveContinuumHF), or in a H-like
  potential with effective charge (solveContinuumZeff, or analytically with
  solveContinuumZeffAnalytic). Orbitals are energy-normalised (see
  DiracODE::solveContinuum). Each solve appends to orbitals (see clear).
*/
class ContinuumOrbitals {

public:
  /*!
    @brief Takes existing Wavefunction object to construct.

    @details
    Stores the grid, alpha, and a pointer to its HF object (may be nullptr if
    HF not set; only the Zeff methods may then be used).

    @warning The HF pointer is stored: @p wf must outlive this object.
  */
  ContinuumOrbitals(const Wavefunction &wf);
  /*!
    @brief Construct with Hartree-Fock pointer.

    @details
    Pointer is stored: must outlive this object. Must not be nullptr.
  */
  ContinuumOrbitals(const HF::HartreeFock *hf);

  ContinuumOrbitals &operator=(const ContinuumOrbitals &) = default;
  ContinuumOrbitals(const ContinuumOrbitals &) = default;
  ~ContinuumOrbitals() = default;

  /*!
    @brief Solves continuum states with energy ec between min/max l.

    @details
    Solves for each kappa with min_l <= l <= max_l, in the HF local potential
    (vnuc + vdir + vrad). For non-local methods (HartreeFock, ApproxHF),
    exchange with the full core is then included iteratively (damped; at most
    50 iterations, until |dF|^2/|F|^2 < 1e-6; no warning if not converged).

    @param ec             Continuum energy (au, > 0)
    @param min_l          Minimum l
    @param max_l          Maximum l
    @param psi            Ionised (bound) state; optional (nullptr: none)
    @param force_rescale  Rescale large-r part of the local potential, so it
                          falls off no slower than -Zion/r:
                          V(r) -> min(V(r), -Zion/r), with
                          Zion = max(1, -r V(r)) at the last grid point
    @param subtract_self  Subtract the self-interaction of @p psi from the
                          direct potential, \f$ y^0_{\psi\psi}(r) \f$
                          (V^{N-1} potential of the ion)
    @param force_orthog   Orthogonalise to the core states of the same kappa,
                          and to @p psi (to @p psi also at each exchange
                          iteration)
    @param average_tail   Fill the (zeroed) unresolved high-r tail with a
                          local average (see DiracODE::averageTail); done
                          after exchange, before orthogonalisation

    @return 0 (always).

    @note @p force_rescale and @p subtract_self only apply for the
          HartreeFock, Hartree, and ApproxHF methods.
    @warning Requires a Hartree-Fock object: dereferences a null pointer if
             constructed from a Wavefunction with no HF.
  */
  int solveContinuumHF(double ec, int min_l, int max_l,
                       const DiracSpinor *psi = nullptr,
                       bool force_rescale = false, bool subtract_self = true,
                       bool force_orthog = true, bool average_tail = false);

  /*!
    @brief Solves cntm states using simple H-like potential with Z_eff.

    @details
    Potential is \f$ -Z_{\rm eff}/r \f$ (point-like). Usually,
    \f$ Z_{\rm eff} = n\sqrt{2 I_{njl}} \f$, with \f$ I_{njl} \f$ the
    ionisation energy (au) and n the principal quantum number of the initial
    state.

    @param ec            Continuum energy (au, > 0)
    @param min_l         Minimum l
    @param max_l         Maximum l
    @param Z_eff         Effective charge
    @param Fi            Ionised (bound) state; may be nullptr
    @param force_orthog  Orthogonalise to @p Fi (only; not the core)
    @param average_tail  Fill the (zeroed) unresolved high-r tail with a
                         local average (see DiracODE::averageTail)

    @return 0 (always).
  */
  int solveContinuumZeff(double ec, int min_l, int max_l, double Z_eff,
                         const DiracSpinor *Fi, bool force_orthog,
                         bool average_tail = false);

  /*!
    @brief As solveContinuumZeff, but uses exact analytic Dirac-Coulomb states.

    @details
    Uses exact analytic Dirac-Coulomb functions rather than solving the ODE
    numerically (see DiracSpinor::exactHlike_cntm). Evaluated pointwise over
    the entire grid, so there is no tail cut-off (cf. solveContinuumZeff).

    @return 0 (always).

    @warning Relativistic case requires FLINT (see DiracContinuum::available):
             otherwise prints a warning, and the orbitals are NaN (or an
             assert fails, in debug builds).
  */
  int solveContinuumZeffAnalytic(double ec, int min_l, int max_l, double Z_eff,
                                 const DiracSpinor *Fi, bool force_orthog);

  /*!
    @brief Checks orthogonality between cntm and core orbitals.

    @details
    Overlaps with core states of the same kappa; prints each if @p print.

    @return Worst eps: the largest-magnitude overlap (signed); 0 if no HF.
  */
  double check_orthog(bool print = true) const;

  //! Resets (deletes) all orbitals
  void clear();

  //! The continuum orbitals (each solve appends to this list)
  std::vector<DiracSpinor> orbitals{};

private:
  void IncludeExchange(DiracSpinor &F_cntm, const DiracSpinor *F_i,
                       bool force_orthog, const std::vector<double> &vc);

  std::shared_ptr<const Grid> p_rgrid;
  const HF::HartreeFock *p_hf;
  double m_alpha;
};
