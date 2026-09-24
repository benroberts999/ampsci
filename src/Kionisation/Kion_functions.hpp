#pragma once
#include "DiracOperator/Operators/jL.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
class DiracSpinor;
class Grid;
class Wavefunction;
namespace HF {
class HartreeFock;
}
namespace ExternalField {
class TDHFcntm;
}

//! Functions for atomic ionisation form factors
namespace Kion {

//------------------------------------------------------------------------------
//! DM-electron couplings
enum class Coupling { Vector, Scalar, AxialVector, PseudoScalar, Error };

//------------------------------------------------------------------------------
/*!
  @brief Format for output file. (All new code should use xyz; matrix kept 
  for legacy code)

  @details
  - xyz    : For easy 2D interpolation. List formatted with each row 'E q K(E,q)'
  - matrix : Outputs entire matrix in table form. E and q grids printed prior.
*/
enum class OutputFormat { matrix, xyz, Error };

//------------------------------------------------------------------------------
/*!
  @brief Units used in output file

  @details
  - Atomic   : [q] = [1/a_0], [E] = Hartree
  - Particle : [q] = eV, [E] = eV
  - Form factors are defined to be dimensionless
*/
enum class Units { Atomic, Particle, Error };

//------------------------------------------------------------------------------
/*!
  @brief Method used to solve bound/continuum states for form factors

  @details
  - HF           : Real (Hartree-Fock) bound and continuum states (standard
                   method).
  - Zeff         : H-like (Zeff) bound and continuum states, solved
                   numerically with DiracODE.
  - ZeffAnalytic : H-like (Zeff) bound and continuum states, using exact
                   analytic Dirac-Coulomb functions. Relativistic continuum
                   requires FLINT (see DiracContinuum::available).
  - RPA          : Hartree-Fock states, with core-polarisation (RPA)
                   corrections to every amplitude, from the outgoing-wave
                   TDHF. A method for the amplitudes, not the states: the
                   states are those of HF (see calculate_formFactors_RPA).
*/
enum class AtomicMethod { HF, Zeff, ZeffAnalytic, RPA };

//------------------------------------------------------------------------------
/*!
  @brief Parses string (HF, Zeff, ZeffAnalytic, RPA) to AtomicMethod
  (case-insensitive).

  @param in_method  Method name; unknown input warns and defaults to HF.
  @return The corresponding AtomicMethod.
*/
AtomicMethod parseStatesMethod(const std::string &in_method);

//------------------------------------------------------------------------------
//! AtomicMethod to string (HF, Zeff, ZeffAnalytic, RPA)
std::string parseStatesMethod(const AtomicMethod &in_method);

//------------------------------------------------------------------------------
/*!
  @brief Effective charge of an orbital, from its binding energy.

  @details
  \f[ Z_{\rm eff} = n\sqrt{-2\en} \f]
  Same Zeff as used by DarkARC (see arXiv:1912.08204).

  @param en  Orbital energy (binding energy, negative), in au.
  @param n   Principal quantum number.
  @return Effective charge.
*/
inline double Zeff_nonrel(double en, int n) {
  return n * std::sqrt(std::abs(2.0 * en));
}

//------------------------------------------------------------------------------
/*!
  @brief Checks the radial grid is dense enough for the continuum states and
  momentum transfers required.

  @details
  Two independent checks, both of which print advice (larger num_points, or a
  different loglinear b) when they fail:
  - q: the grid spacing near \f$ r\sim a_0 \f$ must resolve the oscillations
    of \f$ e^{i\vb{q}\cdot\vb{r}} \f$ at @p qmax. Only a rough guide: high q
    may contribute negligibly, in which case error there does not matter.
  - E: the grid must resolve the continuum oscillations out to rmax at
    @p Emax (see DiracODE::RequiredContinuumGrid).

  @param Emax   Maximum continuum state energy, in au.
  @param qmax   Maximum momentum transfer, in au.
  @param rgrid  Radial grid to be checked (loglinear expected).
  @param alpha  Fine-structure constant (as used by the Hartree-Fock).
  @return False if either check fails; calculations may then be inaccurate.
*/
bool check_radial_grid(double Emax, double qmax, const Grid &rgrid,
                       double alpha = PhysConst::alpha);

//! The q part of check_radial_grid only: grid spacing near r ~ a0 must
//! resolve exp(iq.r) at qmax (au). Returns false (with advice) if not.
bool check_radial_grid_q(double qmax, const Grid &rgrid);

//! The E part of check_radial_grid only: grid must resolve the continuum
//! oscillations out to rmax at Emax (au). Returns false (with advice) if not.
bool check_radial_grid_E(double Emax, const Grid &rgrid,
                         double alpha = PhysConst::alpha);

//------------------------------------------------------------------------------
//! The 13 form factors of one bound orbital (or their sum), in the fixed
//! order: {V_T, V_E, V_M, V_L, X, A_T, A_E, A_M, A_L, Y, Z, S, P}. A factor
//! that is not calculated is left empty (0x0).
using FormFactorSet = std::array<LinAlg::Matrix<double>, 13>;

//------------------------------------------------------------------------------
//! Reduced matrix elements <e||h||a> of one (bound, continuum) channel for
//! each operator of the multipole set, in the order of multipole_operators():
//! {t, E, M, L, t5, E5, M5, L5, S, S5}. Real (imaginary part zero) without
//! RPA; complex (outgoing-wave amplitudes) with RPA. Zero if not calculated.
using ChannelAmplitudes = std::array<std::complex<double>, 10>;

//------------------------------------------------------------------------------
//! Options for the RPA (core polarisation) form factors
struct RPAOptions {
  //! Maximum RPA iterations per solve; 1 gives the first-order correction
  int max_its{60};
  //! RPA convergence target
  double eps{1.0e-10};
  //! An RPA solve whose final eps is above this (or nan) is discarded: the
  //! bare (no-RPA) amplitude is used for that (E, q, K, operator)
  double eps_fail{1.0e-3};
  //! The RPA is solved only for energy transfers E up to this (au); above
  //! it the factors are the bare ones (see calculate_formFactors_RPA)
  double E_max{1.0e99};
  //! The RPA is solved only for momentum transfers q up to this (au); in the
  //! diagonal (q = E/c) case this is a limit on E
  double q_max{1.0e99};
  //! The RPA is solved only for multipoles of rank K up to this; the higher
  //! multipoles contribute their bare factors
  int K_max{999};
};

//------------------------------------------------------------------------------
//! Bare and RPA form factors of every core orbital, from
//! calculate_formFactors_RPA()
struct FormFactorsRPA {
  //! Bare (Hartree-Fock) factors of each core orbital, indexed as the core
  std::vector<FormFactorSet> bare{};
  //! Factors including RPA (the total, not the correction), indexed as core
  std::vector<FormFactorSet> rpa{};
  //! Worst RPA eps over K and operators at each (E, q); zero where no
  //! orbital is ionised, or where the RPA was not solved (outside the
  //! RPAOptions limits)
  LinAlg::Matrix<double> eps{};
  //! Number of RPA solves made (one per (E, K, operator, q) with a non-zero
  //! amplitude, within the limits)
  std::size_t n_solves{0};
  //! How many of those did not converge (eps above RPAOptions::eps_fail, or
  //! nan); zero for a first-order solve, whose eps is not tested
  std::size_t n_failed{0};
};

//------------------------------------------------------------------------------
/*!
  @brief Small helper to allocate/size the requested factors.

  @details
  Each requested factor is allocated (E_steps x q_steps) and zeroed; the
  others are left empty (0x0), and are skipped by every function that takes
  a FormFactorSet. The interference terms X and Y follow their vector and
  axial parts (spatial only); Z requires vector, axial, and spatial.

  @param E_steps        Number of energy-transfer grid points (rows).
  @param q_steps        Number of momentum-transfer grid points (columns).
  @param vectorQ        Include the vector factors.
  @param axialQ         Include the axial-vector factors.
  @param scalarQ        Include the scalar factor.
  @param pseudoscalarQ  Include the pseudoscalar factor.
  @param spatialQ       Include the spatial (E, M, L) components; if false,
                        only the temporal components are allocated.
  @return The allocated set, in the fixed FormFactorSet order.
*/
FormFactorSet allocate_formFactors(std::size_t E_steps, std::size_t q_steps,
                                   bool vectorQ, bool axialQ, bool scalarQ,
                                   bool pseudoscalarQ, bool spatialQ);

//------------------------------------------------------------------------------
/*!
  @brief Constructs the required multipole operator set for the form factors (at w=q=0)

  @details
  Returned in the fixed order {t, E, M, L, t5, E5, M5, L5, S, S5}, with
  nullptr for those not requested: vector temporal t, electric E, magnetic M,
  longitudinal L; axial, their \f$ \gamma^5 \f$ partners; scalar S;
  pseudoscalar S5.

  Each is constructed at rank 0 and zero frequency: updateRank() then
  updateFrequency() must be called before use (see
  DiracOperator::MultipoleOperator for the units of the frequency: qc, or E
  in the diagonal case).

  @param grid           Radial grid on which the operators act.
  @param low_q          Use the low-momentum (long-wavelength) form.
  @param jK_tab         Precomputed spherical Bessel table (may be nullptr).
  @param vectorQ        Build the vector operators.
  @param axialQ         Build the axial-vector operators.
  @param scalarQ        Build the scalar operator.
  @param pseudoscalarQ  Build the pseudoscalar operator.
  @param spatialQ       Build the spatial (E, M, L) components.
  @return The operator set; entries not requested are nullptr.
*/
std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10>
multipole_operators(const Grid &grid, bool low_q,
                    const SphericalBessel::JL_table *jK_tab, bool vectorQ,
                    bool axialQ, bool scalarQ, bool pseudoscalarQ,
                    bool spatialQ);

//------------------------------------------------------------------------------
/*!
  @brief Adds the contribution of one (bound, continuum) channel to the form
  factors at (iE, iq).

  @details
  With weight \f$ w = (2K+1)x_{\rm occ} \f$ (@p tkp1_x), the squared terms
  get \f$ w|A|^2 \f$ and the interference terms
  \f$ w\,{\rm Re}(A\,A'^*) \f$: X from (t, L), Y from (t5, L5), and Z from
  (E5, M) - (E, M5). For real (bare) amplitudes this is the plain product.

  @param K_factors  Factors to add to; empty (not allocated) ones are skipped.
  @param iE         Energy-transfer grid index.
  @param iq         Momentum-transfer grid index.
  @param tkp1_x     Weight (2K+1) times the occupation fraction.
  @param A          Channel amplitudes, in multipole_operators() order.
*/
void accumulate_formFactors(FormFactorSet *K_factors, std::size_t iE,
                            std::size_t iq, double tkp1_x,
                            const ChannelAmplitudes &A);

//------------------------------------------------------------------------------
/*!
  @brief The bound state of each core orbital as used in the matrix elements:
  the real (HF) orbital, or its H-like (Zeff) version.

  @details
  For AtomicMethod::Zeff the H-like state is solved numerically (DiracODE) in
  the pointlike -Zeff/r potential, for AtomicMethod::ZeffAnalytic it is the
  exact Dirac-Coulomb function; otherwise the HF orbital itself. Zeff is
  that of each orbital from its binding energy (Zeff_nonrel(), as DarkARC). Every state carries the
  occupation of the real orbital. The Zeff states have the H-like energy,
  not the HF one: continuum energies must be taken from the real orbital.

  @param vHF            Hartree-Fock; its core defines the orbitals.
  @param method  States method (see AtomicMethod).
  @return One bound state per core orbital, indexed as the core.
*/
std::vector<DiracSpinor> model_bound_states(const HF::HartreeFock &vHF,
                                            AtomicMethod method);

//------------------------------------------------------------------------------
/*!
  @brief Adds the multipole sum of one bound orbital at one energy transfer
  to the form factors: every rank K, momentum transfer, and continuum state.

  @details
  For each K in [@p Kmin, @p Kmax] and each (column, qc) of @p q_columns, sets
  the rank and frequency of every (non-null) operator, forms the channel
  amplitudes $ 
edmatel{e}{h}{a} $ with each continuum state e, and
  accumulates them with weight $ (2K+1)x_{
m occ} $ (see
  accumulate_formFactors()) into row @p iE of @p K_factors. The operators
  take qc (or E itself in the diagonal, massless-absorption, case) as their
  frequency; see DiracOperator::MultipoleOperator.

  @param K_factors  Factors to add to; empty (not allocated) ones are skipped.
  @param iE         Row (energy-transfer index) of @p K_factors to add to.
  @param Fa         Bound orbital (as used in the matrix elements).
  @param occ_frac   Its occupation fraction.
  @param continuum  Continuum states of the ejected electron (all l, kappa).
  @param multipoles Operator set (multipole_operators()); modified in place
                    (rank and frequency), so each thread needs its own.
  @param Kmin       Minimum multipolarity K.
  @param Kmax       Maximum multipolarity K.
  @param q_columns  The (column index, frequency qc) pairs to evaluate.
*/
void accumulate_multipole_sum(
  FormFactorSet *K_factors, std::size_t iE, const DiracSpinor &Fa,
  double occ_frac, const std::vector<DiracSpinor> &continuum,
  std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10> &multipoles,
  int Kmin, int Kmax,
  const std::vector<std::pair<std::size_t, double>> &q_columns);

//------------------------------------------------------------------------------
/*!
  @brief Calculates all 13 form factors (V, A, S, P) for every core orbital.

  @details
  The continuum states of each orbital are solved at every E for which the
  ejected electron energy \f$ \en_c = E + \en_a \f$ lies in
  (@p ec_min, @p ec_max]. The continuum l are those reached from the orbital
  by multipoles of rank up to @p Kmax (both parities),
  \f$ j_e = j_a \pm K_{\rm max} \f$, \f$ l_e = j_e \pm 1/2 \f$, clipped to
  @p lc_minmax if given. Parallel over the (orbital, E) pairs.

  The bound state used in the matrix elements of each orbital is
  @p bound_states[ia]: the Hartree-Fock orbital itself, or (for the H-like
  methods) its Zeff version, from model_bound_states(). With @p method
  Zeff or ZeffAnalytic the continuum is also H-like, in the pointlike
  potential of the same charge as that model state (Zeff_nonrel() of the HF
  binding energy), solved numerically (DiracODE) or with the exact
  analytic Dirac-Coulomb functions. The
  occupation x_i is that of the bound state used; continuum energies always
  use the real (HF) orbital.

  @param vHF            Hartree-Fock potential; its core defines the orbitals.
  @param bound_states   Bound state of each core orbital as used in the
                        matrix elements, indexed as the core.
  @param lc_minmax      Optional limits on the continuum orbital l.
  @param ec_min         Minimum ejected electron energy, in au.
  @param ec_max         Maximum ejected electron energy, in au.
  @param force_rescale  Rescale V(r) at large r for the continuum states.
  @param hole_particle  Solve the continuum in the V^(N-1) potential of the
                        hole (include the hole-particle interaction).
  @param force_orthog   Orthogonalise the continuum states to the core.
  @param Egrid          Energy transfer grid, in au.
  @param qgrid          Momentum transfer grid, in au (size 1 if diagonal).
  @param diagonal_Eq    Momentum transfer set equal to the energy transfer
                        (absorption of a massless particle).
  @param low_q          Use the low-q form of the operators.
  @param jK_tab         Precomputed spherical Bessel table.
  @param Kmin           Minimum multipolarity K.
  @param Kmax           Maximum multipolarity K.
  @param vectorQ        Calculate the vector factors.
  @param axialQ         Calculate the axial-vector factors.
  @param scalarQ        Calculate the scalar factor.
  @param pseudoscalarQ  Calculate the pseudoscalar factor.
  @param spatialQ       Calculate the spatial (E, M, L) components.
  @param method         Continuum states:
                        - AtomicMethod::HF           : Hartree-Fock
                        - AtomicMethod::Zeff         : H-like, DiracODE
                        - AtomicMethod::ZeffAnalytic : H-like, analytic
                        (AtomicMethod::RPA is not a states method, and is
                        treated as HF here; see calculate_formFactors_RPA.)
  @return One FormFactorSet per core orbital, indexed as the core (zero for
          an orbital that is not ionised anywhere on the E grid). See
          allocate_formFactors() for which factors are calculated.

  @note @p force_rescale and @p hole_particle have no effect for Zeff states.
*/
std::vector<FormFactorSet> calculate_formFactors(
  const HF::HartreeFock *vHF, const std::vector<DiracSpinor> &bound_states,
  const std::optional<std::array<int, 2>> &lc_minmax, double ec_min,
  double ec_max, bool force_rescale, bool hole_particle, bool force_orthog,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, AtomicMethod method = AtomicMethod::HF);

//------------------------------------------------------------------------------
//! Adds the factors of each orbital of dK to those of K_nk (same orbitals,
//! same allocation pattern; factors not allocated in either are skipped)
void add_formFactors(std::vector<FormFactorSet> *K_nk,
                     const std::vector<FormFactorSet> &dK);

//------------------------------------------------------------------------------
//! One ionisation channel: a hole in a core orbital, and the ejected
//! (continuum) state
struct IonisationChannel {
  //! Index of the hole orbital in the core
  std::size_t hole_index{};
  //! The hole (core) orbital
  const DiracSpinor *hole{nullptr};
  //! The ejected (continuum) state
  const DiracSpinor *ejected{nullptr};
};

//------------------------------------------------------------------------------
/*!
  @brief Range of continuum l reached from a bound orbital by multipoles of
  rank up to Kmax.

  @details
  Either parity: \f$ j_e = j_a \pm K_{\rm max} \f$, \f$ l_e = j_e \pm 1/2 \f$,
  clipped to @p lc_minmax if given.

  @param Fa         Bound (core) orbital.
  @param Kmax       Maximum multipolarity.
  @param lc_minmax  Optional limits {lc_min, lc_max} on the continuum l.
  @return {lc_min, lc_max}.
*/
std::pair<int, int>
continuum_l_range(const DiracSpinor &Fa, int Kmax,
                  const std::optional<std::array<int, 2>> &lc_minmax);

//------------------------------------------------------------------------------
//! A core orbital ionised at one energy transfer, with the continuum states
//! of its ejected electron
struct IonisedOrbital {
  //! Index of the orbital in the core
  std::size_t core_index{};
  //! Continuum states of the ejected electron
  ContinuumOrbitals ejected;
};

//------------------------------------------------------------------------------
/*!
  @brief The core orbitals ionised by an energy transfer omega, each with the
  Hartree-Fock continuum states of its ejected electron.

  @details
  An orbital is ionised if its ejected electron energy
  \f$ \en_c = \omega + \en_a \f$ lies in (@p ec_min, @p ec_max]. Its continuum
  states are solved at that energy, with l from continuum_l_range(). In core
  order; parallel over the orbitals.

  @param vHF            Hartree-Fock potential; its core defines the orbitals.
  @param omega          Energy transfer, in au.
  @param ec_min         Minimum ejected electron energy, in au.
  @param ec_max         Maximum ejected electron energy, in au.
  @param Kmax           Maximum multipolarity (sets the continuum l range).
  @param lc_minmax      Optional limits on the continuum orbital l.
  @param force_rescale  Rescale V(r) at large r for the continuum states.
  @param hole_particle  Solve the continuum in the V^(N-1) potential of the
                        hole.
  @param force_orthog   Orthogonalise the continuum states to the core.
  @return The ionised orbitals and their continuum states; empty if none is
          ionised.
*/
std::vector<IonisedOrbital> solve_ionised_orbitals_at_omega(
  const HF::HartreeFock *vHF, double omega, double ec_min, double ec_max,
  int Kmax, const std::optional<std::array<int, 2>> &lc_minmax,
  bool force_rescale, bool hole_particle, bool force_orthog);

//------------------------------------------------------------------------------
/*!
  @brief Every (hole orbital, ejected state) channel at one energy transfer.

  @details
  Ordered by orbital (as @p ionised), then by continuum state. The channels
  point into @p core and @p ionised, which must outlive them.

  @param core     Core orbitals.
  @param ionised  The ionised orbitals and their continuum states
                  (solve_ionised_orbitals_at_omega()).
  @return The channel list.
*/
std::vector<IonisationChannel>
construct_channels(const std::vector<DiracSpinor> &core,
                   const std::vector<IonisedOrbital> &ionised);

//------------------------------------------------------------------------------
//! Indices of the requested (non-null) operators of a multipole_operators()
//! set, in its order
std::vector<std::size_t> active_operators(
  const std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10>
    &operators);

//------------------------------------------------------------------------------
/*!
  @brief One RPA solve: the bare and RPA amplitudes of every channel, for one
  operator at one (E, K, q).

  @details
  @p h is the operator @p rpa was constructed with, with its rank and
  frequency already set. Solves the TDHF at @p omega (warm starting from the
  solver's current state), then fills amplitude @p i_op of each channel:
  bare \f$ \redmatel{e}{h}{a} \f$ in @p A_bare, and RPA
  \f$ \redmatel{e}{h + \delta V}{a} \f$ (ExternalField::TDHFcntm::dV_complex)
  in @p A_rpa. Channels for which the operator is zero by selection rules are
  left as they are.

  If the solve did not converge (eps above RPAOptions::eps_fail, or nan), the
  RPA amplitude is the bare one, and the solver is cleared so that the next
  solve does not warm start from the failed state. Convergence is not tested
  for a first-order solve (RPAOptions::max_its of 1): its eps is the size of
  the correction, not a convergence measure. An operator with no non-zero
  bare amplitude in any channel (a transverse multipole at K = 0) is not
  solved: the RPA amplitudes are the bare ones, and nothing is returned.

  @param h            Multipole operator (rank and frequency set).
  @param i_op         Its index in the multipole_operators() set.
  @param rpa          TDHF solver for @p h.
  @param omega        Energy transfer, in au.
  @param rpa_options  Iterations and convergence limits (see RPAOptions).
  @param channels     The (hole, ejected) channels (construct_channels()).
  @param A_bare       Bare amplitudes, one per channel; entry @p i_op set.
  @param A_rpa        RPA amplitudes, one per channel; entry @p i_op set.
  @param print        Print the RPA iterations of the solve.
  @return The eps of the solve (see ExternalField::TDHF::last_eps); nullopt
          if the RPA was not solved.
*/
std::optional<double> solve_channel_amplitudes(
  const DiracOperator::TensorOperator &h, std::size_t i_op,
  ExternalField::TDHFcntm *rpa, double omega, const RPAOptions &rpa_options,
  const std::vector<IonisationChannel> &channels,
  std::vector<ChannelAmplitudes> *A_bare, std::vector<ChannelAmplitudes> *A_rpa,
  bool print = false);

//------------------------------------------------------------------------------
/*!
  @brief Bare and RPA (core polarisation) form factors of every core orbital,
  with the RPA solved at every point of the grids, from the outgoing-wave
  TDHF.

  @details
  The bare factors are exactly those of calculate_formFactors() (Hartree-Fock
  states only). For the RPA, each channel amplitude
  \f$ \redmatel{e}{h}{a} \f$ is replaced by the complex outgoing-wave
  amplitude \f$ \redmatel{e}{h + \delta V}{a} \f$ (see
  ExternalField::TDHFcntm::dV_complex), and the factors are accumulated as in
  accumulate_formFactors().

  One RPA solve is required per (E, K, operator, q); it serves every core
  orbital at once, which is why the orbital loop is inside. The RPA includes
  every open channel; the ec limits apply to the output only.

  Parallel over the coarsest axis that keeps every thread busy. Over the
  energies at which some orbital is ionised, when there are at least as many
  as threads: each thread solves its energies in full, with its own
  continuum states and solver, and the q points in one warm-start chain.
  Otherwise, within each (E, K, operator) block, over q when there are at
  least as many as threads (each thread owns a solver, and its solves warm
  start from the neighbouring q); else q runs serially with the solver's own
  parallelism. Prints one summary line at the start, then a progress bar.

  The RPA is solved at every (E, q) given: the E and q limits of
  @p rpa_options are not applied here (see calculate_formFactors_RPA(), which
  restricts the grids to the region within them; what the formFactors module
  uses). The K limit is: multipoles above RPAOptions::K_max contribute their
  bare amplitudes, with no solve.

  @param vHF            Hartree-Fock potential; its core defines the orbitals.
  @param lc_minmax      Optional limits on the continuum orbital l.
  @param ec_min         Minimum ejected electron energy, in au.
  @param ec_max         Maximum ejected electron energy, in au.
  @param force_rescale  Rescale V(r) at large r for the continuum states.
  @param hole_particle  Solve the continuum in the V^(N-1) potential of the
                        hole (required here; see warning).
  @param force_orthog   Orthogonalise the continuum states to the core.
  @param Egrid          Energy transfer grid, in au.
  @param qgrid          Momentum transfer grid, in au (size 1 if diagonal).
  @param diagonal_Eq    Momentum transfer set equal to the energy transfer.
  @param low_q          Use the low-q form of the operators.
  @param jK_tab         Precomputed spherical Bessel table.
  @param Kmin           Minimum multipolarity K.
  @param Kmax           Maximum multipolarity K.
  @param vectorQ        Calculate the vector factors.
  @param axialQ         Calculate the axial-vector factors.
  @param scalarQ        Calculate the scalar factor.
  @param pseudoscalarQ  Calculate the pseudoscalar factor.
  @param spatialQ       Calculate the spatial (E, M, L) components.
  @param rpa_options    RPA iterations, convergence target, the eps above
                        which a solve is discarded, and the K limit (see
                        RPAOptions); the E and q limits are not applied.
  @return The bare and RPA factors of each core orbital, and the worst RPA
          eps at each (E, q).

  @note An unconverged solve (eps above RPAOptions::eps_fail, or nan) is not
        used: the bare amplitude is taken for that (E, K, operator, q), and
        the solver is cleared so the next q does not warm start from it (see
        solve_channel_amplitudes()). FormFactorsRPA::eps records the worst
        eps at each (E, q), for diagnostics.

  @warning The continuum states must be those of the residual ion
           (@p hole_particle = true) for the RPA amplitude to be consistent;
           see ExternalField::TDHFcntm::dV_complex.
*/
FormFactorsRPA solve_formFactors_RPA(
  const HF::HartreeFock *vHF,
  const std::optional<std::array<int, 2>> &lc_minmax, double ec_min,
  double ec_max, bool force_rescale, bool hole_particle, bool force_orthog,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, const RPAOptions &rpa_options = {});

//------------------------------------------------------------------------------
/*!
  @brief Calculates the form factors for every core orbital, without and with
  RPA (core polarisation), the RPA solved only within the limits of
  RPAOptions.

  @details
  The bare factors are those of calculate_formFactors() (Hartree-Fock states)
  over the full grids. The RPA is solved (solve_formFactors_RPA()) only on
  the region of the grids with \f$ E \le E_{\rm max} \f$ and
  \f$ q \le q_{\rm max} \f$ (RPAOptions::E_max, q_max), and its factors
  replace the bare ones there; elsewhere the RPA factors are the bare ones.
  Multipoles above RPAOptions::K_max are bare everywhere (applied inside
  solve_formFactors_RPA()). In the diagonal case (q = E/c) the q limit is a
  limit on E. With no E or q limit in effect the RPA is solved over the full
  grids directly, and the bare factors are calculated once only.

  The limits exist because the RPA is expensive, and hard to converge, at
  high E, q, and K (many open channels; rapidly oscillating operators) where
  its effect on the factors is small: they confine the work to the low-E,
  low-q, low-K region where the core polarisation matters.

  @param vHF            Hartree-Fock potential; its core defines the orbitals.
  @param lc_minmax      Optional limits on the continuum orbital l.
  @param ec_min         Minimum ejected electron energy, in au.
  @param ec_max         Maximum ejected electron energy, in au.
  @param force_rescale  Rescale V(r) at large r for the continuum states.
  @param hole_particle  Solve the continuum in the V^(N-1) potential of the
                        hole (required here; see warning).
  @param force_orthog   Orthogonalise the continuum states to the core.
  @param Egrid          Energy transfer grid, in au.
  @param qgrid          Momentum transfer grid, in au (size 1 if diagonal).
  @param diagonal_Eq    Momentum transfer set equal to the energy transfer.
  @param low_q          Use the low-q form of the operators.
  @param jK_tab         Precomputed spherical Bessel table.
  @param Kmin           Minimum multipolarity K.
  @param Kmax           Maximum multipolarity K.
  @param vectorQ        Calculate the vector factors.
  @param axialQ         Calculate the axial-vector factors.
  @param scalarQ        Calculate the scalar factor.
  @param pseudoscalarQ  Calculate the pseudoscalar factor.
  @param spatialQ       Calculate the spatial (E, M, L) components.
  @param rpa_options    RPA iterations, convergence target, the eps above
                        which a solve is discarded, and the E, q, K limits
                        (see RPAOptions).
  @return The bare and RPA factors of each core orbital, and the worst RPA
          eps at each (E, q): zero outside the RPA region.

  @note An unconverged solve within the region is treated as in
        solve_formFactors_RPA(): the bare amplitude is used there, and
        FormFactorsRPA::eps records it.

  @warning The continuum states must be those of the residual ion
           (@p hole_particle = true) for the RPA amplitude to be consistent;
           see ExternalField::TDHFcntm::dV_complex.
*/
FormFactorsRPA calculate_formFactors_RPA(
  const HF::HartreeFock *vHF,
  const std::optional<std::array<int, 2>> &lc_minmax, double ec_min,
  double ec_max, bool force_rescale, bool hole_particle, bool force_orthog,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, const RPAOptions &rpa_options = {});

//------------------------------------------------------------------------------
//! A failed RPA solve: final eps above eps_fail, or nan
inline bool rpa_failed(double eps, double eps_fail) {
  return std::isnan(eps) || eps > eps_fail;
}

//------------------------------------------------------------------------------
/*!
  @brief Counts the points of an RPA eps table at which the solve failed, and
  those of them with a converged neighbour along the columns.

  @details
  A point fails if rpa_failed(); its neighbours are the adjacent columns of
  the same row (the q grid for the form factors, the energy grid for
  photoionisation). Those with at least one converged neighbour are the
  points interpolate_failed_rpa() corrects.

  @param eps       Final RPA eps of each point (rows x columns).
  @param eps_fail  Failure threshold; as RPAOptions::eps_fail.
  @return {number of failed points, number of those with a converged
          neighbour}.
*/
std::pair<std::size_t, std::size_t>
count_failed_rpa(const LinAlg::Matrix<double> &eps, double eps_fail);

//------------------------------------------------------------------------------
/*!
  @brief Corrects the points at which the RPA solve failed, by interpolating
  the relative RPA shift from the neighbouring columns.

  @details
  Where a solve fails the bare value is used, which leaves a step in an
  otherwise smooth function. The relative shift \f$ R = K_{\rm RPA}/K_{\rm
  bare} \f$ varies slowly along the columns (q for the form factors, the
  photon energy for photoionisation), so it is taken from the converged
  neighbours and applied to the bare value,
  \f$ K_{\rm RPA}(i,j) = R\,K_{\rm bare}(i,j) \f$:
  \f[
  \begin{align}
    R &= \frac{1}{2}(R_{j-1} + R_{j+1}), \\
    R &= 1 + \frac{1}{2}(R_{j\mp1} - 1),
  \end{align}
  \f]
  the mean of the two when both have converged; half the relative correction
  when only one has, since the shift is then unconstrained on the other side.
  A point with no converged neighbour (in a run of adjacent failures, a
  resonance say) is left as it is, as is one whose converged neighbours have
  a zero bare value. Points that fail are never used as neighbours, so the
  result does not depend on the order of the corrections.

  @param eps       Final RPA eps of each point (rows x columns).
  @param eps_fail  Failure threshold (see rpa_failed()).
  @param K_bare    Bare (no-RPA) values, same shape as @p eps.
  @param K_rpa     RPA values, same shape; corrected in place.

  @note Corrects nothing with a single column, since there is no neighbour
        to interpolate from. Only meaningful for an iterated solve: after a
        first-order solve (RPAOptions::max_its of 1) eps is the size of the
        correction, not a convergence measure.
*/
void interpolate_failed_rpa(const LinAlg::Matrix<double> &eps, double eps_fail,
                            const LinAlg::Matrix<double> &K_bare,
                            LinAlg::Matrix<double> *K_rpa);

//------------------------------------------------------------------------------
/*!
  @brief Corrects the (E,q) points at which the RPA solve failed, for every
  factor of every orbital, by interpolating in q.

  @details
  interpolate_failed_rpa() (the matrix form) applied to each allocated
  factor, with the neighbours along q; see there for the formula.

  @param K_rpa    Bare and RPA factors, and the eps of each (E,q); the RPA
                  factors are updated in place.
  @param eps_fail Failure threshold; as RPAOptions::eps_fail.
  @return {number of (E,q) points that failed, number of those corrected};
          the rest keep the no-RPA value (see count_failed_rpa()).
*/
std::pair<std::size_t, std::size_t>
interpolate_failed_rpa(FormFactorsRPA *K_rpa, double eps_fail);

//------------------------------------------------------------------------------
/*!
  @brief Calculates the ionisation factor K(E,q) for one core state, using the
  standard (single multipole operator) method. New code should use calculate_formFactors

  @details
  \f[ 
    K(E,q) = \sum_{L,e} (2L+1)\,x_{\rm occ}\,|\redmatel{e}{j_L}{a}|^2 
  \f]
  summed over the multipoles L up to @p max_L and the continuum states e with
  \f$ l_e \f$ within @p max_L of \f$ l_a \f$. The continuum energy is
  \f$ \en_c = E + \en_a \f$; grid points at which this is not positive (or
  exceeds @p ec_cut) are left zero. Parallelised over E or q, whichever is
  larger.

  @note Should be equivilant to temporal component of calculate_formFactors.
  Prefer calculate_formFactors() or calculate_formFactors_RPA() for new code.

  @param vHF            Hartree-Fock potential.
  @param Fnk            Bound (core) orbital being ionised.
  @param max_L          Maximum multipolarity L.
  @param Egrid          Energy transfer grid, in au.
  @param jl             Operator providing the reduced matrix elements; its
                        q grid sets the columns.
  @param force_rescale  Rescale V(r) at large r for the continuum states.
  @param hole_particle  Solve the continuum in the V^(N-1) potential.
  @param force_orthog   Orthogonalise the continuum states to the core.
  @param zeff_cont      Use H-like (Zeff) continuum states.
  @param zeff_bound     Use an H-like (Zeff) bound state in the matrix
                        elements (the HF energy is still used).
  @param ec_cut         Maximum continuum energy, in au.
  @return K(E,q) as a matrix: each row a new E, each column a new q.
*/
LinAlg::Matrix<double>
calculateK_nk(const HF::HartreeFock *vHF, const DiracSpinor &Fnk, int max_L,
              const Grid &Egrid, const DiracOperator::jL *jl,
              bool force_rescale, bool hole_particle, bool force_orthog,
              bool zeff_cont, bool zeff_bound, double ec_cut = 1.0e99);

//------------------------------------------------------------------------------
/*!
  @brief Writes output file in matrix form.

  @details
  The E and q grids are printed first, then the entire matrix in table form:
  each row is a new E, each column a new q.

  @note Kept for legacy: new code should use XYZ format

  @param K           Factor to write, K(E,q).
  @param E_grid      Energy transfer grid, in au.
  @param q_grid      Momentum transfer grid, in au.
  @param filename    Output file name.
  @param num_digits  Digits printed for each value.
  @param units       Units for E and q in the output (K is dimensionless).
*/
void write_to_file_matrix(const LinAlg::Matrix<double> &K,
                          const std::vector<double> &E_grid,
                          const std::vector<double> &q_grid,
                          const std::string &filename, int num_digits = 5,
                          Units units = Units::Particle);

//------------------------------------------------------------------------------
/*!
  @brief Writes output file in 'xyz' form: for easy 2D interpolation.

  @details
  A header of column descriptions, then one row per (E, q) point:
  'E q K_1(E,q) ... K_n(E,q)'. A blank line separates each energy (which
  gnuplot requires, and pyplot ignores). Factors that are empty are skipped,
  along with their column.

  @param filename      Output file name.
  @param E_grid        Energy transfer grid, in au.
  @param q_grid        Momentum transfer grid, in au.
  @param titles        Short column header of each factor.
  @param descriptions  Longer description of each factor, for the header.
  @param factors       Factors to write; must match @p titles in size.
  @param units         Units for E and q (K is dimensionless).
  @param num_digits    Digits printed for each value (clamped to 3 - 16).
  @param diagonal      Momentum transfer equal to the energy transfer: q is
                       written as alpha*E, and @p q_grid is not used.
*/
void write_to_file_xyz(const std::string &filename,
                       const std::vector<double> &E_grid,
                       const std::vector<double> &q_grid,
                       const std::vector<std::string> &titles,
                       const std::vector<std::string> &descriptions,
                       std::vector<LinAlg::Matrix_view<const double>> factors,
                       Units units = Units::Particle, int num_digits = 6,
                       bool diagonal = false);

//------------------------------------------------------------------------------
//! As write_to_file_xyz, for a FormFactorSet (empty factors are skipped)
void write_to_file_xyz_13(
  const std::string &filename, const std::vector<double> &E_grid,
  const std::vector<double> &q_grid, const std::vector<std::string> &titles,
  const std::vector<std::string> &descriptions, const FormFactorSet &K_factors,
  Units units = Units::Particle, int num_digits = 6, bool diagonal = false);

} // namespace Kion
