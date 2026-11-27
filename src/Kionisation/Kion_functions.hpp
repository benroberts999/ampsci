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
  different loglinear b) when they fail (nothing printed if @p print is
  false):
  - q: the grid spacing near \f$ r\sim a_0 \f$ must resolve the oscillations
    of \f$ e^{i\vb{q}\cdot\vb{r}} \f$ at @p qmax. Only a rough guide: high q
    may contribute negligibly, in which case error there does not matter.
  - E: the grid must resolve the continuum oscillations out to rmax at
    @p Emax (see DiracODE::RequiredContinuumGrid).

  @param Emax   Maximum continuum state energy, in au.
  @param qmax   Maximum momentum transfer, in au.
  @param rgrid  Radial grid to be checked (loglinear expected).
  @param alpha  Fine-structure constant (as used by the Hartree-Fock).
  @param print  Print the result of each check, and the advice.
  @return The num_points (same r0, rmax, b) with which both checks pass: that
          of @p rgrid if they pass already; more if either fails, in which
          case calculations on @p rgrid may be inaccurate.
*/
std::size_t check_radial_grid(double Emax, double qmax, const Grid &rgrid,
                              double alpha = PhysConst::alpha,
                              bool print = true);

//! The q part of check_radial_grid only: grid spacing near r ~ a0 must
//! resolve exp(iq.r) at qmax (au). Returns the num_points required (same r0,
//! rmax, b): that of the grid if it passes; more (with advice) if not.
std::size_t check_radial_grid_q(double qmax, const Grid &rgrid,
                                bool print = true);

//! The E part of check_radial_grid only: grid must resolve the continuum
//! oscillations out to rmax at Emax (au). Returns the num_points required
//! (same r0, rmax, b): that of the grid if it passes; more (with advice) if
//! not.
std::size_t check_radial_grid_E(double Emax, const Grid &rgrid,
                                double alpha = PhysConst::alpha,
                                bool print = true);

//------------------------------------------------------------------------------
/*!
  @brief Checks the momentum-transfer grid is dense enough to resolve the
  Bethe ridge.

  @details
  At energy transfer E, the form factor of orbital a peaks at the momentum
  of the ejected electron (the Bethe ridge),
  \f[
    q \approx p_f, \qquad
    p_f^2 = \en_c\,(2 + \alpha^2\en_c), \qquad \en_c = E + \en_a,
  \f]
  with a width in q set by the momentum spread of the bound orbital, taken
  as
  \f[
    w_a = \sqrt{2|\en_a|}
  \f]
  (the full width at half maximum for a hydrogen-like 1s state, to within a
  few percent). The check passes if the q step at the ridge is at most
  half of this, for every orbital and every energy whose ridge lies within
  the q grid. Only a rough guide. Prints the step required for the narrowest
  ridge (the most loosely bound orbital), then for each decade of energy
  (E <= 100 eV, <= 1000 eV, ...; labelled by its upper boundary, or by the
  largest energy for the last) the q step at the ridge, for the worst case in
  that decade, and the number of orbitals for which the check fails; if it
  fails, also a warning with the lowest energy at which it fails and the
  q_steps required.

  @param core   Core orbitals.
  @param Egrid  Energy transfer grid, in au.
  @param qgrid  Momentum transfer grid, in au (any order).
  @param alpha  Fine-structure constant (as used by the Hartree-Fock).
  @return The number of q points with which the check passes, for a
          logarithmic grid over the same range: the size of @p qgrid if it
          passes already.
*/
std::size_t check_ridge_q_grid(const std::vector<DiracSpinor> &core,
                               const std::vector<double> &Egrid,
                               const std::vector<double> &qgrid,
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
  double eps{1.0e-8};
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
  //! The RPA is solved only at every E_stride-th point of the energy grid
  //! (within the limits); it is interpolated to the points between (see
  //! calculate_formFactors_RPA). 1: every point
  std::size_t E_stride{1};
  //! As E_stride, for the momentum transfer grid
  std::size_t q_stride{1};
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
  //! RPAOptions limits, or a point skipped by the E or q stride)
  LinAlg::Matrix<double> eps{};
  //! Number of RPA solves made (one per (E, K, operator, q) with a non-zero
  //! amplitude, within the limits)
  std::size_t n_solves{0};
  //! How many of those did not converge (eps above RPAOptions::eps_fail, or
  //! nan) or were too large (check_rpa_spikes()); zero for a
  //! first-order solve, whose eps is not tested
  std::size_t n_failed{0};
  //! How many of the failed solves had a converged neighbour (same K and
  //! operator, adjacent E or q) to interpolate the RPA shift from
  std::size_t n_interpolated{0};
  //! Energy transfers at which the RPA was solved, in au (interpolated
  //! between, within the RPA limits); empty if it was not solved
  std::vector<double> Egrid_solved{};
  //! As Egrid_solved, for the momentum transfers
  std::vector<double> qgrid_solved{};
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
  amplitudes \f$ \redmatel{e}{h}{a} \f$ with each continuum state e, and
  accumulates them with weight \f$ (2K+1)x_{\rm occ} \f$ (see
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
  ejected electron energy \f$ \en_c = E + \en_a \f$ is above
  @p ec_min. The continuum l are those reached from the orbital
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
  @param print_time     Print the time taken on completion (false when the
                        caller times a larger calculation).
  @return One FormFactorSet per core orbital, indexed as the core (zero for
          an orbital that is not ionised anywhere on the E grid). See
          allocate_formFactors() for which factors are calculated.

  @note @p force_rescale and @p hole_particle have no effect for Zeff states.
*/
std::vector<FormFactorSet> calculate_formFactors(
  const HF::HartreeFock *vHF, const std::vector<DiracSpinor> &bound_states,
  const std::optional<std::array<int, 2>> &lc_minmax, double ec_min,
  bool force_rescale, bool hole_particle, bool force_orthog,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, AtomicMethod method = AtomicMethod::HF,
  bool print_time = true);

//------------------------------------------------------------------------------
//! Adds the factors of each orbital of dK to those of K_nk (same orbitals,
//! same allocation pattern; factors not allocated in either are skipped)
void add_formFactors(std::vector<FormFactorSet> *K_nk,
                     const std::vector<FormFactorSet> &dK);

//------------------------------------------------------------------------------
//! Copies the factors of each orbital of K_part, calculated at some of the
//! energies, into those of K_nk (same orbitals, same momentum transfers): row
//! j of K_part goes to row E_indices[j] of K_nk. Factors not allocated in
//! either are skipped
void insert_formFactors(std::vector<FormFactorSet> *K_nk,
                        const std::vector<FormFactorSet> &K_part,
                        const std::vector<std::size_t> &E_indices);

//------------------------------------------------------------------------------
/*!
  @brief Splits a list of energies into decades.

  @details
  The first decade is every energy up to @p first_decade_max; each one after
  spans a factor of 10 (up to 10, 100, ... times @p first_decade_max). Decades
  with no energies are left out.

  @param energies          The energies (any order).
  @param first_decade_max  Upper limit of the first decade (same units).
  @return For each decade, in ascending order of energy, the indices into
          @p energies of its energies (in the order of @p energies).
*/
std::vector<std::vector<std::size_t>>
energy_decades(const std::vector<double> &energies, double first_decade_max);

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
  @brief The core orbitals ionised by each of a set of energy transfers, each
  with the Hartree-Fock continuum states of its ejected electron.

  @details
  An orbital is ionised if its ejected electron energy
  \f$ \en_c = \omega + \en_a \f$ is above @p ec_min. Its continuum
  states are solved at that energy, with l from continuum_l_range(), the
  unresolved high-r tail locally averaged (DiracODE::averageTail). In core
  order at each energy; parallel over (energy, orbital, l).

  @param vHF            Hartree-Fock potential; its core defines the orbitals.
  @param omegas         Energy transfers, in au.
  @param ec_min         Minimum ejected electron energy, in au.
  @param Kmax           Maximum multipolarity (sets the continuum l range).
  @param lc_minmax      Optional limits on the continuum orbital l.
  @param force_rescale  Rescale V(r) at large r for the continuum states.
  @param hole_particle  Solve the continuum in the V^(N-1) potential of the
                        hole.
  @param force_orthog   Orthogonalise the continuum states to the core.
  @return At each energy (as @p omegas), the ionised orbitals and their
          continuum states; empty if none is ionised.
*/
std::vector<std::vector<IonisedOrbital>> solve_ionised_orbitals(
  const HF::HartreeFock *vHF, const std::vector<double> &omegas, double ec_min,
  int Kmax, const std::optional<std::array<int, 2>> &lc_minmax,
  bool force_rescale, bool hole_particle, bool force_orthog);

//------------------------------------------------------------------------------
/*!
  @brief Every (hole orbital, ejected state) channel at one energy transfer.

  @details
  Ordered by orbital (as @p ionised), then by continuum state. The channels
  point into @p core and @p ionised, which must outlive them.

  @param core     Core orbitals.
  @param ionised  The ionised orbitals and their continuum states at one
                  energy (solve_ionised_orbitals()).
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
  @brief True if an RPA solve changed the amplitudes of some orbital by more
  than a factor of 10: such a solve is spurious.

  @details
  The size of the amplitudes of a hole orbital is the root of the sum of
  \f$ |\redmatel{e}{h}{a}|^2 \f$ over its channels (what enters its form
  factor), so that a channel whose bare amplitude happens to vanish does not
  matter. Tested for each hole orbital, with RPA against without. A nan RPA
  amplitude also counts as too large.

  @param i_op      Operator index in the multipole_operators() set.
  @param channels  The (hole, ejected) channels (construct_channels()).
  @param A_bare    Bare amplitudes, one per channel.
  @param A_rpa     RPA amplitudes, one per channel.
  @return True if the RPA amplitudes of any hole orbital are too large.
*/
bool check_rpa_spikes(std::size_t i_op,
                      const std::vector<IonisationChannel> &channels,
                      const std::vector<ChannelAmplitudes> &A_bare,
                      const std::vector<ChannelAmplitudes> &A_rpa);

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
  the correction, not a convergence measure. A solve whose amplitudes are
  too large (check_rpa_spikes()) is discarded in the same way, and
  nan is returned for its eps. An operator with no non-zero
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
  @return The eps of the solve (see ExternalField::TDHF::last_eps): nan if
          the solve was discarded as too large; nullopt if the RPA was not
          solved.
*/
std::optional<double> solve_channel_amplitudes(
  const DiracOperator::TensorOperator &h, std::size_t i_op,
  ExternalField::TDHFcntm *rpa, double omega, const RPAOptions &rpa_options,
  const std::vector<IonisationChannel> &channels,
  std::vector<ChannelAmplitudes> *A_bare, std::vector<ChannelAmplitudes> *A_rpa,
  bool print = false);

//------------------------------------------------------------------------------
/*!
  @brief Bare and RPA (core polarisation) form factors of every core orbital
  from the multipoles Kmin to Kmax, with the RPA solved at every point of the
  grids, from the outgoing-wave TDHF.

  @details
  The bare factors are those of calculate_formFactors() (Hartree-Fock states
  only) for the same multipoles, so that rpa - bare is the RPA shift of those
  multipoles. For the RPA, each channel amplitude
  \f$ \redmatel{e}{h}{a} \f$ is replaced by the complex outgoing-wave
  amplitude \f$ \redmatel{e}{h + \delta V}{a} \f$ (see
  ExternalField::TDHFcntm::dV_complex), and the factors are accumulated as in
  accumulate_formFactors().

  One RPA solve is required per (E, K, operator, q); it serves every core
  orbital at once, which is why the orbital loop is inside. The RPA includes
  every open channel; the ec limits apply to the output only.

  Parallel over chains of solves (see the note below). Prints a progress bar
  over the solves, and nothing else.

  The RPA is solved at every (E, q) given and every K from Kmin to Kmax: the
  E, q and K limits of @p rpa_options are not applied here (see
  calculate_formFactors_RPA(), which restricts the grids and the K range to
  within them; what the formFactors module uses).

  @param vHF            Hartree-Fock potential; its core defines the orbitals.
  @param lc_minmax      Optional limits on the continuum orbital l.
  @param ec_min         Minimum ejected electron energy, in au.
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
  @param rpa_options    RPA iterations, convergence target, and the eps above
                        which a solve is discarded (see RPAOptions); the E, q
                        and K limits are not applied.
  @return The bare and RPA factors of each core orbital from the multipoles
          Kmin to Kmax, and the worst RPA eps at each (E, q).

  @note An unconverged solve (eps above RPAOptions::eps_fail, or nan), or one
        whose amplitudes are too large (check_rpa_spikes()), is not
        used: the bare amplitude is taken for that (E, K, operator, q), and
        the solver is cleared so the next q does not warm start from it (see
        solve_channel_amplitudes()). The factors are accumulated per K for
        the RPA-solved K, and each (K, factor) block is corrected at its
        failed (E, q) points from the converged neighbours in E and q
        (interpolate_failed_rpa_2d()) before the sum over K; not after a
        first-order solve. FormFactorsRPA::eps records the worst eps at each
        (E, q), and n_solves, n_failed, n_interpolated the counts.

  @note Parallel over chains, each one (E, K, operator) and a run of
        consecutive q (warm starts along it). The energies are taken in
        batches of the thread count, the continuum states of a batch held
        together, and the chains of a batch are one dynamically scheduled
        pool; the run of q is chosen to give a few chains per thread.

  @warning The continuum states must be those of the residual ion
           (@p hole_particle = true) for the RPA amplitude to be consistent;
           see ExternalField::TDHFcntm::dV_complex.
*/
FormFactorsRPA solve_formFactors_RPA(
  const HF::HartreeFock *vHF,
  const std::optional<std::array<int, 2>> &lc_minmax, double ec_min,
  bool force_rescale, bool hole_particle, bool force_orthog,
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
  for every multipole over the full grids. The RPA is solved
  (solve_formFactors_RPA()) only for the multipoles up to RPAOptions::K_max,
  on the region of the grids with \f$ E \le E_{\rm max} \f$ and
  \f$ q \le q_{\rm max} \f$ (RPAOptions::E_max, q_max); the RPA factors
  are the bare ones plus the RPA shift of that block (its rpa - bare) on the
  region, and the bare ones elsewhere. In the diagonal case (q = E/c) the q
  limit is a limit on E.

  With a stride (RPAOptions::E_stride, q_stride) the RPA is solved only at
  some points of the region: every stride-th energy and momentum transfer
  (stride_points()), the last of each, and the first and last energy at which
  each orbital is ionised (so that no factor is interpolated across its
  threshold). Failed solves are corrected on those points, as without a
  stride. The RPA correction is then interpolated to every point of the
  region (interpolate_rpa()). For the direct factors it is the
  relative shift \f$ R = K_{\rm RPA}/K_{\rm bare} \f$ of the RPA
  multipoles that is interpolated, since it varies slowly where the factors
  themselves change by orders of magnitude,
  \f[
    K_{\rm RPA}(E,q) = K_{\rm bare}(E,q) + [R(E,q) - 1]\,
                       K^{K \le K_{\rm max}}_{\rm bare}(E,q);
  \f]
  for the cross terms (X, Y, Z), which may change sign, it is the shift
  \f$ K_{\rm RPA} - K_{\rm bare} \f$ itself. The bare factors are
  calculated in two parts, the RPA multipoles and the higher ones, so that
  those of the RPA multipoles alone are known at every point.

  Prints "Calculating RPA form factors:" with a progress bar over the RPA
  solves, then "RPA-excluded contributions:" with a progress bar over each
  part of the bare calculation.

  The limits exist because the RPA is expensive, and hard to converge, at
  high E, q, and K (many open channels; rapidly oscillating operators) where
  its effect on the factors is small: they confine the work to the low-E,
  low-q, low-K region where the core polarisation matters.

  The RPA block may be solved on its own radial grid (@p vHF_rpa,
  @p jK_tab_rpa): the RPA correction is formed from the rpa and bare factors
  of the block itself, both on that grid, so the grid need only be dense
  enough for the E and q within the RPA limits, not for the full grids.

  @param vHF            Hartree-Fock potential; its core defines the orbitals.
  @param lc_minmax      Optional limits on the continuum orbital l.
  @param ec_min         Minimum ejected electron energy, in au.
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
  @param vHF_rpa        Hartree-Fock potential (same atom) on the radial grid
                        the RPA block is solved on; nullptr: @p vHF.
  @param jK_tab_rpa     Spherical Bessel table on the grid of @p vHF_rpa; must
                        hold every q within the RPA limits. nullptr:
                        @p jK_tab.
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
  bool force_rescale, bool hole_particle, bool force_orthog,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, const RPAOptions &rpa_options = {},
  const HF::HartreeFock *vHF_rpa = nullptr,
  const SphericalBessel::JL_table *jK_tab_rpa = nullptr);

//------------------------------------------------------------------------------
/*!
  @brief Adds the RPA shift, solved at some points of a region of the grids,
  to the form factors at every point of the region.

  @details
  The RPA factors of @p factors must hold the bare ones on entry. If every
  point of the region was solved, the shift (rpa - bare of @p solved) is
  added as it is. Otherwise the RPA correction is interpolated from the
  solved points, bilinearly in ln E and ln q (interpolate_log_bilinear()):
  for the direct factors the relative shift
  \f$ R = K_{\rm RPA}/K_{\rm bare} \f$, which is then applied to the bare
  factor of the same multipoles at each point (@p bare_rpa_K); for the cross
  terms (X, Y, Z), which may change sign, the shift itself. See
  calculate_formFactors_RPA().

  Also sets the eps of the solved points (it is left as it is at the others),
  and copies the counts of solves.

  @param solved      Bare and RPA factors at the solved points only, as from
                     solve_formFactors_RPA(): indexed (E_solved, q_solved).
  @param bare_rpa_K  Bare factors of the same multipoles, over the full grids.
  @param Egrid       Full energy transfer grid.
  @param qgrid       Full momentum transfer grid.
  @param E_region    The energies of the region: indices into @p Egrid, in
                     ascending order of energy.
  @param q_region    As @p E_region, for the momentum transfer.
  @param E_solved    The solved energies: indices into @p E_region,
                     ascending; first and last points of the region included.
  @param q_solved    As @p E_solved, for the momentum transfer.
  @param factors     Factors over the full grids; rpa, eps and the counts are
                     updated.
*/
void interpolate_rpa(const FormFactorsRPA &solved,
                     const std::vector<FormFactorSet> &bare_rpa_K,
                     const std::vector<double> &Egrid,
                     const std::vector<double> &qgrid,
                     const std::vector<std::size_t> &E_region,
                     const std::vector<std::size_t> &q_region,
                     const std::vector<std::size_t> &E_solved,
                     const std::vector<std::size_t> &q_solved,
                     FormFactorsRPA *factors);

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
  @brief Corrects the (E, q) points at which the RPA solve failed, by
  interpolating the relative RPA shift from the neighbours in E and q.

  @details
  As interpolate_failed_rpa(), with the neighbours on both axes: the
  adjacent rows (E) and columns (q). The shift \f$ R = K_{\rm RPA}/K_{\rm
  bare} \f$ is the mean over the converged neighbours with a non-zero bare
  value (up to four); with a single one, half its relative correction. A
  point with none is left as it is. Failed points are never used as
  neighbours. With a single column (the diagonal, q = E/c, case) the
  neighbours are along E.

  @param eps       Final RPA eps of each point (E rows x q columns); a point
                   that was not solved must not count as failed (eps of 0,
                   or negative).
  @param eps_fail  Failure threshold (see rpa_failed()).
  @param K_bare    Bare (no-RPA) values, same shape as @p eps.
  @param K_rpa     RPA values, same shape; corrected in place.
*/
void interpolate_failed_rpa_2d(const LinAlg::Matrix<double> &eps,
                               double eps_fail,
                               const LinAlg::Matrix<double> &K_bare,
                               LinAlg::Matrix<double> *K_rpa);

//------------------------------------------------------------------------------
/*!
  @brief The operators (indices into ChannelAmplitudes) that a form factor is
  built from.

  @details
  The direct factors depend on one amplitude, the cross terms on two (X: t
  and L; Y: t5 and L5) or four (Z: E, M, E5, M5); see
  accumulate_formFactors(). A factor is affected by a failed RPA solve of
  any of its operators.

  @param i_factor  Index into FormFactorSet.
  @return The operator indices.
*/
std::vector<std::size_t> factor_operators(std::size_t i_factor);

//------------------------------------------------------------------------------
/*!
  @brief The points of a grid kept by a stride: every stride-th point, the
  last point, and any points that must be kept.

  @param num_points  Number of points of the grid.
  @param stride      Every stride-th point is kept, from the first (1, or 0:
                     every point).
  @param keep        Points that are kept whatever the stride.
  @return The kept points, in ascending order, without repeats.
*/
std::vector<std::size_t>
stride_points(std::size_t num_points, std::size_t stride,
              const std::vector<std::size_t> &keep = {});

//------------------------------------------------------------------------------
/*!
  @brief The smallest grid value at or above a limit: the grid point that
  covers the limit.

  @details
  Used to extend a limit given as input (e.g., the RPA E and q limits) to a
  point of the grid, so that the range up to the limit is covered in full.

  @param grid   Grid values (any order).
  @param limit  The limit.
  @return The smallest value of @p grid that is >= @p limit; @p limit itself
          if every value is below it.
*/
double grid_point_at_or_above(const std::vector<double> &grid, double limit);

//------------------------------------------------------------------------------
/*!
  @brief Interpolates a table known on some rows and columns of a grid to
  every row and column: bilinear in the logarithms of the grid values.

  @details
  Between the known rows \f$ x_0 \le x \le x_1 \f$ and columns
  \f$ y_0 \le y \le y_1 \f$:
  \f[
    f(x,y) = (1-u)(1-v) f_{00} + (1-u) v f_{01} + u (1-v) f_{10} + u v f_{11},
  \f]
  \f[
    u = \frac{\ln(x/x_0)}{\ln(x_1/x_0)},
  \f]
  and v likewise in y. A known point is returned unchanged. Outside the range
  of the known rows (columns) the value of the nearest one is used; with a
  single known row (column) the table is constant in that direction.

  @param x_known  Grid values of the known rows: positive, ascending.
  @param y_known  Grid values of the known columns: positive, ascending.
  @param f_known  The table on the known rows and columns.
  @param x        Grid values of every row (positive).
  @param y        Grid values of every column (positive).
  @return The table on every row and column.
*/
LinAlg::Matrix<double> interpolate_log_bilinear(
  const std::vector<double> &x_known, const std::vector<double> &y_known,
  const LinAlg::Matrix<double> &f_known, const std::vector<double> &x,
  const std::vector<double> &y);

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
  @param radial_grid   Parameters of the radial grid, as text for the header
                       ('r grid' line; left out if empty).
  @param rpa_E_grid    Energy transfers at which the RPA was solved, in au,
                       ascending ('RPA E grid' line; left out if empty).
  @param rpa_q_grid    As @p rpa_E_grid, for the momentum transfers ('RPA q
                       grid' line; left out if empty, or if @p diagonal).
*/
void write_to_file_xyz(const std::string &filename,
                       const std::vector<double> &E_grid,
                       const std::vector<double> &q_grid,
                       const std::vector<std::string> &titles,
                       const std::vector<std::string> &descriptions,
                       std::vector<LinAlg::Matrix_view<const double>> factors,
                       Units units = Units::Particle, int num_digits = 6,
                       bool diagonal = false,
                       const std::string &radial_grid = "",
                       const std::vector<double> &rpa_E_grid = {},
                       const std::vector<double> &rpa_q_grid = {});

//------------------------------------------------------------------------------
//! As write_to_file_xyz, for a FormFactorSet (empty factors are skipped)
void write_to_file_xyz_13(
  const std::string &filename, const std::vector<double> &E_grid,
  const std::vector<double> &q_grid, const std::vector<std::string> &titles,
  const std::vector<std::string> &descriptions, const FormFactorSet &K_factors,
  Units units = Units::Particle, int num_digits = 6, bool diagonal = false,
  const std::string &radial_grid = "",
  const std::vector<double> &rpa_E_grid = {},
  const std::vector<double> &rpa_q_grid = {});

} // namespace Kion
