#pragma once
#include "Kionisation/Kion_functions.hpp"
#include "Maths/SphericalBessel.hpp"
#include <array>
#include <cstddef>
#include <vector>
class DiracSpinor;
namespace HF {
class HartreeFock;
}

//! Bethe-ridge (plane-wave) completion of the truncated multipole sum
namespace Kion {

//------------------------------------------------------------------------------
/*!
  @brief Momentum-space (Fourier transformed) bound orbital, for the all-K
  plane-wave response.

  @details
  The momentum-space orbital,
  \f[
    \tilde\psi_{am}(\vb{p}) = \int d^3r\, e^{-i\vb{p}\cdot\vb{r}}\,
      \psi_{am}(\vb{r})
    = (-i)^l \begin{pmatrix}
        \tilde f_a(p)\, \ket{\kappa m} \\
        -s_\kappa\, \tilde g_a(p)\, \ket{-\kappa, m}
      \end{pmatrix}
  \f]
  (the spinors are functions of \f$ \hat{\vb{p}} \f$, and
  \f$ s_\kappa = \kappa/|\kappa| \f$), has the radial functions
  \f[
  \begin{align}
    \tilde f_a(p) &= 4\pi\int f_a(r)\, j_l(pr)\, r\, dr, \\
    \tilde g_a(p) &= 4\pi\int g_a(r)\, j_{\tilde l}(pr)\, r\, dr,
  \end{align}
  \f]
  with \f$ \tilde l = l(-\kappa) \f$, normalised as
  \f[
    \int (\tilde f_a^2 + \tilde g_a^2)\, p^2\, dp = (2\pi)^3.
  \f]
  The grid is logarithmic, and ends where the density
  \f$ \tilde f_a^2 + \tilde g_a^2 \f$ has fallen below a set fraction of
  its peak (see momentum_orbital()): the density is zero beyond p.back().
*/
struct MomentumOrbital {
  //! Dirac quantum number of the orbital
  int kappa{};
  //! Binding energy (au); that of the real (HF) orbital
  double en{};
  //! Number of electrons in the orbital: (2j+1) times the occupation fraction
  double num_electrons{};
  //! Momentum grid (au), logarithmic
  std::vector<double> p{};
  //! Transformed large component, f~_a, at each p
  std::vector<double> F{};
  //! Transformed small component, g~_a, at each p
  std::vector<double> G{};
};

//------------------------------------------------------------------------------
/*!
  @brief Fourier transforms a bound orbital onto a logarithmic momentum grid
  (MomentumOrbital).

  @details
  The grid is logarithmic, @p points_per_decade points per factor of 10 in
  p, running upwards from @p p_min, and stops (at most at @p p_max) once the density has fallen below
  @p density_cut times its peak, beyond the peak. The radial integrals run
  over the extent of the orbital (its max_pt()), with the cell-averaged
  Bessel functions of SphericalBessel::fillBesselVec_kr, so the high-p tail
  (from small r) is integrated correctly on a coarse grid.

  @param Fa                 Bound orbital (as used in the matrix elements);
                            its occupation sets the electron count.
  @param en                 Binding energy to record (au): that of the real
                            orbital, which may differ from Fa.en() for a
                            Zeff model state.
  @param p_min              First momentum (au).
  @param p_max              Largest momentum considered (au).
  @param points_per_decade  Grid points per factor of 10 in p.
  @param density_cut        Stop once density < density_cut times the peak.
  @return The momentum-space orbital.
*/
MomentumOrbital momentum_orbital(const DiracSpinor &Fa, double en,
                                 double p_min = 1.0e-3, double p_max = 1.0e4,
                                 std::size_t points_per_decade = 150,
                                 double density_cut = 1.0e-12);

//------------------------------------------------------------------------------
/*!
  @brief Integrand of the all-K plane-wave response at bound-electron
  momentum p: the Dirac traces of every factor.

  @details
  With z along q, \f$ \vb{p}_f = \vb{p} + \vb{q} \f$, and the angle between
  p and q fixed by energy conservation,
  \f[
    \cos\theta_{pq} = \frac{p_f^2 - p^2 - q^2}{2pq},
  \f]
  returns the traces over the Dirac indices
  \f[
    T = {\rm Tr}\big[(E_f + c\,\vb{\alpha}\cdot\vb{p}_f + \beta mc^2)\,
        \Gamma\,\rho\,\Gamma'^\dagger\big]
  \f]
  for the momentum-space density matrix of the shell without its
  \f$ N_a/(8\pi) \f$ prefactor (2x2 blocks; see planewave_formFactors()),
  \f[
    \rho = \begin{pmatrix}
      \tilde f_a^2 & s_\kappa\,\tilde f_a\tilde g_a\,\vb{\sigma}\cdot\hat{\vb{p}} \\
      s_\kappa\,\tilde f_a\tilde g_a\,\vb{\sigma}\cdot\hat{\vb{p}} & \tilde g_a^2
    \end{pmatrix},
  \f]
  where \f$ s_\kappa = \kappa/|\kappa| \f$. The traces are returned in the
  order {V: 00, 33, 11+22, 03; A: 00, 33, 11+22, 03; Im VA 12; S; P}.
  See planewave_formFactors() for the vertices and the use.

  @param p      Bound-electron momentum (au).
  @param F      Transformed large component, f~_a, at p.
  @param G      Transformed small component, g~_a, at p.
  @param kappa  Dirac quantum number of the orbital.
  @param pf     Ejected-electron momentum, p_f (au).
  @param q      Momentum transfer (au).
  @param ef     Ejected-electron kinetic energy (au).
  @param alpha  Fine-structure constant.
  @return The 11 traces.
*/
std::array<double, 11> planewave_traces(double p, double F, double G, int kappa,
                                        double pf, double q, double ef,
                                        double alpha);

//------------------------------------------------------------------------------
/*!
  @brief The 13 form factors of one orbital with a free (plane-wave)
  ejected electron, summed over all multipoles, at energy transfer E and
  momentum transfer q.

  @details
  The ejected electron is free (the potential does not act during the
  collision), with kinetic energy, momentum, and total energy
  \f[
  \begin{align}
    \en_f &= E + \en_a, \\
    p_f &= \sqrt{\en_f(2 + \alpha^2\en_f)}, \\
    E_f &= \en_f + mc^2,
  \end{align}
  \f]
  and the bound electron has the momentum distribution of the orbital
  (MomentumOrbital). The factors of the orbital are the Cartesian
  components (z along q) of its response tensor,
  \f[
    R^{\mu\nu}_a = x_a \sum_{m_a}\sum_f J^\mu_{fa}\, J^{\nu *}_{fa},
  \f]
  summed over the final plane waves at energy \f$ \en_f \f$ and over the
  closed shell. The sum reduces to a single integral over the
  bound-electron momentum p, the angle between p and q being fixed by
  energy conservation (planewave_traces()):
  \f[
    R^{\mu\nu}_a = \frac{1}{8\pi^2 c^2 q}\int_{|p_f-q|}^{p_f+q} p\,dp\;
      {\rm Tr}\big[(E_f + c\,\vb{\alpha}\cdot\vb{p}_f + \beta mc^2)\,
      \Gamma^\mu \rho_a(p)\, \Gamma^{\nu\dagger}\big],
  \f]
  with \f$ \vb{p}_f = \vb{p} + \vb{q} \f$, the trace over the Dirac
  indices, and
  \f[
    \rho_a = \frac{N_a}{8\pi}\begin{pmatrix}
      \tilde f_a^2 & s_\kappa\,\tilde f_a\tilde g_a\,\vb{\sigma}\cdot\hat{\vb{p}} \\
      s_\kappa\,\tilde f_a\tilde g_a\,\vb{\sigma}\cdot\hat{\vb{p}} & \tilde g_a^2
    \end{pmatrix}
  \f]
  the momentum-space density matrix of the shell (\f$ N_a \f$ electrons).
  The vertex \f$ \Gamma^\mu = \gamma^0\gamma^\mu\tilde\gamma \f$ is that
  of the operator:
  - vector: \f$ (1, \vb{\alpha}) \f$
  - axial: \f$ (\gamma^5, \vb{\Sigma}) \f$
  - scalar: \f$ \beta \f$
  - pseudoscalar: \f$ i\beta\gamma^5 \f$

  The factors are the components
  \f[
  \begin{align}
    V_T &= R^{00}, \\
    V_L &= R^{33}, \\
    V_E + V_M &= R^{11} + R^{22}, \\
    X &= -{\rm Re}\,R^{03}, \\
    A_T &= R^{00}_A, \\
    A_L &= R^{33}_A, \\
    A_E + A_M &= R^{11}_A + R^{22}_A, \\
    Y &= {\rm Re}\,R^{03}_A, \\
    Z &= 2\,{\rm Im}\,R^{12}_{VA}, \\
    S &= R_{SS}, \\
    P &= R_{PP},
  \end{align}
  \f]
  where Z is the antisymmetric transverse vector-axial interference. With
  the couplings in the vertex, \f$ \gamma^0\gamma^\mu(c_V - c_A\gamma^5) \f$,
  the tensor is \f$ c_V^2 R_V + c_A^2 R_A - c_V c_A (R_{VA} + R_{AV}) \f$,
  and its components are the coefficients of the spin-summed cross
  section:
  \f[
  \begin{align}
    R^{00} &= c_V^2 V_T + c_A^2 A_T, \\
    R^{33} &= c_V^2 V_L + c_A^2 A_L, \\
    R^{11} + R^{22} &= c_V^2 (V_E + V_M) + c_A^2 (A_E + A_M), \\
    -{\rm Re}\,R^{03} &= c_V^2 X - c_A^2 Y, \\
    -{\rm Im}\,R^{12} &= c_V c_A Z,
  \end{align}
  \f]
  (for an isotropic target \f$ R^{21}_{VA} = -R^{12}_{VA} \f$, so the VA and
  AV terms contribute \f$ {\rm Im}\,R^{12}_{VA} \f$ each), and
  \f$ R = c_S^2 S + c_P^2 P \f$ for the scalar-pseudoscalar vertex
  \f$ c_S\gamma^0 + i c_P\gamma^0\gamma^5 \f$. For an electron at rest the
  factors reduce to the free-electron values, e.g.
  \f[
  \begin{align}
    V_L &= (E/qc)^2\, V_T, \\
    X &= -(E/qc)\, V_T, \\
    A_L &= V_T.
  \end{align}
  \f]

  The electric and magnetic multipoles are not separately defined by the
  Cartesian tensor (only their sum enters a spin-summed cross section): the
  whole transverse response is returned as the electric factor, and the
  magnetic factors are zero.

  @param orb    Momentum-space orbital (momentum_orbital()).
  @param E      Energy transfer (au).
  @param q      Momentum transfer (au).
  @param alpha  Fine-structure constant.
  @return The factors, in FormFactorSet order; all zero if the orbital is not
          ionised (\f$ \en_f \le 0 \f$).
*/
std::array<double, 13> planewave_formFactors(const MomentumOrbital &orb,
                                             double E, double q, double alpha);

//------------------------------------------------------------------------------
/*!
  @brief Fraction of the norm of \f$ e^{i\vb{q}\cdot\vb{r}}\psi_a \f$ carried by
  the multipoles K <= Kmax.

  @details
  \f[
    c_a(q) = \sum_{K=0}^{K_{\rm max}} (2K+1)\int (f_a^2 + g_a^2)\,
             j_K(qr)^2\, dr,
  \f]
  which tends to 1 by the identity
  \f[
    \sum_{K=0}^{\infty} (2K+1)\, j_K(x)^2 = 1.
  \f]
  The remainder 1 - c is the norm of the part of the state that the
  multipoles above Kmax carry, and bounds what their omission can cost: the
  plane-wave completion is skipped where it is negligible.

  @param Fa      Bound orbital.
  @param Kmax    Maximum multipolarity of the sum.
  @param iq      Momentum-transfer index in the Bessel table.
  @param jK_tab  Spherical Bessel table, on the orbital's grid.
  @return The captured fraction.
*/
double captured_fraction(const DiracSpinor &Fa, int Kmax, std::size_t iq,
                         const SphericalBessel::JL_table &jK_tab);

//------------------------------------------------------------------------------
/*!
  @brief The Bethe-ridge correction: completes the multipole sum of every core
  orbital to all K with plane waves.

  @details
  The multipole sum truncated at Kmax misses the response near the
  quasi-free (Bethe) ridge,
  \f[
    E \approx \frac{q^2}{2m},
  \f]
  where multipoles up to \f$ K \sim q\,r_a \f$ contribute. There the ejected
  electron is fast and its high partial waves are undistorted, so the
  missing part is supplied by plane waves. Per orbital and factor, the
  returned correction is
  \f[
    \Delta K_a = K_a^{\rm PW}({\rm all}\ K) - K_a^{\rm PW}(K \le K_{\rm max}),
  \f]
  the all-K plane-wave response (planewave_formFactors()) minus the same
  multipole sum evaluated with free Dirac spherical waves (DiracODE::freeDirac()), i.e.
  exactly the plane-wave multipoles above Kmax. Added to the distorted-wave
  factors of calculate_formFactors(), the total is
  \f[
    K = K^{\rm DW}(K \le K_{\rm max}) + \Delta K.
  \f]
  Off the ridge the
  correction vanishes by itself; on the ridge the result is independent of
  Kmax once the low multipoles, which carry the distortion, are included.

  The plane-wave sum always runs K = 0..Kmax over the full continuum l
  range, with the same operators, Bessel table, bound states, and
  continuum truncation (max_pt) as the distorted-wave sum, so that the
  K <= Kmax parts cancel to grid accuracy. It has no hole-particle or
  orthogonalisation (distorted-wave physics only). Bound-bound (Pauli)
  leakage of the plane waves is confined to K <= j_a + j_b and cancels as
  long as Kmax >= 2 j_max of the core (checked). The squared factors are
  clamped at zero (they are sums of squares; a negative value is noise from
  the subtraction). The whole transverse correction goes into the electric
  factors (see planewave_formFactors()).

  The correction is computed only where it can matter: for each orbital,
  at the momentum transfers where the multipoles above Kmax carry more than
  @p ridge_eps of the norm (captured_fraction()), at the energies where
  the orbital is ionised (ejected energy within ( @p ec_min, @p ec_max ),
  and, at each energy, only where |p_f - q| lies within the momentum grid
  of the orbital (elsewhere, far off the ridge, the plane-wave response is
  zero by construction).
  Parallel over the (orbital, E) pairs. Prints, per orbital, the momentum
  transfer above which the correction is active.

  @param vHF            Hartree-Fock; its core defines the orbitals.
  @param bound_states   Bound state of each orbital as used in the matrix
                        elements (model_bound_states()); must be the same
                        as used for the distorted-wave factors.
  @param ec_min         Minimum ejected electron energy, in au.
  @param ec_max         Maximum ejected electron energy, in au.
  @param Egrid          Energy transfer grid, in au.
  @param qgrid          Momentum transfer grid, in au.
  @param jK_tab         Spherical Bessel table on @p qgrid (as used for the
                        distorted-wave factors).
  @param Kmax           Maximum multipolarity of the distorted-wave sum.
  @param vectorQ        Vector factors.
  @param axialQ         Axial-vector factors.
  @param scalarQ        Scalar factor.
  @param pseudoscalarQ  Pseudoscalar factor.
  @param spatialQ       Spatial (E, M, L) components.
  @param ridge_eps      Skip (orbital, q) where the multipoles above Kmax
                        carry less than this fraction of the norm.
  @return The correction of each core orbital, indexed as the core, with the
          allocation pattern of allocate_formFactors().

  @warning Requires Kmax >= 2 j_max(core): if not, prints a warning and
           returns a zero correction. Not defined for the diagonal (q = E/c)
           case, which never approaches the ridge.
*/
std::vector<FormFactorSet> calculate_ridge_correction(
  const HF::HartreeFock *vHF, const std::vector<DiracSpinor> &bound_states,
  double ec_min, double ec_max, const std::vector<double> &Egrid,
  const std::vector<double> &qgrid, const SphericalBessel::JL_table &jK_tab,
  int Kmax, bool vectorQ, bool axialQ, bool scalarQ, bool pseudoscalarQ,
  bool spatialQ, double ridge_eps);

} // namespace Kion
