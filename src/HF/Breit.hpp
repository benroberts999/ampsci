#pragma once
#include "Coulomb/meTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <cassert>
#include <utility>
#include <vector>

namespace HF {

//==============================================================================
namespace Breit_gb {

struct single_k_mop {
  // Class to hold the Breit-Coulomb integrals, for single k
public:
  single_k_mop() {}
  single_k_mop(const DiracSpinor &Fi, const DiracSpinor &Fj, int k) {
    calculate(Fi, Fj, k);
  }

  void calculate(const DiracSpinor &Fi, const DiracSpinor &Fj, int k);
  std::vector<double> b0_minus{}, bi_minus{};
  std::vector<double> b0_plus{}, bi_plus{};
  std::vector<double> g0_minus{}, gi_minus{};
  std::vector<double> g0_plus{}, gi_plus{};
};

struct single_k_n {
  // Class to hold the Breit-Coulomb integrals, for single k
public:
  single_k_n() {}
  single_k_n(const DiracSpinor &Fi, const DiracSpinor &Fj, int k) {
    calculate(Fi, Fj, k);
  }

  void calculate(const DiracSpinor &Fi, const DiracSpinor &Fj, int k);
  std::vector<double> g{};

private:
  std::vector<double> gi{};
};

} // namespace Breit_gb

namespace Breit_gh_freqdep {

struct single_k_mop_freq {
  // Class to hold the frequency-dependent Breit-Coulomb integrals, for single k
public:
  single_k_mop_freq() {}
  single_k_mop_freq(const DiracSpinor &Fi, const DiracSpinor &Fj, int k,
                    const double w) {
    calculate(Fi, Fj, k, w);
  }

  void calculate(const DiracSpinor &Fi, const DiracSpinor &Fj, int k,
                 const double w);
  std::vector<double> g0_minus_freqw{}, gi_minus_freqw{};
  std::vector<double> g0_plus_freqw{}, gi_plus_freqw{};
  std::vector<double> h0_minus_freqw{}, hi_minus_freqw{};
  std::vector<double> h0_plus_freqw{}, hi_plus_freqw{};
  std::vector<double> v1_freqw{}, v2_freqw{}, v3_freqw{}, v4_freqw{};
};

struct single_k_n_freq {
  // Class to hold the Breit-Coulomb integrals, for single k
public:
  single_k_n_freq() {}
  single_k_n_freq(const DiracSpinor &Fi, const DiracSpinor &Fj, int k,
                  const double w) {
    calculate(Fi, Fj, k, w);
  }

  void calculate(const DiracSpinor &Fi, const DiracSpinor &Fj, int k,
                 const double w);
  std::vector<double> g{};

private:
  std::vector<double> gi{};
};

} // namespace Breit_gh_freqdep

//==============================================================================
//! Breit potentials for one- (Hartree-Fock Breit) and two-body Breit integrals
class Breit {
  // overall scaling factor
  double m_scale;

  // Scaling factor for frequency in frequency-dependent Breit (default 0.0 = static)
  double m_lambda_f;

  // Scaling factors for each term (mainly for tests). {M,N}=Gaunt; {O,P} retarded
  double m_M{1.0}, m_N{1.0}, m_O{1.0}, m_P{1.0};

  // For speedy lookup, when only full integral is required
  std::vector<Coulomb::meTable<Breit_gb::single_k_mop>> m_gb{};
  std::vector<Coulomb::meTable<Breit_gb::single_k_n>> m_gb_N{};

public:
  /*!
    @brief Constructs Breit interaction operator

    @details
    Initializes the Breit operator with an overall scaling factor. Individual term
    scaling factors (M, N, O, P) default to 1.0.
    Also scaling for frequency dependence (lambda_f); default to 0.0.

    @param scale    Overall scaling factor for Breit contributions (default 1.0)
    @param lambda_f Scaling factor for frequency-dependent Breit (default 0.0)

    @note If lambda_f is zero (default), frequency-independent form will be used.
    If 1, then will be frequency-dependent.
    Can be any value (e.g., 0.5 to test linearity/scaling of f-dependent results).
  */
  Breit(double scale = 1.0, double lambda_f = 0.0)
    : m_scale(scale), m_lambda_f(lambda_f) {}

  /*!
    @brief Update all scaling factors

    @details
    Updates the overall scaling factor, individual term scaling factors (M, N, O, P),
    and the frequency scaling factor for frequency-dependent Breit calculations.

    @param t_scale     Overall scaling factor (default 1.0)
    @param t_M        Scaling for M term (Gaunt part, default 1.0)
    @param t_N        Scaling for N term (Gaunt part, default 1.0)
    @param t_O        Scaling for O term (retarded part, default 1.0)
    @param t_P        Scaling for P term (retarded part, default 1.0)

    @note Does not update lambda_f (f-dependent scaling)
  */
  void update_scale(double t_scale = 1.0, double t_M = 1.0, double t_N = 1.0,
                    double t_O = 1.0, double t_P = 1.0) {
    m_scale = t_scale;
    m_M = t_M;
    m_N = t_N;
    m_O = t_O;
    m_P = t_P;
  }

  /*!
    @brief Update frequency scaling factor

    @details
    Sets the scaling factor for the frequency in frequency-dependent Breit calculations.
    The frequency used in the integrals is multiplied by this factor.

    @param lambda_f Frequency scaling factor (should be > 0; default 1.0)

    @note Setting lambda_f to very small values approaches the static Breit limit.
  */
  void update_lambda_f(double lambda_f) { m_lambda_f = lambda_f; }

  /*!
    @brief Selection rule check for Breit integrals

    @details
    Tests whether the four-body Breit integral B^k_{vwxy} is nonzero based on
    angular momentum selection rules.

    @param k Multipolarity (angular momentum rank of the interaction)
    @param v (w,x,y): electron states

    @return True if the integral is potentially nonzero (angular momentum
            selection rules are satisfied); false if the integral vanishes.

  */
  static bool Bk_SR(int k, const DiracSpinor &v, const DiracSpinor &w,
                    const DiracSpinor &x, const DiracSpinor &y) {

    const auto have_mop = Angular::Ck_kk_SR(k, v.kappa(), x.kappa()) &&
                          Angular::Ck_kk_SR(k, w.kappa(), y.kappa()) &&
                          v != x && w != y;

    const auto have_n = Angular::Ck_kk_SR(k, -v.kappa(), x.kappa()) &&
                        Angular::Ck_kk_SR(k, -w.kappa(), y.kappa()) &&
                        (v.kappa() + x.kappa() != 0) &&
                        (w.kappa() + y.kappa() != 0);

    return (have_mop || have_n);
  }

  /*!
    @brief Determine valid multipolarity range for Breit integrals

    @details
    Minimum and maximum allowed multipolarity k for a four-body
    Breit integral.

    @param a (b,c,d) electron states

    @return A pair {k_min, k_max} giving the valid multipolarity range.
            k_min is at least 1 (forbidden to have k=0 for HF Breit).
            Returns {k_max, k_min} (invalid range) if no valid k exists.

    @note This is the static method analogue of Coulomb::k_minmax_tj, adapted
          for the four-body Breit structure.
    
    @warning Discounts k=0 - only valid within Hartree-Fock??
  */
  static std::pair<int, int> k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b,
                                      const DiracSpinor &c,
                                      const DiracSpinor &d) {
    const auto [k1, k2] = Coulomb::k_minmax_tj(a.twoj(), c.twoj());
    const auto [k3, k4] = Coulomb::k_minmax_tj(b.twoj(), d.twoj());
    return {std::max({1, k1, k3}), std::min(k2, k4)};
  }

  /*!
    @brief Determine valid multipolarity range for Breit integrals from quantum numbers

    @details
    As @ref k_minmax but for 2*j (doesn't require DiracSpinors)
  */
  static std::pair<int, int> k_minmax_tj(int tja, int tjb, int tjc, int tjd) {
    const auto [k1, k2] = Coulomb::k_minmax_tj(tja, tjc);
    const auto [k3, k4] = Coulomb::k_minmax_tj(tjb, tjd);
    return {std::max({1, k1, k3}), std::min(k2, k4)};
  }

  /*!
    @brief Precompute Breit integral lookup tables for rapid evaluation

    @details
    Pre-calculates and stores reduced Breit integrals for all unique pairs of
    basis orbitals and all multipolarity values k, enabling
    much faster integral lookups via the Bk_abcd_2() family of functions.
    This trades memory for speed, allowing rapid evaluation in HF iterations.

    @param basis    The set of basis orbitals to use
    @param t_max_k  Maximum multipolarity k to compute (default 99).
                    Actual maximum used is the minimum of this value and
                    the physical constraint k_minmax() for the given basis.

    @warning This function uses substantial memory — use with caution for large 
    calculations. Have mostly found the speedup is not worth the memory cost, so 
    this is generally not used anymore.

    @note Must be called once before using the faster variants Bk_abcd_2(),
          BPk_abcd_2(), or BWk_abcd_2(). Calling it multiple times will
          recompute and overwrite previous results.
    
    @note ONLY for static Breit
  */
  void fill_gb(const std::vector<DiracSpinor> &basis, int t_max_k = 99);

  //! Returns the overall scaling factor
  double scale_factor() const { return m_scale; };

  /*!
    @brief Calculates Breit contribution with automatic frequency dependence

    @details
    Computes the Hartree-Fock Breit interaction V_br*Fa for a valence electron
    interacting with the core. This is the direct Breit contribution from all
    core electrons.

    Will be static version is lambda = 0, otherwise, frequency-dependent

    @param Fa   Valence electron state
    @param core Core electron states

    @return The Breit-Hartree-Fock potential applied to Fa, computed at the
            appropriate frequency regime based on m_lambda_f.

    @note For frequency-dependent calculations, this may be substantially slower
          due to spherical Bessel function evaluation in the radial integrals.
  */
  DiracSpinor VbrFa(const DiracSpinor &Fa,
                    const std::vector<DiracSpinor> &core) const;

  /*!
    @brief Calculates Breit potential (wrapper, equivalent to VbrFa)

    @details
    Backward compatibility wrapper. Now simply calls VbrFa(), which automatically
    handles both static and frequency-dependent regimes based on m_lambda_f.

    @param Fa   Valence state
    @param core Core electron states

    @return Equivalent to VbrFa(Fa, core).

    @note Kept for backward compatibility. New code should use VbrFa() directly.
  */
  [[deprecated]] DiracSpinor
  VbrFa_freqw(const DiracSpinor &Fa,
              const std::vector<DiracSpinor> &core) const;

  /*!
    @brief Breit-TDHF: Breit correction to the TDHF correction to Hartree-Fock

    @details
    Calculates the Breit correction to the TDHF potential, dV.
    This represents the response of the Breit field to the core electron perturbations.

    @param kappa Dirac quantum number of the resulting state/projection
    @param K     Multipolarity (rank) of the RPA operator
    @param Fa    Electron state (acting on this)
    @param Fb    Core state undergoing perturbation
    @param Xbeta X perturbation to core state: from @ref ExternalField::TDHF
    @param Ybeta Y perturbation to core state: from @ref ExternalField::TDHF

    @return The reduced RPA correction dV_Br*Fa
  */
  DiracSpinor dV_Br(int kappa, int K, const DiracSpinor &Fa,
                    const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                    const DiracSpinor &Ybeta) const;

  /*!
    @brief Reduced Breit two-body matrix element with automatic frequency dependence

    @details
    Calculates the reduced two-body Breit matrix element B^k_{abcd}, the Breit
    analogue of the Coulomb Q^k integral. The Breit interaction is the sum of
    Gaunt (instantaneous magnetic) and retarded (transverse photon) contributions,
    decomposed into multipole ranks k.

    Internally checks m_lambda_f: if zero (or close to zero), uses the static
    (frequency-independent) formulation. Otherwise, uses frequency-dependent
    formulation with frequency w = lambda_f * alpha * |E_a - E_c|.

    The individual term scaling factors (m_M, m_N, m_O, m_P) and overall
    scaling m_scale modulate the contributions:
    - M, N terms arise from Gaunt (instantaneous) interaction
    - O, P terms arise from the retarded (photon propagation) contribution

    @param k  Multipolarity (angular momentum rank of the interaction)
    @param Fa (Fb,Fc,Fd) electron states

    @return The reduced Breit matrix element B^k_{abcd}, evaluated at the
            appropriate frequency based on m_lambda_f.

    @note Selection rules depend on angular momentum quantum numbers via Bk_SR().
          Use k_minmax() to determine the valid multipolarity range before calling.

    @note This function re-computes all four terms (m, n, o, p) each call.
          For repeated calls with the same basis, call fill_gb() once, then use
          the faster variant @ref Bk_abcd_2() instead.
  */
  double Bk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Reduced frequency-dependent Breit two-body matrix element

    @details
    Calculates the reduced two-body Breit matrix element evaluated at the frequency
    corresponding to the energy difference between electrons a and c.
    The frequency w is scaled by lambda_f.

    @param k  Multipolarity (angular momentum rank)
    @param A  First electron state
    @param Fb Second electron state
    @param Fc Third electron state
    @param Fd Fourth electron state

    @return B^k_abcd evaluated at w = lambda_f * alpha * |E_a - E_c|

  */
  [[deprecated]] double Bk_abcd_eac_freqw(int k, const DiracSpinor &Fa,
                                          const DiracSpinor &Fb,
                                          const DiracSpinor &Fc,
                                          const DiracSpinor &Fd) const;

  /*!
    @brief Reduced frequency-dependent Breit two-body matrix element

    @details
    Calculates the reduced two-body Breit matrix element evaluated at the frequency
    corresponding to the energy difference between electrons b and d.
    The frequency w is scaled by lambda_f.

    @param k  Multipolarity (angular momentum rank)
    @param Fa First electron state
    @param Fb Second electron state
    @param Fc Third electron state
    @param Fd Fourth electron state

    @return B^k_abcd evaluated at w = lambda_f * alpha * |E_b - E_d|
  */
  [[deprecated]] double Bk_abcd_ebd_freqw(int k, const DiracSpinor &Fa,
                                          const DiracSpinor &Fb,
                                          const DiracSpinor &Fc,
                                          const DiracSpinor &Fd) const;

  /*!
    @brief Reduced exchange Breit two-body matrix element

    @details
    Calculates the reduced two-body exchange Breit matrix element P(B)^k_{abcd},
    the static Breit analogue of the Coulomb P^k integral.
  */
  double BPk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                  const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Reduced anti-symmetrised Breit two-body matrix element

    @details
    Calculates the reduced anti-symmetrised Breit matrix element W(B)^k_{abcd}
    for fermionic systems, the static Breit analogue of the Coulomb W^k integral.
  */
  double BWk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                  const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Direct Breit one-body potential operator

    @details
    Computes the radial function of the direct part of the reduced Breit
    operator acting on an orbital with quantum number kappa_v. This is defined
    such that the two-body matrix element factorises as:

    \f[
      B^k_{abcd} = \langle a | B^k_v(b,c,d) \rangle
    \f]

    See @ref Bk_abcd()
  */
  DiracSpinor Bkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                      const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Frequency-dependent Breit one-body potential operator

    @details
    Frequency-dependent breit analogue of @ref Bkv_bcd()
  */
  DiracSpinor Bkv_bcd_freqw(int k, int kappa_v, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd,
                            const double w) const;

  // Should add two-body matrix elements for frequency-dependent Breit at some point.
  // As well as writing function for direct contribution to HF

  /*!
    @brief Exchange Breit one-body potential operator

    @details
    See @ref Bkv_bcd() and @ref BPk_abcd()
  */
  DiracSpinor BPkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                       const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Frequency-dependent exchange Breit one-body potential operator

    @details
    Frequency-dependent version of BPkv_bcd(). See @ref Bkv_bcd_freqw().
  */
  DiracSpinor BPkv_bcd_freqw(int k, int kappa_v, const DiracSpinor &Fb,
                             const DiracSpinor &Fc, const DiracSpinor &Fd,
                             const double w) const;

  /*!
    @brief Anti-symmetrised Breit one-body potential operator

    @details
    See @ref Bkv_bcd() and @ref BWk_abcd()
  */
  DiracSpinor BWkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                       const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Frequency-dependent anti-symmetrised Breit one-body potential operator

    @details
    Frequency-dependent version of BWkv_bcd(). See @ref Bkv_bcd_freqw().
  */
  DiracSpinor BWkv_bcd_freqw(int k, int kappa_v, const DiracSpinor &Fb,
                             const DiracSpinor &Fc, const DiracSpinor &Fd,
                             const double w) const;

  /*!
    @brief Reduced static Breit matrix element (tabulated, fast lookup)

    @details
    A faster implementation of Bk_abcd() that looks up pre-tabulated
    integrals computed by fill_gb(). Reduces computational cost, at significant
    memory cost.

    The result is the same as Bk_abcd() for the same inputs, provided
    fill_gb() has been called at least once on a basis containing Fa, Fb, Fc, Fd.

    @param k  Multipolarity (angular momentum rank of the interaction)
    @param Fa (Fb,Fc,Fd) electron states

    @return The reduced Breit matrix element B^k_{abcd}, retrieved from
            the cached m_gb tables.

    @note REQUIRES: fill_gb() must have been called before this function
          can be used. If fill_gb() has not been called, this function returns
          zero or garbage data.
          This function assumes the spinors Fa, Fb, Fc, Fd were members of
          the basis passed to fill_gb(); using spinors outside that basis
          is undefined.

    @warning Do not use this function unless fill_gb() has been successfully called.
  */
  double Bk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                   const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Reduced exchange Breit matrix element (tabulated, fast lookup)

    @details
    A much faster implementation of BPk_abcd() that looks up pre-tabulated
    exchange integrals computed by fill_gb(). Reduces computational cost from
    O(radial grid size) per call to O(1) lookup.

    The result is identical to BPk_abcd() for the same inputs, provided
    fill_gb() has been called at least once on a basis containing all four
    electron states.

    @param k  Multipolarity (angular momentum rank of the interaction)
    @param Fa (Fb,Fc,Fd) electron states

    @return The reduced exchange Breit matrix element P(B)^k_{abcd}, retrieved
            from the cached m_gb tables.

    @note REQUIRES: fill_gb() must have been called before this function
          can be used. The spinors must be members of the basis provided to fill_gb().

    @warning Do not use this function unless fill_gb() has been successfully called.
  */
  double BPk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief Reduced anti-symmetrised Breit matrix element (tabulated, fast lookup)

    @details
    A much faster implementation of BWk_abcd() that looks up pre-tabulated
    integrals computed by fill_gb(). Combines the tabulated direct and exchange
    terms for O(1) evaluation.

    The result is identical to BWk_abcd() for the same inputs, provided
    fill_gb() has been called at least once on a basis containing all four
    electron states.

    @param k  Multipolarity (angular momentum rank of the interaction)
    @param Fa (Fb,Fc,Fd) electron states

    @return The reduced anti-symmetrised Breit matrix element
            W(B)^k_{abcd} = Bk_abcd_2(...) + BPk_abcd_2(...), retrieved
            from cached m_gb tables.

    @note REQUIRES: fill_gb() must have been called before this function
          can be used. The spinors must be members of the basis provided to fill_gb().

    @warning Do not use this function unless fill_gb() has been successfully called.
  */
  double BWk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) const;

  /*!
    @brief The one-body Breit (Breit-Hartree-Fock) correction to second-order energy

    @details
    Calculates the one-body Breit (Breit-Hartree-Fock) correction to second-order 
    energy for a single-valence atom. This is included automatically if Breit is included into Hartree-Fock.
    (The lowest-order Breit correction <v|V_Br|v> not included.)
    

    @param Fv      The valence electron state of interest
    @param holes   Occupied core electron states
    @param excited Virtual excited electron states

    @return The one-body Breit correction to second-order energy: delta E^(2, 1).

    @note This represents only the one-particle (Hartree-Fock) contribution.
          The two-particle (Sigma) contribution is computed by de2().
          For Breit-Coulomb self-consistent HF, this is already accounted for
          and should not be added separately.

    @warning May significantly overcount or undercount energy shifts if
             Breit-Coulomb HF is not self-consistent. Use primarily as
             a correction term in perturbative Breit-on-Coulomb calculations.
  */
  double de2_HF(const DiracSpinor &Fv, const std::vector<DiracSpinor> &holes,
                const std::vector<DiracSpinor> &excited) const;

  /*!
    @brief The two-body Breit correction to second-order energy

    @details
    Calculates the two-body Breit correction to second-order 
    energy for a single-valence atom.
    This is _not_ included automatically if Breit is included into Hartree-Fock, 
    since it requires modification of two-body Coulomb integals.
    
    @note This is just the energy shift - Breit can be included fully into MBPT 
    calculations self-consistantly another way (see @ref MBPT)
  */
  double de2(const DiracSpinor &Fv, const std::vector<DiracSpinor> &holes,
             const std::vector<DiracSpinor> &excited) const;
};

} // namespace HF
