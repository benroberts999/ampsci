#pragma once
#include "HF/HartreeFock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/RadPot.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
namespace SplineBasis {
struct Parameters;
}

//==============================================================================
/*!
@brief Stores Wavefunction (set of valence orbitals, grid, HF etc.)
@details
\par Construction:
  - Set of GridParameters [see Maths/Grid]
  - Set of Nuclear::Nucleus [see Physics/NuclearPotentials]
  - var_alpha = \f$\lambda\f$, \f$\alpha = \lambda\alpha_0\f$

*/
class Wavefunction {

public:
  //! Construct with a Grid [shared resource], a nucleus (isotope data etc.),
  //! and (optional) fractional variation in alpha  [alpha = var_alpha *
  //! alpha_0, alpha_0=~1/137]
  Wavefunction(std::shared_ptr<const Grid> grid,
               const Nuclear::Nucleus &nucleus, double var_alpha = 1.0);
  //! As above, but Grid is constructed here using given parameters
  Wavefunction(const GridParameters &gridparams,
               const Nuclear::Nucleus &nucleus, double var_alpha = 1.0);

  //! User-defined copy-constructor. Note: Does not copy Sigma
  Wavefunction(const Wavefunction &wf); // XXX make sigma copyable!?
  //! Deleted, as sigma cannot be copied. This will be fixed
  Wavefunction &operator=(const Wavefunction &) = delete;
  ~Wavefunction() = default;

private:
  // Radial grid
  std::shared_ptr<const Grid> rgrid;
  // Internal value for alpha (alpha = var_alpha * alpha_0, alpha_0=~1/137)
  double m_alpha;
  // Holds nuclear parameters (isotope, charge distro etc.)
  Nuclear::Nucleus m_nucleus;
  // Valence (single-particle) orbitals
  std::vector<DiracSpinor> m_valence{};
  // Basis, daigonalised over HF core. Used for MBPT
  std::vector<DiracSpinor> m_basis{};
  // Sprectrum: like basis, but includes Sigma (correlations).
  std::vector<DiracSpinor> m_spectrum{};
  // Nuclear potential // here AND hf?
  std::vector<double> m_vnuc{};
  // Hartree-Fock potential
  std::optional<HF::HartreeFock> m_HF{std::nullopt};
  // Correlation potential; for now unique_ptr; prefer std::optional
  std::unique_ptr<MBPT::CorrelationPotential> m_Sigma{nullptr};
  // Core configuration (non-rel terms)
  std::string m_core_string = "";

public:
  //! Returns a const reference to the radial grid
  const Grid &grid() const { return *rgrid; };
  //! Copy of shared_ptr to grid [shared resource] - used when we want to
  //! construct a new object that shares this grid
  std::shared_ptr<const Grid> grid_sptr() const { return rgrid; };

  //! Local value of fine-structure constant.
  double alpha() const { return m_alpha; }

  const Nuclear::Nucleus &nucleus() const { return m_nucleus; }
  int Znuc() const { return m_nucleus.z(); }
  int Anuc() const { return m_nucleus.a(); }
  double get_rrms() const { return m_nucleus.r_rms(); }

  //! Core orbitals (frozen HF core)
  const std::vector<DiracSpinor> &core() const {
    static const auto empty = std::vector<DiracSpinor>{}; //?
    return m_HF ? m_HF->core() : empty;
  }

  //! Valence orbitals (HF or Brueckner orbitals)
  const std::vector<DiracSpinor> &valence() const { return m_valence; }
  std::vector<DiracSpinor> &valence() { return m_valence; }

  //! Basis, eigenstates of HF potential. Used for MBPT. Includes Breit and
  //! QED (if they are included), but not correlations
  const std::vector<DiracSpinor> &basis() const { return m_basis; }
  std::vector<DiracSpinor> &basis() { return m_basis; }

  //! Sprectrum: like basis, but includes correlations
  const std::vector<DiracSpinor> &spectrum() const { return m_spectrum; }
  std::vector<DiracSpinor> &spectrum() { return m_spectrum; }

  //! Nuclear potential. Only provide const version, since HF and WF version of
  //! vnuc must be kept in sync
  const std::vector<double> &vnuc() const { return m_vnuc; }

  //! Returns ptr to Hartree Fock (class)
  const HF::HartreeFock *vHF() const { return m_HF ? &*m_HF : nullptr; }
  HF::HartreeFock *vHF() { return m_HF ? &*m_HF : nullptr; }

  //! Local part of potential, e.g., Vl = Vnuc + Vdir + Vrad_el(l) - can be
  //! l-dependent. Returns a copy
  std::vector<double> vlocal(int l = 0) const;

  //! QED Magnetic form factor. May return empty vector. Not typically
  //! l-dependent, but may be in future. Returns a copy
  std::vector<double> Hmag(int l = 0) const;

  //! Pointer to QED radiative potnential. May be nullptr
  const QED::RadPot *vrad() const { return m_HF ? m_HF->Vrad() : nullptr; }
  QED::RadPot *vrad() { return m_HF ? m_HF->Vrad() : nullptr; }

  //! Returns ptr to (const) Correlation Potential, Sigma
  const MBPT::CorrelationPotential *Sigma() const { return m_Sigma.get(); }
  MBPT::CorrelationPotential *Sigma() { return m_Sigma.get(); }

  //----------------------------------

  //! Number of electrons in the core
  int Ncore() const;

  //! Finds requested state; returns nullptr if not found
  //! @details is_valence is optional out-parameter; tells you where orb was
  //! found
  const DiracSpinor *getState(int n, int k) const;
  //! As above, but takes 'short symbol' (e.g., 6s+, 6p-)
  const DiracSpinor *getState(std::string_view state) const;

  //! Returns energy location of the "Fermi Level", - energy half way between core/valence. Defined: 0.5*( max(e_core) + min(e_valence)). Should be -ve
  double FermiLevel() const;

  //! Energy gap between lowest valence + highest core state: e(v) - e(c)
  //! [should be positive]
  double energy_gap() const;

  //! Returns full core configuration
  std::string coreConfiguration() const { return m_core_string; }

  //! Returns core configuration, in nice output notation
  std::string coreConfiguration_nice() const {
    return AtomData::niceCoreOutput(m_core_string);
  }

  //! String of atom info (e.g., "Cs, Z=55, A=133")
  std::string atom() const {
    return AtomData::atomicSymbol(m_nucleus.z()) +
           ", Z=" + std::to_string(m_nucleus.z()) +
           " A=" + std::to_string(m_nucleus.a());
  }
  //! e.g., "Cs"
  std::string atomicSymbol() const {
    return AtomData::atomicSymbol(m_nucleus.z());
  }

  //! E.g., Cs in V^N-1, gives Cs-i
  std::string identity() const {
    const auto zionRoman = qip::int_to_roman(Zion());
    return AtomData::atomicSymbol(m_nucleus.z()) + zionRoman;
  }

  //! Effective charge (for core) = Z-N_core
  int Zion() const { return Znuc() - Ncore(); }

  //! Prints table of core orbitals + energies etc.
  void printCore() const;
  //! @brief Prints table of valence orbitals + energies etc.
  //! @details Can optionally give it any list of orbitals to print
  void printValence(const std::vector<DiracSpinor> &tmp_orbitals = {}) const;

  //! Prints table of Basis/Spectrum orbitals, compares to HF orbitals
  void printBasis(const std::vector<DiracSpinor> &the_basis) const;

  //! Check if a state is in the core (or valence) list
  bool isInCore(int n, int k) const;
  bool isInValence(int n, int k) const;

  //! Calculates rho(r) = sum_c psi^2(r) for core states, c={n,k,m}
  std::vector<double> coreDensity() const;

  //! Calculates HF core energy (doesn't include magnetic QED?)
  double coreEnergyHF() const;

  //------------------------------------------------------------------

  //! Initialises HF object and populates core orbitals (does not solve HF
  //! equations)
  void set_HF(const std::string &method = "HartreeFock",
              const double x_Breit = 0.0, const std::string &in_core = "",
              double eps_HF = 1.0e-13, bool print = true);

  //! Performs hartree-Fock procedure for core
  void solve_core(bool print = true);

  //! This version will first set_HF(), then solve_core()
  void solve_core(const std::string &method, const double x_Breit = 0.0,
                  const std::string &in_core = "", double eps_HF = 1.0e-13,
                  bool print = true);

  //! Performs hartree-Fock procedure for valence: note: poplulates valnece
  void solve_valence(const std::string &in_valence_str = "",
                     const bool print = true);

  //! Forms Bruckner valence orbitals: (H_hf + Sigma)|nk> = e|nk>. Replaces
  //! existing valence states
  void hartreeFockBrueckner(const bool print = true);

  //! First, fits Sigma to energies, then forms fitted Brueckner orbitals
  void fitSigma_hfBrueckner(const std::string &valence_list,
                            const std::vector<double> &fit_energies);

  //! Second-order MBPT energy shifts, calculates + prints
  void SOEnergyShift();

  //! OLD: deprecated
  void radiativePotential(QED::RadPot::Scale s, double rcut, double scale_rN,
                          const std::vector<double> &x_spd,
                          bool do_readwrite = true, bool print = true);

  //! Calculates radiative potential, adds to HF potential
  void radiativePotential(const IO::InputBlock &qed_input, bool do_readwrite,
                          bool print);

  //! Calculates + populates basis [see BSplineBasis]
  void formBasis(const SplineBasis::Parameters &params);

  //! Calculates + populates Spectrum [see BSplineBasis]
  void formSpectrum(const SplineBasis::Parameters &params);

  //! Forms + stores correlation potential Sigma
  void formSigma(
      const int nmin_core = 1, const bool form_matrix = true,
      const double r0 = 1.0e-4, const double rmax = 30.0, const int stride = 4,
      const bool each_valence = false, const bool include_G = false,
      const std::vector<double> &lambdas = {},
      const std::vector<double> &fk = {}, const std::vector<double> &etak = {},
      const std::string &in_fname = "", const std::string &out_fname = "",
      const std::string &ladder_file = "", const bool FeynmanQ = false,
      const bool ScreeningQ = false, const bool holeParticleQ = false,
      const int lmax = 6, const bool GreenBasis = false,
      const bool PolBasis = false, const double omre = -0.2, double w0 = 0.01,
      double wratio = 1.5,
      const std::optional<IO::InputBlock> &ek = std::nullopt);

  void copySigma(const MBPT::CorrelationPotential *const Sigma) {
    if (Sigma != nullptr)
      m_Sigma = std::make_unique<MBPT::CorrelationPotential>(*Sigma);
  }

  //! Allows extra potential to be added to Vnuc (updates both in Wavefunction
  // _and_ HartreeFock)
  void update_Vnuc(const std::vector<double> &v_new) {
    /// nb: two versions of Vnuc...
    m_vnuc = v_new;
    if (m_HF) {
      m_HF->vnuc() = v_new;
    }
  }

  //! @brief Returns [min,max] r values for which the core density (given l) is
  //! larger than cutoff (= eps*max_value)
  //! @details Returns the r values (au) for which the value of rho =
  //! \sum|psi^2|(r) drops below cutoff. Sum goes over all m for given l.
  //! Cut-off defined as eps*max, where max is maximum value for rho(r).
  //! Set l<0 to get for all l (entire core)
  std::tuple<double, double> lminmax_core_range(int l, double eps = 0.0) const;

  //! Returns <a|H|b> for Hamiltonian H (inludes Rad.pot, NOT sigma or Breit)
  double Hab(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  double Hab(const DiracSpinor &Fa, const DiracSpinor &dFa,
             const DiracSpinor &Fb, const DiracSpinor &dFb) const;

private:
  // Creates set of blank core orbitals
  std::vector<DiracSpinor> determineCore(const std::string &str_core_in);
};
