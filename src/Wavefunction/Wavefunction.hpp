#pragma once
#include "HF/HartreeFock.hpp" // forward decl..
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp" // NonRelSEConfig
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp" //PhysConst::alpha
#include "Physics/RadPot.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
namespace SplineBasis {
struct Parameters;
}

//******************************************************************************
/*!
@brief Stores Wavefunction (set of core+valence orbitals, grid etc.)
@details
\par Construction:
  - Set of GridParameters [see Maths/Grid]
  - Set of Nuclear::Parameters [see Physics/NuclearPotentials]
  - var_alpha = \f$\lambda\f$, \f$\alpha = \lambda\alpha_0\f$

*/
class Wavefunction {

public:
  Wavefunction(const GridParameters &gridparams,
               const Nuclear::Parameters &nuc_params, double var_alpha = 1.0);

  //! User-defined copy-constructor. Note: Does not copy HF or Sigma
  Wavefunction(const Wavefunction &wf);
  Wavefunction &operator=(const Wavefunction &) = delete;
  ~Wavefunction() = default;

public:
  //! "Frozen" Core orbitals
  std::vector<DiracSpinor> core{};
  //! Valence (single-particle) orbitals
  std::vector<DiracSpinor> valence{};
  //! Basis, daigonalised over HF core. Used for MBPT
  std::vector<DiracSpinor> basis{};
  //! Sprectrum: like basis, but includes Sigma (correlations).
  std::vector<DiracSpinor> spectrum{};
  //! Radial grid
  std::shared_ptr<const Grid> rgrid;
  //! Internal value for alpha (alpha = var_alpha * alpha_0, alpha_0=~1/137)
  const double alpha;

private:
  const Nuclear::Parameters m_nuclear;
  std::unique_ptr<HF::HartreeFock> m_pHF{nullptr};
  std::unique_ptr<MBPT::CorrelationPotential> m_Sigma{nullptr};

public:
  //! Nuclear potential
  std::vector<double> vnuc{};
  //! Direct/local part of the electron potential
  std::vector<double> vdir{};
  //! QED radiative potential
  std::unique_ptr<QED::RadPot> qed{nullptr};

private:
  // Core configuration (non-rel terms)
  std::vector<AtomData::NonRelSEConfig> m_core_configs = {};
  int num_core_electrons = 0; // Nc = N - M
  std::string m_core_string = "";

public: // const methods: "views" into WF object
  // Rule is: if function is single-line, define here. Else, in .cpp
  int Znuc() const { return m_nuclear.z; }
  int Anuc() const { return m_nuclear.a; }

  //! Number of electrons in the core
  int Ncore() const { return num_core_electrons; }
  const Nuclear::Parameters &get_nuclearParameters() const { return m_nuclear; }
  double get_rrms() const { return m_nuclear.r_rms; }

  bool exclude_exchangeQ() const {
    if (m_pHF == nullptr)
      return true;
    return m_pHF->excludeExchangeQ();
  }

  //! Returns ptr to (const) Correlation Potential, Sigma
  const MBPT::CorrelationPotential *getSigma() const { return m_Sigma.get(); }
  //! Returns ptr to (const) Hartree Fock (class)
  const HF::HartreeFock *getHF() const { return m_pHF.get(); }

  //! Finds requested state; returns nullptr if not found
  //! @details is_valence is optional out-parameter; tells you where orb was
  //! found
  const DiracSpinor *getState(int n, int k, bool *is_valence = nullptr) const;
  //! As above, but takes 'short symbol' (e.g., 6s+, 6p-)
  const DiracSpinor *getState(std::string_view state,
                              bool *is_valence = nullptr) const;

  //! Returns energy location of the "core-valence gap", 0.5*( max(e_core) +
  //! min(e_valence)) - energy half way between core/valence
  double en_coreval_gap() const;

  //! Energy gap between lowest valence + highest core state
  double energy_gap() const {
    const auto c =
        std::max_element(cbegin(core), cend(core), DiracSpinor::comp_en);
    const auto v =
        std::min_element(cbegin(valence), cend(valence), DiracSpinor::comp_en);
    if (c != cend(core) && v != cend(valence))
      return v->en() - c->en();
    return 0.0;
  }

  //! Returns full core configuration
  std::string coreConfiguration() const { return m_core_string; }
  //! Returns core configuration, in nice output notation
  std::string coreConfiguration_nice() const {
    return AtomData::niceCoreOutput(m_core_string);
  }
  //! Outputs screen-friendly nuclear parameters
  std::string nuclearParams() const;

  //! String of atom info (e.g., "Cs, Z=55, A=133")
  std::string atom() const {
    return AtomData::atomicSymbol(m_nuclear.z) +
           ", Z=" + std::to_string(m_nuclear.z) +
           " A=" + std::to_string(m_nuclear.a);
  }
  //! e.g., "Cs"
  std::string atomicSymbol() const {
    return AtomData::atomicSymbol(m_nuclear.z);
  }
  //! Effective charge (for core) = Z-N_core
  int Zion() const { return Znuc() - Ncore(); }
  //! E.g., Cs in V^N-1, gives Cs-i
  std::string identity() const {
    const auto zionRoman = AtomData::int_to_roman(Zion());
    return AtomData::atomicSymbol(m_nuclear.z) + zionRoman;
  }

  //! Prints table of core orbitals + energies etc. Optionally sorted by energy
  void printCore(bool sorted = true) const;
  //! @brief Prints table of valence orbitals + energies etc. Optionally sorted
  //! by energy
  //! @details Can optionally give it any list of orbitals to print
  void printValence(bool sorted = true,
                    const std::vector<DiracSpinor> &tmp_orbitals = {}) const;
  //! Prints table of Basis/Spectrum orbitals, compares to HF orbitals
  void printBasis(const std::vector<DiracSpinor> &the_basis,
                  bool sorted = false) const;

  bool isInCore(int n, int k) const;
  bool isInValence(int n, int k) const;

  //! Largest n for core states (optional: for given kappa, otherwise overall)
  int maxCore_n(int ka_in = 0) const;
  //! Largest l for core states
  int maxCore_l() const;

  //! Calculates rho(r) = sum_c psi^2(r) for core states, c={n,k,m}
  std::vector<double> coreDensity() const;

  //! Performs hartree-Fock procedure for core: note: poplulates core
  void solve_core(const std::string &method = "HartreeFock",
                  const double x_Breit = 0.0, const std::string &in_core = "",
                  double eps_HF = 0, bool print = true);

  //! Calculates HF core energy (doesn't include magnetic QED?)
  auto coreEnergyHF() const;

  //! Performs hartree-Fock procedure for valence: note: poplulates valnece
  void solve_valence(const std::string &in_valence_str = "",
                     const bool print = true);
  //! Solves new local valence (e.g., Kohn-Sham): note: poplulates valence
  void localValence(const std::string &in_valence_str, bool list_each = false);
  //! Forms Bruckner valence orbitals: (H_hf + Sigma)|nk> = e|nk>.
  void hartreeFockBrueckner(const bool print = true);
  //! First, fits Sigma to energies, then forms fitted Brueckner orbitals
  void fitSigma_hfBrueckner(const std::string &valence_list,
                            const std::vector<double> &fit_energies);
  //! Second-order MBPT energy shifts, calculates + prints
  void SOEnergyShift();

  //! Calculates radiative potential. Stores in vnuc, and Hmag
  void radiativePotential(QED::RadPot::Scale s, double rcut, double scale_rN,
                          const std::vector<double> &x_spd,
                          bool do_readwrite = true, bool print = true);

  //! Calculates + populates basis [see BSplineBasis]
  void formBasis(const SplineBasis::Parameters &params);

  //! Calculates + populates Spectrum [see BSplineBasis]
  void formSpectrum(const SplineBasis::Parameters &params);

  //! Forms + stores correlation potential Sigma
  void formSigma(const int nmin_core = 1, const bool form_matrix = true,
                 const double r0 = 1.0e-4, const double rmax = 30.0,
                 const int stride = 4, const bool each_valence = false,
                 const bool include_G = false,
                 const std::vector<double> &lambdas = {},
                 const std::vector<double> &fk = {},
                 const std::string &in_fname = "",
                 const std::string &out_fname = "", const bool FeynmanQ = false,
                 const bool ScreeningQ = false,
                 const bool holeParticleQ = false, const int lmax = 6,
                 const bool GreenBasis = false, const bool PolBasis = false,
                 const double omre = -0.2, double w0 = 0.01,
                 double wratio = 1.5,
                 const std::optional<IO::InputBlock> &ek = std::nullopt);
  void copySigma(const MBPT::CorrelationPotential *const Sigma) {
    if (Sigma != nullptr)
      m_Sigma = std::make_unique<MBPT::CorrelationPotential>(*Sigma);
  }

  //! @brief Solves Dirac bound state problem, with optional 'extra' potential
  //! log_eps is log_10(convergence_target).
  void solveDirac(DiracSpinor &psi, double e_a, const std::vector<double> &vex,
                  int log_eps = 0) const;
  void solveDirac(DiracSpinor &psi, double e_a = 0, int log_eps = 0) const;

  //! Populates core orbitals accorind to given core string (+solves)
  void solveLocalCore(const std::string &str_core_in, int log_dele_or = 0);
  //! Adds new valence orbtial (+solves using vdir)
  void solveNewValence(int n, int k, double en_a = 0, int log_dele_or = 0);

  //! energy guess for a core state with n,l; quite rough, but good enough
  double enGuessCore(int n, int l) const;
  //! energy guess for a valence state with n,l
  double enGuessVal(int n, int ka) const;

  void add_to_Vdir(const std::vector<double> &dv) { //
    /// XXX Fix: two versions of Vdir...
    qip::add(&vdir, dv);
    if (m_pHF) {
      m_pHF->update_Vdir(vdir);
    }
  }

  //! (approximately) OrthoNormalises a set of any orbitals.
  //! @details Note: only updates orbs, not energies
  static void orthonormaliseOrbitals(std::vector<DiracSpinor> &in_orbs,
                                     int num_its = 1);
  //! (exactly) OrthoNormalises psi_v against of any orbitals.
  static void orthonormaliseWrt(DiracSpinor &psi_v,
                                const std::vector<DiracSpinor> &in_orbs);
  //! (exactly) OrthoGonalises psi_v against of any orbitals (no Norm).
  static void orthogonaliseWrt(DiracSpinor &psi_v,
                               const std::vector<DiracSpinor> &in_orbs);

  //! @brief Returns [min,max] r values for which the core density (given l) is
  //! larger than cutoff (= eps*max_value)
  //! @details Returns the r values (au) for which the value of rho =
  //! \sum|psi^2|(r) drops below cutoff. Sum goes over all m for given l.
  //! Cut-off defined as eps*max, where max is maximum value for rho(r).
  //! Set l<0 to get for all l (entire core)
  std::tuple<double, double> lminmax_core_range(int l, double eps = 0.0) const;

  //! Local potential, e.g., Vl = Vnuc + Vdir + Vrad_el(l) - can be l-dependent
  std::vector<double> get_Vlocal(int l = 0) const;
  std::vector<double> get_Hmag(int l = 0) const;

  //! Returns <a|H|b> for Hamiltonian H (inludes Rad.pot, NOT sigma or Breit)
  double Hab(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  double Hab(const DiracSpinor &Fa, const DiracSpinor &dFa,
             const DiracSpinor &Fb, const DiracSpinor &dFb) const;

private:
  void determineCore(const std::string &str_core_in);
  static std::vector<std::size_t>
  sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                   bool do_sort = false);
};
