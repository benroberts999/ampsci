#pragma once
#include "CI/CI_Integrals.hpp"
#include "CI/CSF.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Potentials/RadPot.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "json/json.hpp"
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
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
  Also holds the Hartree-Fock potential (and core), basis, spectrum,
  correlation potential (Sigma), and CI solutions.

  \par Construction:
    - Set of GridParameters [see Maths/Grid]
    - Nuclear::Nucleus [see Potentials/NuclearPotentials]
    - var_alpha = \f$\lambda\f$, \f$\alpha = \lambda\alpha_0\f$
    - run_label:  Optional label for output identity - for distinguishing
      outputs with different parameters
*/
class Wavefunction {

public:
  /*!
    @brief Construct with a Grid [shared resource] and a nucleus.

    @details
    Nucleus: isotope data etc. Optional fractional variation in alpha:
    alpha = var_alpha * alpha_0, alpha_0 =~ 1/137. Forms the nuclear
    potential; orbitals are formed later (see set_HF, solve_core).

    @param grid       Radial grid (shared resource)
    @param nucleus    Nuclear parameters (isotope, charge distribution)
    @param var_alpha  Ratio alpha/alpha_0
    @param run_label  Optional label for output identity (see identity)
  */
  Wavefunction(std::shared_ptr<const Grid> grid,
               const Nuclear::Nucleus &nucleus, double var_alpha = 1.0,
               const std::string &run_label = "");
  //! As above, but Grid is constructed here using given parameters
  Wavefunction(const GridParameters &gridparams,
               const Nuclear::Nucleus &nucleus, double var_alpha = 1.0,
               const std::string &run_label = "");

  //! Default grid and nucleus (see GridParameters, Nuclear::Nucleus)
  Wavefunction() : Wavefunction(GridParameters{}, Nuclear::Nucleus{}) {}

private:
  // Radial grid
  std::shared_ptr<const Grid> rgrid;
  // Internal value for alpha (alpha = var_alpha * alpha_0, alpha_0=~1/137)
  double m_alpha;
  std::string m_run_label;
  // Holds nuclear parameters (isotope, charge distro etc.)
  Nuclear::Nucleus m_nucleus;
  // Valence (single-particle) orbitals
  std::vector<DiracSpinor> m_valence{};
  std::vector<DiracSpinor> m_hf_valence{};
  // Basis, daigonalised over HF core. Used for MBPT
  std::vector<DiracSpinor> m_basis{};
  // Sprectrum: like basis, but includes Sigma (correlations).
  std::vector<DiracSpinor> m_spectrum{};
  // Nuclear potential // here AND hf?
  std::vector<double> m_vnuc{};
  // Hartree-Fock potential
  std::optional<HF::HartreeFock> m_HF{std::nullopt};
  // Correlation potential; for now unique_ptr; prefer std::optional
  std::optional<MBPT::CorrelationPotential> m_Sigma{};
  // Core configuration (non-rel terms)
  std::string m_core_string = "";
  std::string m_aboveFermi_core_string = "";

  // CI solutions, and the integrals used to construct the CI Hamiltonians
  CI::Solutions m_CI{};

public:
  //! Returns a const reference to the radial grid
  const Grid &grid() const { return *rgrid; };
  /*!
    @brief Copy of shared_ptr to grid [shared resource].

    @details
    Used when we want to construct a new object that shares this grid.
  */
  std::shared_ptr<const Grid> grid_sptr() const { return rgrid; };

  //! Local value of fine-structure constant.
  double alpha() const { return m_alpha; }

  //! Variation in alpha^2 : x = (alpha/alpha_0)^2 - 1
  double dalpha2() const {
    // (alpha/alpha_0)^2 -1
    return (m_alpha * m_alpha / PhysConst::alpha2) - 1.0;
  }

  //! Returns Nuclear::nucleus object (contains nuc. parameters)
  const Nuclear::Nucleus &nucleus() const { return m_nucleus; }
  //! Nuclear charge, Z
  int Znuc() const { return m_nucleus.z(); }
  //! Nuclear mass number, A
  int Anuc() const { return m_nucleus.a(); }
  //! Nuclear rms charge radius, in fm (femtometres)
  double get_rrms() const { return m_nucleus.r_rms(); }

  //! Core orbitals (frozen HF core); empty if HF not set
  const std::vector<DiracSpinor> &core() const {
    static const auto empty = std::vector<DiracSpinor>{}; //?
    return m_HF ? m_HF->core() : empty;
  }

  //! Valence orbitals (HF or Brueckner orbitals)
  const std::vector<DiracSpinor> &valence() const { return m_valence; }
  std::vector<DiracSpinor> &valence() { return m_valence; }

  //! HF valence orbitals: copy of valence, as last solved by solve_valence
  const std::vector<DiracSpinor> &hf_valence() const { return m_hf_valence; }

  /*!
    @brief Basis, eigenstates of HF potential. Used for MBPT.

    @details
    Includes Breit and QED (if they are included), but not correlations.
  */
  const std::vector<DiracSpinor> &basis() const { return m_basis; }
  std::vector<DiracSpinor> &basis() { return m_basis; }

  //! Spectrum: like basis, but includes correlations
  const std::vector<DiracSpinor> &spectrum() const { return m_spectrum; }
  std::vector<DiracSpinor> &spectrum() { return m_spectrum; }

  //! CI solutions, one for each {J, parity}; empty unless CI was run
  const std::vector<CI::PsiJPi> &CIwfs() const { return m_CI.levels; }

  /*!
    @brief CI solutions for given J and parity; nullptr if not found.

    @details
    @p J is integer (CI is for two-electron systems); @p parity is +1 or -1.
  */
  const CI::PsiJPi *CIwf(int J, int parity) const {
    for (const auto &ci_wf : m_CI.levels) {
      if (ci_wf.twoJ() == 2 * J && ci_wf.parity() == parity)
        return &ci_wf;
    }
    return nullptr;
  }

  /*!
    @brief Integral tables used to construct the CI Hamiltonians.

    @details
    Empty unless CI was run; see CI::Integrals.
  */
  const CI::Integrals &CI_integrals() const { return m_CI.integrals; }

  /*!
    @brief Nuclear potential.

    @details
    Only provide const version, since HF and WF version of vnuc must be kept
    in sync (see update_Vnuc).
  */
  const std::vector<double> &vnuc() const { return m_vnuc; }

  //! Returns ptr to Hartree Fock (class); nullptr if not set
  const HF::HartreeFock *vHF() const { return m_HF ? &*m_HF : nullptr; }
  HF::HartreeFock *vHF() { return m_HF ? &*m_HF : nullptr; }

  /*!
    @brief Local part of potential, e.g., Vl = Vnuc + Vdir + Vrad_el(l).

    @details
    Can be l-dependent. Returns a copy. Just Vnuc if HF is not set.
  */
  std::vector<double> vlocal(int l = 0) const;

  /*!
    @brief QED Magnetic form factor.

    @details
    May return empty vector (no QED, or HF not set). Not typically
    l-dependent, but may be in future. Returns a copy.
  */
  std::vector<double> Hmag(int l = 0) const;

  //! Pointer to QED radiative potential. May be nullptr
  const QED::RadPot *vrad() const { return m_HF ? m_HF->Vrad() : nullptr; }
  QED::RadPot *vrad() { return m_HF ? m_HF->Vrad() : nullptr; }

  //! Returns ptr to (const) Correlation Potential, Sigma; nullptr if none
  const MBPT::CorrelationPotential *Sigma() const {
    return m_Sigma ? &*m_Sigma : nullptr;
  }
  MBPT::CorrelationPotential *Sigma() { return m_Sigma ? &*m_Sigma : nullptr; }

  //----------------------------------

  //! Number of electrons in the core
  int Ncore() const;

  /*!
    @brief Finds requested state (n, kappa); returns nullptr if not found.

    @details
    Searches valence first, then core.
  */
  const DiracSpinor *getState(int n, int k) const;
  //! As above, but takes 'short symbol' (e.g., 6s+, 6p-)
  const DiracSpinor *getState(std::string_view state) const;

  /*!
    @brief Returns energy location of the "Fermi Level".

    @details
    Energy half way between core/valence. Defined:
    0.5*( max(e_core) + min(e_valence)), plus a small offset (1e-5; or 1e-3
    if there are no valence states, in which case min(e_valence) is taken as
    zero). If there is no core: 1.1*min(e_valence). Should be -ve.
  */
  double FermiLevel() const;

  /*!
    @brief Energy gap between lowest valence + highest core state.

    @details
    e(v) - e(c) [should be positive]. Zero if there are no core or no valence
    states.
  */
  double energy_gap() const;

  //! Returns full core configuration
  std::string coreConfiguration() const { return m_core_string; }

  //! Returns core configuration, in nice output notation
  std::string coreConfiguration_nice() const {
    return AtomData::niceCoreOutput(m_core_string);
  }

  //! String of atom info (e.g., "Cs, Z=55 A=133")
  std::string atom() const {
    return AtomData::atomicSymbol(m_nucleus.z()) +
           ", Z=" + std::to_string(m_nucleus.z()) +
           " A=" + std::to_string(m_nucleus.a());
  }

  //! e.g., "Cs"
  std::string atomicSymbol() const {
    return AtomData::atomicSymbol(m_nucleus.z());
  }

  /*!
    @brief Atomic symbol, including core ionisation degree and run_label.

    @details
    Used as the base name for output and cache files (Sigma, Qk, CI, etc.).
    Form: symbol + Zion + Breit/QED flag ("b", "q", or "bq") + alpha
    variation ("_x", with x = (alpha/alpha_0)^2 to 4 significant figures, if
    dalpha2() > 1e-5) + "_run_label" (if not empty). E.g., "Cs1",
    "Cs1bq_label".

    @warning Alpha variation is only labelled for alpha > alpha_0: runs with
    alpha < alpha_0 get the same identity (and so read the same cache files)
    as alpha = alpha_0.
  */
  std::string identity() const;

  //! Optional label for output identity (see identity)
  const std::string &run_label() const { return m_run_label; };

  //! 0 for neutral, 1 for singly-ionised etc. (num_val: valence electrons)
  int ion_degree(int num_val) const { return Zion() - num_val; }

  //! I for neutral, II for singly-ionised etc. (num_val: valence electrons)
  std::string ion_symbol(int num_val) const {
    return qip::int_to_roman(Zion() - num_val + 1);
  }

  //! Effective charge (for core) = Z-N_core
  int Zion() const { return Znuc() - Ncore(); }

  //! Prints table of core orbitals + energies etc.
  void printCore() const;
  /*!
    @brief Prints table of valence orbitals + energies etc.

    @details
    Can optionally give it any list of orbitals to print (default: valence).
  */
  void printValence(const std::vector<DiracSpinor> &tmp_orbitals = {}) const;

  //! Prints table of Basis/Spectrum orbitals, compares to HF orbitals
  void printBasis(const std::vector<DiracSpinor> &the_basis) const;

  //! Check if a state (n, kappa) is in the core list
  bool isInCore(int n, int k) const;
  //! Check if a state (n, kappa) is in the valence list
  bool isInValence(int n, int k) const;

  //! Calculates rho(r) = sum_c psi^2(r) for core states, c={n,k,m}
  std::vector<double> coreDensity() const;

  /*!
    @brief Calculates HF core energy.

    @details
    Sum of core orbital energies (weighted by occupation), minus the
    double-counted Coulomb (direct + exchange) interaction. QED (if included
    in the core) enters through the orbital energies. Zero if HF is not set.

    @note The double-counting correction is for Coulomb only: Breit (if
          included) enters only through the orbital energies, without the
          factor 1/2.
  */
  double coreEnergyHF() const;

  //------------------------------------------------------------------

  /*!
    @brief Initialises HF object and populates core orbitals.

    @details
    Does not solve HF equations (see solve_core). Replaces any existing HF
    object, including any radiative potential (QED is added later; see
    radiativePotential).

    @param method        Mean-field method (HartreeFock, Hartree, KohnSham..)
    @param in_core       Core configuration, "[X],nLm,nLm,...", e.g., "[Xe]"
                         or "[Xe],6s1". A ':' may replace one of the commas:
                         shells after it are in the core, but may also be
                         valence states (see solve_valence)
    @param breit_params  Breit parameters (nullopt: no Breit). HartreeFock
                         method only: a warning is printed otherwise
    @param eps_HF        HF convergence goal (if >= 1, taken as the exponent:
                         10^-eps_HF)
    @param print         Print method details

    @warning Aborts if the core configuration is invalid, or has more than Z
             electrons.
  */
  void set_HF(HF::Method method = HF::Method::HartreeFock,
              const std::string &in_core = "",
              std::optional<HF::Breit::Params> breit_params = std::nullopt,
              double eps_HF = 1.0e-13, bool print = true);

  //! Performs Hartree-Fock procedure for core (nothing if HF not set)
  void solve_core(bool print = true);

  //! This version will first set_HF(), then solve_core()
  void solve_core(HF::Method method, const std::string &in_core = "",
                  std::optional<HF::Breit::Params> breit_params = std::nullopt,
                  double eps_HF = 1.0e-13, bool print = true);
  //! @overload
  void solve_core(const std::string &method, const std::string &in_core = "",
                  std::optional<HF::Breit::Params> breit_params = std::nullopt,
                  double eps_HF = 1.0e-13, bool print = true);

  /*!
    @brief Performs Hartree-Fock procedure for valence; populates valence.

    @details
    Adds the states in @p in_valence_str (basis-string format, e.g., "7sp5d":
    up to n=7 for s and p, n=5 for d) that are not already in valence, and
    not in the core (except those after the ':' in the core string; see
    set_HF). Then solves for all valence states (existing ones are re-solved)
    in the HF potential, without Sigma, and copies them to hf_valence. If HF
    is not set, solves in the nuclear potential (H-like).
  */
  void solve_valence(const std::string &in_valence_str = "",
                     const bool print = true);

  /*!
    @brief Solves for exotic atoms (e.g., muonic), including screening.

    @details
    Resulting states are included in valence (flagged as exotic); the
    screening also updates core.

    Each state is first solved in the nuclear (+QED) potential. With
    @p screening (and a non-empty core), the core HF equations are iterated
    including the direct potential of the first exotic state; the HF object
    (core) is then replaced by the screened one (Vnuc itself unchanged), and
    the valence states are re-solved in the screened potential.

    Note: The exotic states are just added to the valence list, so they can
    be used more simply with all the modules. However, be careful; for
    example, RPA will now be meaningless for these states!

    @param in_exotic_str  States, in basis-string format (e.g., "1s2p")
    @param mass           Mass of the exotic particle, in units of m_e
    @param screening      Include screening (core and exotic states)
    @param print          Print energies and screening iterations

    @warning With @p screening, all valence states are re-solved with the
             exotic mass, including any (electron) valence states already
             present: valence should be empty before calling this.
  */
  void solve_exotic(const std::string &in_exotic_str,
                    double mass = PhysConst::m_muon, bool screening = true,
                    bool print = true);

  /*!
    @brief Forms Brueckner valence orbitals; replaces existing valence states.

    @details
    \f$ (H_{\rm HF} + \Sigma)\ket{nk} = \varepsilon\ket{nk} \f$.
    Does nothing if Sigma has not been formed; prints a warning and does
    nothing if HF is not set.

    @note Sigma is only included for the HartreeFock method (otherwise the
          valence states are re-solved in the local potential).
  */
  void hartreeFockBrueckner(const bool print = true);

  /*!
    @brief First, fits Sigma to energies, then forms fitted Brueckner orbitals.

    @details
    For each of the first fit_energies.size() valence states (in order),
    finds the scaling factor lambda such that the Brueckner energy, with
    lambda*Sigma, matches the given energy (damped iteration; lambda kept
    within [0.25, 2.5]). By default, fitted lambdas are averaged over each
    fine-structure pair (same n and l); separate_lambda = true keeps them
    separate. Then calls hartreeFockBrueckner.

    @param valence_list     Unused
    @param fit_energies     Energies to fit to (au), in the same order as
                            valence; sign ignored (binding energies)
    @param separate_lambda  Keep lambdas of fine-structure pairs separate

    @note Each fitted state needs its own Sigma (same kappa and n; see
          each_valence in formSigma), otherwise its fit has no effect.
    @warning Sigma must have been formed (not checked).
  */
  void fitSigma_hfBrueckner(const std::string &valence_list,
                            const std::vector<double> &fit_energies,
                            bool separate_lambda = false);

  //! OLD: deprecated
  void radiativePotential(QED::RadPot::Scale s, double rcut, double scale_rN,
                          const std::vector<double> &x_spd,
                          bool do_readwrite = true, bool print = true);

  /*!
    @brief Calculates radiative potential, adds to HF potential.

    @details
    Replaces any existing radiative potential in HF. Included in subsequent
    HF solutions: call before solve_core to include QED in the core, after
    it for valence only. Prints a warning and does nothing if HF is not set.

    @param qed_input     QED options (see QED::ConstructRadPot)
    @param do_readwrite  Read/write the potential from/to file
    @param print         Print details
  */
  void radiativePotential(const IO::InputBlock &qed_input, bool do_readwrite,
                          bool print);

  //! Calculates + populates basis [see BSplineBasis]
  void formBasis(const SplineBasis::Parameters &params);

  //! Calculates + populates Spectrum [see BSplineBasis]
  void formSpectrum(const SplineBasis::Parameters &params);

  /*!
    @brief Forms + stores correlation potential Sigma.

    @details
    Constructs Sigma (see MBPT::CorrelationPotential), reading it from file if
    it exists, and forms it for the valence states (see @p each_valence,
    @p ek). Then applies @p lambdas, optionally prints the energy shifts, and
    writes Sigma to file. Requires the basis (formBasis) for the Goldstone
    sums. Does nothing if the core is empty or HF is not set.

    @param nmin_core           Minimum n of core states included in
                               polarisation
    @param r0                  Sub-grid first point (aB)
    @param rmax                Sub-grid last point (aB)
    @param stride              Sub-grid stride: every stride-th point of the
                               full grid
    @param each_valence        Form Sigma for each valence state; else for the
                               first valence state of each kappa. Ignored if
                               @p ek is given
    @param include_G           Include small-component (g) parts of Sigma
    @param include_Breit       Include two-body Breit correction, B2 (only if
                               HF includes Breit)
    @param n_max_breit         Maximum n of excited states in B2 (0: all)
    @param lambdas             Scaling factors, lambda, for each Sigma, in the
                               order formed
    @param fk                  Effective screening factors for exchange, for
                               each k (empty: calculated for Feynman with
                               screening)
    @param etak                Effective hole-particle factors, for each k
                               (Goldstone only)
    @param read_write          Read Sigma from file (if it exists) and write
                               it; false: no file
    @param fname               Sigma file name (empty: identity() plus an
                               extension encoding the method)
    @param FeynmanQ            Feynman method for direct diagrams; else
                               Goldstone
    @param ScreeningQ          Include all-order screening (Feynman only)
    @param hole_particleQ      Include hole-particle interaction (Feynman
                               only)
    @param lmax                Maximum l of internal lines (Feynman)
    @param omre                Real part of frequency for the Im(w)
                               integration (Feynman)
    @param w0                  First non-zero point of the Im(w) grid
                               (Feynman)
    @param wratio              Ratio of the logarithmic Im(w) grid (Feynman)
    @param complex_green       Solve Green's function directly at complex
                               energy (see MBPT::FeynmanOptions)
    @param ek                  State/energy pairs (energies in au): form Sigma
                               for each at the given energy, instead of the
                               valence energies
    @param ladder_file         File with ladder correction, Sigma_L (empty:
                               none)
    @param derivative          Also form dSigma/dE
    @param print_de            Print MBPT(2) energy corrections for valence
    @param fk_both_lines       Apply fk to both Coulomb lines of the exchange
                               diagrams (else outer only)
    @param tail_extrapolation  Extrapolate Sigma beyond the sub-grid as
                               a_eff/r^4 (experimental)

    @note If read from file, the sub-grid and @p include_G are taken from the
          file (see MBPT::CorrelationPotential).
  */
  void formSigma(int nmin_core = 1, double r0 = 1.0e-4, double rmax = 30.0,
                 int stride = 4, bool each_valence = false,
                 bool include_G = false, bool include_Breit = false,
                 int n_max_breit = 0, const std::vector<double> &lambdas = {},
                 const std::vector<double> &fk = {},
                 const std::vector<double> &etak = {}, bool read_write = true,
                 const std::string &fname = "", bool FeynmanQ = false,
                 bool ScreeningQ = false, bool hole_particleQ = false,
                 int lmax = 6, double omre = -0.2, double w0 = 0.01,
                 double wratio = 1.5, bool complex_green = false,
                 const std::optional<IO::InputBlock> &ek = std::nullopt,
                 const std::string &ladder_file = "", bool derivative = false,
                 bool print_de = false, bool fk_both_lines = false,
                 bool tail_extrapolation = false);

  // void correlations(const IO::InputBlock &input);

  /*!
    @brief Copies correlation potential Sigma (nothing if @p Sigma is null).

    @details
    E.g., to re-use Sigma from another Wavefunction (different isotope).

    @note The copy still refers to the Hartree-Fock object it was formed with,
          which must outlive it if it is used to form or write Sigma.
  */
  void copySigma(const MBPT::CorrelationPotential *const Sigma) {
    if (Sigma != nullptr)
      m_Sigma = *Sigma;
  }

  /*!
    @brief Replaces nuclear potential, in both Wavefunction and HartreeFock.

    @details
    Allows extra potential to be added to Vnuc: pass vnuc() + extra. Updates
    both in Wavefunction _and_ HartreeFock, which must be kept in sync.
    Orbitals are not re-solved, and nucleus() is not updated.
  */
  void update_Vnuc(const std::vector<double> &v_new) {
    /// nb: two versions of Vnuc...
    m_vnuc = v_new;
    if (m_HF) {
      m_HF->vnuc() = v_new;
    }
  }

  /*!
    @brief Returns [min,max] r values for which the core density is > cutoff.

    @details
    Core density for given l; cutoff = eps*max_value.
    Returns the r values (au) for which the value of
    \f$ \rho(r) = \sum|\psi|^2 \f$ drops below cutoff. Sum goes over all m for
    given l. Cut-off defined as eps*max, where max is maximum value for
    rho(r). Set l<0 to get for all l (entire core). Returns the full grid
    range if there are no core states with given l (or if eps = 0).
  */
  std::tuple<double, double> lminmax_core_range(int l, double eps = 0.0) const;

  /*!
    @brief Returns <a|H|b> for Hamiltonian H (NOT Sigma, Breit, or exchange).

    @details
    Includes the local potential, vlocal (incl. Rad.pot), and the QED
    magnetic form factor, Hmag. Uses the symmetric form of the derivative
    term (integrated by parts). Zero if @p Fa and @p Fb have different kappa.
  */
  double H0ab(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  /*!
    @brief As above, but with the radial derivatives given.

    @details
    @p dFa and @p dFb are the derivatives (d/dr) of @p Fa and @p Fb; only
    their g components are used.
  */
  double H0ab(const DiracSpinor &Fa, const DiracSpinor &dFa,
              const DiracSpinor &Fb, const DiracSpinor &dFb) const;

  /*!
    @brief Returns <a|H|b>, including exchange, Sigma, and Breit.

    @details
    H0ab, plus HF exchange (if HF is non-local), Sigma (if formed, including
    lambda), and Breit (if included in HF).

    @warning Only meaningful for @p Fa and @p Fb of the same kappa: H0ab is
             then zero, but the other terms are still added (the DiracSpinor
             product does not check kappa).
  */
  double Hab(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  //! Runs the CI+MBPT routines; stores wavefunctions
  void ConfigurationInteraction(const IO::InputBlock &input);

  /*!
    @brief Writes wavefunction information to json file.

    @details
    Grid, nuclear parameters, and core, valence, and exotic (muon) orbitals,
    with a description of each field ("metadata").

    @param out_name  If given (not empty), will print to that file
    @return The json object.
  */
  nlohmann::json
  output_to_json(const std::string &out_name = "ampsci_output.json");

private:
  double H0ab_impl(const DiracSpinor &Fa, std::vector<double> dga,
                   const DiracSpinor &Fb, std::vector<double> dgb) const;

  // Creates set of blank core orbitals
  std::vector<DiracSpinor> determineCore(const std::string &str_core_in);
  bool isInAboveFermiCore(int n, int k) const;
};
