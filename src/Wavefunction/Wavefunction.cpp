#include "Wavefunction/Wavefunction.hpp"
#include "CI/CI.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //just for enum..
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Potentials/Parametric_potentials.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
Wavefunction::Wavefunction(const GridParameters &gridparams,
                           const Nuclear::Nucleus &t_nucleus, double var_alpha,
                           const std::string &run_label)
    : Wavefunction(std::make_shared<const Grid>(gridparams), t_nucleus,
                   var_alpha, run_label) {}

Wavefunction::Wavefunction(std::shared_ptr<const Grid> in_grid,
                           const Nuclear::Nucleus &t_nucleus, double var_alpha,
                           const std::string &run_label)
    : rgrid(std::move(in_grid)),
      m_alpha(PhysConst::alpha * var_alpha),
      m_run_label(run_label),
      m_nucleus(t_nucleus),
      m_vnuc(Nuclear::formPotential(m_nucleus, rgrid->r())) {}

//==============================================================================
std::vector<DiracSpinor>
Wavefunction::determineCore(const std::string &str_core_in)
// Takes in a string list for the core configuration, outputs an int list
// Takes in previous closed shell (noble), + 'rest' (or just the rest)
// E.g:
//   Core of Cs: Xe
//   Core of Gold: Xe 4f14 5d10
// 'rest' is in form nLm : n=n, L=l, m=number of electrons in that nl shell.
{
  std::vector<DiracSpinor> core;

  // std::string a, b;
  std::string str_core = str_core_in;

  const auto fermi_pos = str_core.find(':');
  if (fermi_pos != std::string::npos) {
    m_aboveFermi_core_string = str_core.substr(fermi_pos + 1);
    str_core.replace(fermi_pos, 1, ",");
  }

  const auto core_configs = AtomData::core_parser(str_core);

  bool bad_core = false;
  m_core_string = "";
  for (auto config : core_configs) {
    if (!config.ok()) {
      m_core_string += " **";
      bad_core = true;
    }
    m_core_string += config.symbol();
    if (config != core_configs.back())
      m_core_string += ",";
  }

  if (bad_core) {
    std::cout << "Problem with core: " << str_core_in << " = \n";
    std::cout << m_core_string << "\n";
    std::cout << "In this house, we obey the Pauli exclusion principle!\n";
    std::abort();
  }

  // Count number of electrons in the core
  int num_core_electrons = 0;
  for (const auto &config : core_configs) {
    num_core_electrons += config.num;
  }

  if (num_core_electrons > m_nucleus.z()) {
    std::cout << "Problem with core: " << str_core_in << "\n";
    std::cout << "= " << m_core_string << "\n";
    std::cout << "= " << AtomData::niceCoreOutput(m_core_string) << "\n";
    std::cout << "Too many electrons: N_core=" << num_core_electrons
              << ", Z=" << m_nucleus.z() << "\n";
    std::abort();
  }

  // for (const auto &[n, l, num] : core_configs) {
  // This case of Structured Bindings fails with GCC<=7.4 ?
  for (auto &config : core_configs) {
    const auto n = config.n;
    const auto l = config.l;
    const auto num = config.num;
    if (num == 0)
      continue;
    int k1 = l; // j = l-1/2
    if (k1 != 0) {
      auto &new_Fc = core.emplace_back(n, k1, rgrid);
      new_Fc.occ_frac() = double(num) / (4 * l + 2);
    }
    int k2 = -(l + 1); // j=l+1/2
    auto &new_Fc = core.emplace_back(n, k2, rgrid);
    new_Fc.occ_frac() = double(num) / (4 * l + 2);
  }

  return core;
}

//==============================================================================
void Wavefunction::set_HF(const std::string &method, const double x_Breit,
                          const std::string &in_core, double eps_HF,
                          bool print) {

  auto core = determineCore(in_core);
  const auto qed = std::nullopt; // we add QED (optionally) later - to allow for
                                 // QED into valence, but not core
  m_HF = HF::HartreeFock(rgrid, m_vnuc, std::move(core), qed, m_alpha,
                         HF::parseMethod(method), x_Breit, eps_HF);

  // Move this into HF?
  if (print) {
    // Print some HF into to screen:
    if (method == "Hartree")
      std::cout << "Using Hartree Method (no Exchange)\n";
    else if (method == "ApproxHF")
      std::cout << "Using approximate HF Method (approx Exchange)\n";
    else if (method == "KohnSham") {
      std::cout
          << "Using Kohn-Sham Method.\n"
          << "Note: You should include first valence state into the core:\n"
             "Kohn-Sham is NOT a V^N-1 method!\n";
    } else if (method == "Local") {
      std::cout << "Using local potential\n";
    } else if (method != "HartreeFock") {
      fmt2::styled_print(fg(fmt::color::orange), "\nWARNING\n");
      std::cout << "unkown method: " << method
                << "\nDefaulting to HartreeFock method.\n";
    }

    // Can only include Breit within HF
    if (method == "HartreeFock" && x_Breit != 0.0) {
      std::cout << "Including Breit (scale = " << x_Breit << ")\n";
    } else if (method != "HartreeFock" && x_Breit != 0.0) {
      fmt2::styled_print(fg(fmt::color::orange), "\nWARNING\n");
      std::cout << "can only include Breit in Hartree-Fock "
                   "method. Breit will not be included.\n";
    }
  }
}

//==============================================================================
void Wavefunction::solve_core(bool print) {
  if (m_HF)
    m_HF->solve_core(print);
}

void Wavefunction::solve_core(const std::string &method, const double x_Breit,
                              const std::string &in_core, double eps_HF,
                              bool print) {
  set_HF(method, x_Breit, in_core, eps_HF, print);
  solve_core(print);
}

//==============================================================================
double Wavefunction::coreEnergyHF() const {
  if (!m_HF) {
    return 0.0;
  }
  return m_HF->calculateCoreEnergy();
}

//==============================================================================
void Wavefunction::solve_valence(const std::string &in_valence_str,
                                 const bool print) {

  const auto explicite_val_list = false;
  // (m_HF && m_HF->method() == HF::Method::KohnSham);
  const auto val_lst = explicite_val_list ?
                           AtomData::listOfStates_singlen(in_valence_str) :
                           AtomData::listOfStates_nk(in_valence_str);

  // 1. populate valence states
  for (const auto &[n, k, en] : val_lst) {
    (void)en;
    // if (!isInValence(n, k) && (!isInCore(n, k) || explicite_val_list)) {
    if (!isInValence(n, k) && (!isInCore(n, k) || isInAboveFermiCore(n, k))) {
      // For Kohn-Sham, valence state may be in core
      auto &Fv = m_valence.emplace_back(n, k, rgrid);
      Fv.occ_frac() = 1.0 / Fv.twojp1();
    }
  }

  // 2. Solve HF
  if (m_HF) {
    m_HF->solve_valence(&m_valence, print);
  } else {
    // only for H-like:
    const double z2 = m_nucleus.z() * m_nucleus.z();
    for (auto &Fv : m_valence) {
      const auto e0 = -0.5 * z2 / std::pow(Fv.n(), 2);
      DiracODE::boundState(Fv, e0, vlocal(), {}, m_alpha, 1.0e-14);
    }
  }
  m_hf_valence = m_valence;
}

//==============================================================================
void Wavefunction::radiativePotential(QED::RadPot::Scale scale, double rcut,
                                      double scale_rN,
                                      const std::vector<double> &x_spd,
                                      bool do_readwrite, bool print) {
  if (x_spd.empty())
    return;

  const auto r_N_au =
      std::sqrt(5.0 / 3.0) * scale_rN * m_nucleus.r_rms() / PhysConst::aB_fm;

  auto qed = QED::RadPot(rgrid->r(), Znuc(), r_N_au, rcut, scale, x_spd, print,
                         do_readwrite);

  // If HF already exists, update it to include new qed!
  if (m_HF) {
    m_HF->set_Vrad(std::move(qed));
  } else {
    std::cout << "\nWarning: Can only include QED with a HF method\n";
  }
}

//==============================================================================
void Wavefunction::radiativePotential(const IO::InputBlock &qed_input,
                                      bool do_readwrite, bool print) {

  //re-scale happends inside ConstructRadPot
  const auto r_N_au =
      std::sqrt(5.0 / 3.0) * m_nucleus.r_rms() / PhysConst::aB_fm;

  auto qed = QED::ConstructRadPot(rgrid->r(), Znuc(), r_N_au, qed_input, print,
                                  do_readwrite);

  // If HF already exists, update it to include new qed!
  if (m_HF) {
    m_HF->set_Vrad(std::move(qed));
  } else {
    std::cout << "\nWarning: Can only include QED with a HF method\n";
  }
}

//==============================================================================
bool Wavefunction::isInCore(int n, int k) const {
  const auto find_nk = [n, k](const auto &Fa) {
    return Fa.n() == n && Fa.kappa() == k;
  };
  const auto Fnk = std::find_if(cbegin(core()), cend(core()), find_nk);
  return Fnk != cend(core());
}
//------------------------------------------------------------------------------
bool Wavefunction::isInValence(int n, int k) const {
  const auto find_nk = [n, k](const auto Fa) {
    return Fa.n() == n && Fa.kappa() == k;
  };
  const auto Fnk = std::find_if(cbegin(m_valence), cend(m_valence), find_nk);
  return Fnk != cend(m_valence);
}

//------------------------------------------------------------------------------
bool Wavefunction::isInAboveFermiCore(int n, int k) const {
  std::istringstream ss(m_aboveFermi_core_string);
  std::string each;
  auto sn = std::to_string(n) + AtomData::l_symbol(Angular::l_k(k));
  while (std::getline(ss, each, ',')) {
    if (each.substr(0, sn.length()) == sn)
      return true;
  }
  return false;
}

//==============================================================================
const DiracSpinor *Wavefunction::getState(int n, int k) const {
  const auto find_nk = [n, k](const auto Fa) {
    return Fa.n() == n && Fa.kappa() == k;
  };
  // Try to find in valence first:
  // nb: search valence first
  auto Fnk = std::find_if(cbegin(m_valence), cend(m_valence), find_nk);
  if (Fnk != cend(m_valence)) {
    return &*Fnk;
  }
  // If not in valence, try to find in core:
  Fnk = std::find_if(cbegin(core()), cend(core()), find_nk);
  if (Fnk != cend(core())) {
    return &*Fnk;
  }

  // otherwise, return nope
  return nullptr;
}

const DiracSpinor *Wavefunction::getState(std::string_view state) const {
  // std::pair<int, int> parse_symbol(std::string_view symbol);
  const auto [n, k] = AtomData::parse_symbol(state);
  return getState(n, k);
}

//==============================================================================
double Wavefunction::FermiLevel() const {
  // Find core/valence energy: allows distingush core/valence states
  const auto ec_max =
      core().empty() ?
          0.0 :
          std::max_element(cbegin(core()), cend(core()), DiracSpinor::comp_en)
              ->en();
  const auto ev_min = m_valence.empty() ?
                          0.0 :
                          std::min_element(cbegin(m_valence), cend(m_valence),
                                           DiracSpinor::comp_en)
                              ->en();
  // return 0.5 * (ev_min + ec_max);
  if (core().empty()) {
    return 1.1 * ev_min;
  }
  const auto extra = m_valence.empty() ? 1.0e-3 : 1.0e-5;
  return ec_max + 0.5 * (ev_min - ec_max) + extra;
}
double Wavefunction::energy_gap() const {

  const auto c =
      std::max_element(cbegin(core()), cend(core()), DiracSpinor::comp_en);

  const auto v = std::min_element(cbegin(m_valence), cend(m_valence),
                                  DiracSpinor::comp_en);
  if (c != cend(core()) && v != cend(m_valence))
    return v->en() - c->en();
  return 0.0;
}

//==============================================================================
std::tuple<double, double> Wavefunction::lminmax_core_range(int l,
                                                            double eps) const {
  std::vector<double> rho_l(rgrid->num_points());
  bool found = false;
  for (const auto &Fc : core()) {
    if (Fc.l() == l || l < 0) {
      found = true;
      qip::add(&rho_l, Fc.rho());
    }
  }
  if (!found) {
    return {rgrid->r().front(), rgrid->r().back()};
  }
  // find maximum rho:
  const auto max = std::max_element(rho_l.begin(), rho_l.end());
  // find first position that rho=|psi^2| reaches cut-off
  const auto cut = max != rho_l.end() ? eps * (*max) : 0.0;
  const auto lam = [cut](const auto &v) { return v >= cut; };
  const auto first = std::find_if(rho_l.begin(), rho_l.end(), lam);
  const auto last = std::find_if(rho_l.rbegin(), rho_l.rend(), lam);
  const auto index_first = std::size_t(first - rho_l.begin());
  const auto index_last = std::size_t(rho_l.rend() - last);
  return {rgrid->r()[index_first], rgrid->r()[index_last]};
}

//==============================================================================
void Wavefunction::printCore() const
// prints core orbitals
{
  // int Zion = Znuc() - Ncore();
  std::cout << "Core: " << coreConfiguration_nice() << " V^N";
  if (Zion() != 0) {
    std::cout << "-" << Zion();
  }

  if (Ncore() == 0) {
    std::cout << " - H-like";
  }
  std::cout << "\n";
  if (Ncore() < 1)
    return;

  std::cout
      << "     state  k   Rinf its   eps         En (au)        En (/cm)\n";
  int i = 0;
  for (const auto &phi : core()) {
    printf("%-2i %7s %2i  %5.1f %2i  %5.0e %15.9f %15.3f", i++,
           phi.symbol().c_str(), phi.kappa(), phi.rinf(), phi.its(), phi.eps(),
           phi.en(), phi.en() * PhysConst::Hartree_invcm);
    if (phi.occ_frac() < 1.0) {
      printf("     [%4.2f]\n", phi.occ_frac());
    } else {
      std::cout << "\n";
    }
  }
}

//==============================================================================
void Wavefunction::printValence(
    const std::vector<DiracSpinor> &in_orbitals) const {
  const auto &tmp_orbs = (in_orbitals.empty()) ? m_valence : in_orbitals;
  if (tmp_orbs.empty())
    return;

  std::cout << "Valence: " << atomicSymbol() << ion_symbol(1) << "\n";

  // Find lowest valence energy:
  const auto min_it =
      std::min_element(tmp_orbs.begin(), tmp_orbs.end(), DiracSpinor::comp_en);
  const auto e0 = (min_it == tmp_orbs.end()) ? 0.0 : min_it->en();

  std::cout
      << "     state  "
      << "k   Rinf its   eps         En (au)        En (/cm)   En (/cm)\n";
  int i = 0;
  for (const auto &phi : tmp_orbs) {
    printf("%-2i %7s %2i  %5.1f %2i  %5.0e %15.9f %15.3f", i++,
           phi.symbol().c_str(), phi.kappa(), phi.rinf(), phi.its(), phi.eps(),
           phi.en(), phi.en() * PhysConst::Hartree_invcm);
    printf(" %10.2f\n", (phi.en() - e0) * PhysConst::Hartree_invcm);
  }
}

//==============================================================================
void Wavefunction::printBasis(const std::vector<DiracSpinor> &the_basis) const {
  std::cout
      << "     State   k  R0     Rinf          En(Basis)         En(HF)\n";
  int i = 0;
  for (const auto &phi : the_basis) {

    const auto *hf_phi = getState(phi.n(), phi.kappa());
    if (hf_phi != nullptr) {
      // found HF state
      const auto eps =
          2.0 * (phi.en() - hf_phi->en()) / (phi.en() + hf_phi->en());
      printf("%2i) %7s %2i  %5.0e %5.1f %18.7f  %13.7f  %6.0e\n", i,
             phi.symbol().c_str(), phi.kappa(), phi.r0(), phi.rinf(), phi.en(),
             hf_phi->en(), eps);
    } else {
      printf("%2i) %7s %2i  %5.0e %5.1f %18.7f\n", i, phi.symbol().c_str(),
             phi.kappa(), phi.r0(), phi.rinf(), phi.en());
    }
    ++i;
  }
}

//==============================================================================
std::vector<double> Wavefunction::coreDensity() const {
  std::vector<double> rho(rgrid->num_points(), 0.0);
  for (const auto &phi : core()) {
    auto f = double(phi.twoj() + 1) * phi.occ_frac();
    for (auto i = 0ul; i < rgrid->num_points(); i++) {
      rho[i] += f * (phi.f(i) * phi.f(i) + phi.g(i) * phi.g(i));
    }
  }
  return rho;
}

//==============================================================================
void Wavefunction::formBasis(const SplineBasis::Parameters &params) {
  if (params.n > 0) {
    IO::ChronoTimer t("Basis");
    m_basis = SplineBasis::form_basis(params, *this, false);

    if (params.orthogonalise) {
      std::cout << "Forcing spectrum to be orthog to core:\n";
      for (auto &Fb : m_basis) {
        Fb = DiracSpinor::orthonormaliseWrt(Fb, core());
      }
    }

    std::cout << "Basis/core:\n";
    SplineBasis::check(m_basis, core(), true);
    std::cout << "Basis/valence:\n";
    SplineBasis::check(m_basis, valence(), true);
  }
}
//------------------------------------------------------------------------------
void Wavefunction::formSpectrum(const SplineBasis::Parameters &params) {
  if (params.n > 0) {
    IO::ChronoTimer t("Spectrum");
    m_spectrum = SplineBasis::form_basis(params, *this, true);
  }

  if (params.orthogonalise) {
    std::cout << "Forcing spectrum to be orthog to valence:\n";
    for (auto &Fb : m_spectrum) {
      Fb = DiracSpinor::orthonormaliseWrt(Fb, m_valence);
    }
  }

  std::cout << "Spectrum/core:\n";
  SplineBasis::check(m_spectrum, core(), m_Sigma == std::nullopt);
  std::cout << "Spectrum/valence:\n";
  SplineBasis::check(m_spectrum, valence(), true);
}

//==============================================================================
void Wavefunction::formSigma(
    int nmin_core, int nmin_core_F, double r0, double rmax, int stride,
    bool each_valence, bool include_G, bool include_Breit, int n_max_breit,
    const std::vector<double> &lambdas, const std::vector<double> &fk,
    const std::vector<double> &etak, const std::string &in_fname,
    const std::string &out_fname, bool FeynmanQ, bool ScreeningQ,
    bool holeParticleQ, int lmax, double omre, double w0, double wratio,
    const std::optional<IO::InputBlock> &ek) {
  if (m_valence.empty())
    return;
  std::cout << "\nIncluding correlation potential:\n" << std::flush;

  std::string ext = FeynmanQ ? ".sf" : ".s2";
  if (include_G)
    ext += "g";
  if (FeynmanQ && ScreeningQ)
    ext += "s";
  if (FeynmanQ && holeParticleQ)
    ext += "h";
  if (include_Breit && m_HF->vBreit() && !FeynmanQ)
    ext += "b";

  ext += ".abf";

  const auto ifname = in_fname == "" ? identity() + ext : in_fname + ext;
  const auto ofname = out_fname == "" ? identity() + ext : out_fname + ext;

  const auto method =
      FeynmanQ ? MBPT::SigmaMethod::Feynman : MBPT::SigmaMethod::Goldstone;

  MBPT::Screening screening =
      ScreeningQ ? MBPT::Screening::include : MBPT::Screening::exclude;

  // XXX Update to allow all k
  // MBPT::HoleParticle hp = holeParticleQ ? MBPT::HoleParticle::include_k0 :
  //                                         MBPT::HoleParticle::exclude;

  MBPT::HoleParticle hp =
      holeParticleQ ? MBPT::HoleParticle::include : MBPT::HoleParticle::exclude;

  bool calculate_fk = FeynmanQ && fk.empty();

  m_Sigma = MBPT::CorrelationPotential(
      ifname, &*m_HF, m_basis, r0, rmax, std::size_t(stride), nmin_core, method,
      include_G, include_Breit, n_max_breit,
      MBPT::FeynmanOptions{screening, hp, lmax, omre, w0, wratio}, calculate_fk,
      fk, etak, nmin_core_F);

  std::cout << "\n";

  if (ek)
    each_valence = false;

  // This is for each valence state.... otherwise, just do for lowest??
  if (!m_valence.empty()) {
    if (each_valence) {
      // calculate sigma for each valence state:
      for (const auto &Fv : m_valence) {
        m_Sigma->formSigma(Fv.kappa(), Fv.en(), Fv.n(), &Fv);
      }
    } else if (!ek && m_Sigma->empty()) {
      // calculate sigma for lowest n valence state of each kappa:
      const auto max_ki = DiracSpinor::max_kindex(m_valence);
      for (int ki = 0; ki <= max_ki; ++ki) {
        auto Fv = std::find_if(cbegin(m_valence), cend(m_valence),
                               [ki](auto f) { return f.k_index() == ki; });
        if (Fv != cend(m_valence))
          m_Sigma->formSigma(Fv->kappa(), Fv->en(), Fv->n(), &*Fv);
      }
    } else if (ek && m_Sigma->empty()) {
      // solve at specific energies:
      for (auto &[state, en] : ek->options()) {
        auto [n, k] = AtomData::parse_symbol(state);

        auto Fv = std::find_if(cbegin(m_valence), cend(m_valence),
                               [ki = k](auto f) { return f.k_index() == ki; });
        if (Fv != cend(m_valence)) {
          m_Sigma->formSigma(k, std::stod(en), n, &*Fv);
        } else {
          m_Sigma->formSigma(k, std::stod(en), n);
        }
      }
    }
  }

  if (!lambdas.empty()) {
    std::cout << "Note: be careful order of input lambdas matches Sigma\n";
    // If Sigma diverged from valence, may be diff order to valence...
    // Better to make it explicit!
    m_Sigma->scale_Sigma(lambdas);
    m_Sigma->print_scaling();
  }

  if (out_fname != "false")
    m_Sigma->write(ofname);
}

//==============================================================================
void Wavefunction::hartreeFockBrueckner(const bool print) {
  if (!m_HF) {
    std::cerr << "WARNING 62: Cant call solve_valence before solve_core\n";
    return;
  }
  if (m_Sigma)
    m_HF->solve_valence(&m_valence, print, &*m_Sigma);
}

//==============================================================================
void Wavefunction::fitSigma_hfBrueckner(
    const std::string &, const std::vector<double> &fit_energies) {
  std::cout << "Fitting Sigma for lowest valence states:\n" << std::flush;

  const auto max_its = 30;
  const auto eps_targ = 1.0e-7;

  // XXX Assume the 'fit_to' are in same order as valence!!
  // Must be called after HF, before Bruckenr....

  //
  for (auto i = 0ul; i < fit_energies.size(); ++i) {
    if (i >= m_valence.size())
      break;
    const auto &Fv = m_valence[i];
    const auto e_exp = -std::abs(fit_energies[i]);

    const double en_0 = Fv.en(); // HF value
    auto e_Sig1 = 0.0;
    auto lambda = 1.0;
    const double a_damp = 0.5; // 1 means no damping
    double eps = 1.0;
    int its = 0;
    for (; its <= max_its; its++) {
      auto Fv_l = Fv;
      m_Sigma->scale_Sigma(lambda, Fv_l.kappa(), Fv_l.n());
      m_HF->hf_valence(Fv_l, &*m_Sigma);
      double en_l = Fv_l.en();
      if (its == 0)
        e_Sig1 = en_l;
      eps = std::abs((e_exp - en_l) / e_exp);
      if (eps < eps_targ || its == max_its)
        break;
      const auto r = (e_exp - en_l) / (en_l - en_0);
      const auto new_lambda = lambda * (1.0 + a_damp * r);
      lambda = std::clamp(new_lambda, 0.5, 1.5);
    }
    printf("%4s | %8.1f, %8.1f [%8.1f] : ", Fv.shortSymbol().c_str(),
           en_0 * PhysConst::Hartree_invcm, e_Sig1 * PhysConst::Hartree_invcm,
           e_exp * PhysConst::Hartree_invcm);
    printf("%.0e (%2i); lambda = %.4f", eps, its, lambda);
    if (eps > 1.0e-6)
      std::cout << " **";
    if (eps > 1.0e-5)
      std::cout << "***";
    std::cout << "\n";
  }
  std::cout << "\n";

  hartreeFockBrueckner(true);
}

//==============================================================================
std::vector<double> Wavefunction::vlocal(int l) const {
  return m_HF ? m_HF->vlocal(l) : m_vnuc;
}

//==============================================================================
std::vector<double> Wavefunction::Hmag(int l) const {
  return m_HF ? m_HF->Hmag(l) : std::vector<double>{};
}

//==============================================================================
// Number of electrons in the core
int Wavefunction::Ncore() const {
  auto count_electrons = [](int count, const DiracSpinor &Fc) {
    return count + Fc.num_electrons();
  };
  return std::accumulate(core().cbegin(), core().cend(), 0, count_electrons);
}

//==============================================================================
std::string Wavefunction::identity() const {
  const auto lab = m_run_label == "" ? "" : "_" + m_run_label;

  const std::string bq = vrad() && m_HF && m_HF->vBreit() ? "bq" :
                         m_HF && m_HF->vBreit()           ? "b" :
                         vrad()                           ? "q" :
                                                            "";
  const std::string va2 =
      dalpha2() > 1.0e-5 ? fmt::format("_{:1.4g}", 1.0 + dalpha2()) : "";
  return atomicSymbol() + std::to_string(Zion()) + bq + va2 + lab;
}

//==============================================================================
double Wavefunction::Hab(const DiracSpinor &Fa, const DiracSpinor &Fb) const {

  auto H0 = H0ab(Fa, Fb);

  const auto excl_exch = vHF() ? vHF()->excludeExchangeQ() : true;
  const auto sigmaQ = Sigma() != nullptr;
  const auto VBr = vHF() ? vHF()->vBreit() : nullptr;

  if (!excl_exch)
    H0 += Fa * vHF()->vexFa(Fb); // what about non-HF method?
  if (sigmaQ)
    H0 += Fa * (*Sigma())(Fb);
  if (VBr)
    H0 += Fa * VBr->VbrFa(Fb, core());
  return H0;
}

//==============================================================================
double Wavefunction::H0ab(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
  if (Fa.kappa() != Fb.kappa())
    return 0.0;
  const auto the_same = &Fa == &Fb;
  const auto dga = NumCalc::derivative(Fa.g(), rgrid->drdu(), rgrid->du(), 1);
  const auto dgb =
      the_same ? dga :
                 NumCalc::derivative(Fb.g(), rgrid->drdu(), rgrid->du(), 1);
  return H0ab_impl(Fa, dga, Fb, dgb);
}

//==============================================================================
double Wavefunction::H0ab(const DiracSpinor &Fa, const DiracSpinor &dFa,
                          const DiracSpinor &Fb, const DiracSpinor &dFb) const {
  return H0ab_impl(Fa, dFa.g(), Fb, dFb.g());
}

//==============================================================================
double Wavefunction::H0ab_impl(const DiracSpinor &Fa, std::vector<double> dga,
                               const DiracSpinor &Fb,
                               std::vector<double> dgb) const {
  // as above, but for when derivatives are already known
  if (Fa.kappa() != Fb.kappa())
    return 0.0;
  const auto kappa = Fa.kappa();
  const auto max = std::min(Fa.max_pt(), Fb.max_pt());
  const auto min = std::max(Fa.min_pt(), Fb.min_pt());
  const auto &drdu = Fa.grid().drdu();

  // auto dga = dFa.g();
  // auto dgb = dFb.g();

  for (std::size_t i = min; i < max; i++) {
    const auto r = Fa.grid().r(i);
    dga[i] -= (kappa * Fa.g(i) / r);
    dgb[i] -= (kappa * Fb.g(i) / r);
  }

  const auto D1m2 = NumCalc::integrate(1.0, min, max, Fa.f(), dgb, drdu) +
                    NumCalc::integrate(1.0, min, max, Fb.f(), dga, drdu);

  const auto Sab = NumCalc::integrate(1.0, min, max, Fa.g(), Fb.g(), drdu);

  const auto &v = vlocal(Fa.l());
  const auto Vab = NumCalc::integrate(1.0, min, max, Fa.f(), Fb.f(), v, drdu) +
                   NumCalc::integrate(1.0, min, max, Fa.g(), Fb.g(), v, drdu);

  const auto &Hmaga = Hmag(Fa.l());

  const auto H_mag =
      Hmaga.empty() ?
          0.0 :
          NumCalc::integrate(1.0, min, max, Fa.f(), Fb.g(), Hmaga, drdu) +
              NumCalc::integrate(1.0, min, max, Fa.g(), Fb.f(), Hmaga, drdu);
  const auto c = 1.0 / m_alpha;

  return (Vab - H_mag - c * (D1m2 + 2.0 * c * Sab)) * Fa.grid().du();
}

//==============================================================================
void Wavefunction::ConfigurationInteraction(const IO::InputBlock &input) {
  std::cout << "\n========================================================\n"
            << "Configuration Interaction:\n";
  IO::ChronoTimer t("CI");
  m_CIwfs = CI::configuration_interaction(input, *this);
  std::cout << "\n";
}

//==============================================================================
void Wavefunction::solve_exotic(const std::string &in_exotic_str, double mass,
                                bool print) {

  using namespace qip::overloads;

  const auto states = AtomData::listOfStates_nk(in_exotic_str);

  if (print) {
    std::cout << "\n---------------------------------------------\n";
    std::cout << "Exotic " << AtomData::atomicSymbol(this->Znuc()) << "\n";
    fmt::print("M = {:.8f} m_e = {:.12f} MeV\n\n", mass,
               mass * PhysConst::m_e_MeV);

    std::cout << "Energies - without screening:\n";
    std::cout << "nk    Rinf  eps    R_rms (a0)    E (au)            E "
                 "(keV)\n";
  }

  std::optional<DiracSpinor> Fscreen;
  for (const auto [n, kappa, x_en] : states) {

    // initial energy guess:
    const auto e0 =
        mass * AtomData::diracen(this->Znuc(), n, kappa, this->alpha());

    // nuclear potential (+QED). Note: only Uehling is really OK here.
    const auto v0 = this->vnuc() + (this->vrad() ? this->vrad()->Vel() :
                                                   std::vector<double>{});
    const auto Fnk = DiracODE::boundState(
        n, kappa, e0, this->grid_sptr(), v0, this->Hmag(), this->alpha(),
        1.0e-14, nullptr, nullptr, double(this->Znuc()), mass);

    const auto R_rms =
        std::sqrt(Fnk * (this->grid().r() * this->grid().r() * Fnk));

    this->valence().push_back(Fnk);

    if (!Fscreen) {
      Fscreen = Fnk;
    }

    if (print) {
      fmt::print("{:4s} {:5.2f}  {:5.0e}  {:.5e}  {:.9e}  {:.9e}\n",
                 Fnk.shortSymbol(), Fnk.rinf(), Fnk.eps(), R_rms, Fnk.en(),
                 Fnk.en() * PhysConst::Hartree_eV / 1.0e3);
    }
  }

  if (!this->core().empty() && Fscreen) {
    // make a _copy_ of the Hartree-Fock object for doing muon+electron HF
    // We include muon screening by temporarily updating Vnuc
    // (Can't add it to the core, since we must exclude exchange!)
    auto hf = *this->vHF();
    auto Fmu = *Fscreen;

    // Iterate screening of Hartree-Fock including muon
    if (print) {
      std::cout << "\nHartree-Fock for screening:\n";
      std::cout << "(Including the exotic " << Fmu.shortSymbol()
                << " state into direct part of HF)\n";
      fmt::print("it {:8s}  {:8s}  {:8s}\n", "HF_core", "electron", "exotic");
    }
    for (int it = 1; it < 100; ++it) {
      // Direct potential due to single muon:
      const auto Vmu = Coulomb::yk_ab(0, Fmu, Fmu);
      // Convenient: just add it to the nuclear potential:
      // Note: modify the copy of HF only. Be careful
      hf.vnuc() = this->vnuc() + Vmu;
      const auto E0 = hf.calculateCoreEnergy();
      const auto [elhf_eps, core_its, symbol] = hf.solve_core(false);
      const auto E1 = hf.calculateCoreEnergy();
      const auto eps_core = std::abs(E1 / E0 - 1.0);
      const auto v =
          hf.vdir() + this->vnuc() +
          (this->vrad() ? this->vrad()->Vel() : std::vector<double>{});
      const auto en0 = Fmu.en();
      DiracODE::boundState(Fmu, en0, v, this->Hmag(), this->alpha(), 1.0e-14,
                           nullptr, nullptr, double(this->Znuc()), mass);

      const auto eps_mu = std::abs(Fmu.en() / en0 - 1.0);

      if (print) {
        fmt::print("{:<2} {:.1e}   {:.1e}   {:.1e}\n", it, elhf_eps, eps_core,
                   eps_mu);
      }
      // check for convergance:
      if (it > 2 && std::max(eps_mu, eps_core) < 1.0e-13)
        break;
    }

    // Not really needed, but update HF in wavefunction:
    // Remember to not include muon in Vnuc
    hf.vnuc() = this->vnuc();
    *this->vHF() = hf;

    if (print) {
      std::cout << "\nCore (including screening by exotic "
                << Fscreen->shortSymbol() << "):\n";
      this->printCore();
      printf("E_c = %.6f\n\n", this->coreEnergyHF());
    }

    // Re-solve the muon states in the screened HF potential (just direct part)
    if (print) {
      std::cout << "Exotic energies - with screening:\n";
      std::cout << "nk    Rinf  eps    R_rms (a0)    E (au)            E "
                   "(keV)\n";
    }
    for (auto &Fnk : this->valence()) {

      const auto v =
          hf.vdir() + this->vnuc() +
          (this->vrad() ? this->vrad()->Vel() : std::vector<double>{});
      DiracODE::boundState(Fnk, Fnk.en(), this->vlocal(), this->Hmag(),
                           this->alpha(), 1.0e-14, nullptr, nullptr,
                           double(this->Znuc()), mass);

      if (print) {
        const auto R_rms =
            std::sqrt(Fnk * (this->grid().r() * this->grid().r() * Fnk));
        fmt::print("{:4s} {:5.2f}  {:5.0e}  {:.5e}  {:.9e}  {:.9e}\n",
                   Fnk.shortSymbol(), Fnk.rinf(), Fnk.eps(), R_rms, Fnk.en(),
                   Fnk.en() * PhysConst::Hartree_eV / 1.0e3);
      }
    }
  }
}