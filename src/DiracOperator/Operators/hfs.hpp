#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "qip/Maths.hpp"
#include <functional>

namespace DiracOperator {

//==============================================================================
//! Auxillary Functions for hyperfine operatrs; F(r) [nuclear distribution] and similar
//! @details See \ref DiracOperator::hfs for operator
namespace Hyperfine {

//! Type for radial function F(r,rN) (type alias to save typing)
/*! @details
 - F(r,rN) contains only finite nuclear size correction, not HFS radial function.
 - F(r,rN) -> 1 for r > rN
 - r and rN always in atomic units
*/
using RadialFunction = std::function<double(double r, double r_Nuc)>;

//==============================================================================
/*! Forms the hyperfine radial function: F(r,rN)/r^{k+1}
  @details
  \f[
    t^k_{\rm radial}[r_i] = \frac{F(r_i,r_N)}{r_i^{k+1}}
  \f]

 @param k     Multipole rank (1 = M1, 2 = E2, ...)
 @param rN    Nuclear radius (a.u.); passed to F(r,r_N)
 @param r     Radial grid (a.u.)
 @param hfs_F Finite nuclear magnetisation function F(r,r_N)
 @return      Radial values of t^k(r)

 @note
 - Angular factors and signs are handled elsewhere.
 - F(r,r_N) contains only the finite nuclear-size correction.
 - Default F corresponds to the spherical-ball model.
 - F(r,r_N) -> 1 for r > r_N.
*/
std::vector<double> tk_radial(int k, double rN, const std::vector<double> &r,
                              const RadialFunction &hfs_F);

//==============================================================================
//! Spherical ball model for F(r,rN) [default]. Uniformly distributed point k-poles.
/*! @details
Based on a simple multipole expansion, assuming nucleus is made up from 
uniformly distributed point k-poles. Note: not everyone seems to define this 
the same way!

\f[
  F(r,r_N) = 
  \begin{cases}
    (r/r_N)^{2k+1} & r<r_N \\
    1 & r>r_N \\
  \end{cases}
\f]
*/
RadialFunction sphericalBall_F(int k);

//! Spherical shell F(r): 0 for r<rN, 1 for r>rN
RadialFunction sphericalShell_F();

//! Pointlike F(r): 1
RadialFunction pointlike_F();

//! Volotka single-particle nuclear model: F(r,rN)
/*! @details 
Calculates the BohrWeisskopf magnetisation distribution using the `Volotka' 
single-particle model of Phys. Rev. Lett. 125, 063002 (2020).
Returns a RadialFunction F(r,r_N).
F goes to 1 as r_N->0, and F=1 for r>rN.
Assumes radial nuclear density is a step function.
 
 @param mu     Nuclear magnetic moment (in nuclear magnetons)
 @param I_nuc  Nuclear spin
 @param l_pn   Orbital angular momentum of valence nucleon
 @param gl     Orbital g-factor (1 = proton, 0 = neutron)
 @param print  Print model details
 
 @return Radial function \f$ F_{\mathrm{BW}}(r, r_N) \f$
*/
RadialFunction VolotkaSP_F(double mu, double I_nuc, double l_pn, int gl,
                           bool print = true);

//! Elizarov single-particle magnetisation model [extended Volotka]
/*! @details
Returns the radial magnetisation function F(r,r_N) for the model of
Elizarov, A. A. et al., Opt. Spectrosc. 100, 361 (2006).
Uses nuclear radial density u(r), with:

\f[
  u(r) = 
    \begin{cases}
      u_0 (R-r)^n & \text{u1(r)} \\
      u_0 r^n & \text{u2(r)}
    \end{cases}
\f]

n=0 should correspond to Volotka model.

  @param mu     Nuclear magnetic moment (in nuclear magnetons)
  @param I_nuc  Nuclear spin
  @param l_pn   Orbital angular momentum of valence nucleon
  @param gl     Orbital g-factor (1 = proton, 0 = neutron)
  @param n      Power in usually 0,1,2, but can be anything. 0 => Volotka
  @param R      Nuclear radius
  @param u_option Selects returned radial form: true means u1(r)
  @param print  Print model details
  @return Radial magnetisation function
  @warning Does normalisation by numerical integration; may be unstable. 
  Check that returns to Volotka as n->0
*/
RadialFunction uSP(double mu, double I_nuc, double l_pn, int gl, double n,
                   double R, bool u_option, bool print = true);

//! Volotka single-particle model for doubly-odd nuclei
/*! @details
  From: [Phys. Rev. Lett. 125, 063002 (2020).](http://arxiv.org/abs/2001.01907)

  Total magnetisation:
  \f[
    g F(r) = 0.5 \left[ g_1 F_1(r) + g_2 F_2(r) 
                        + (g_1 F_1(r) - g_2 F_2(r)) K \right]
  \f]

  with

  \f[
    K = \frac{[ I_1(I_1+1) - I_2(I_2+1) ]}{[ I(I+1) ]}
  \f]

  Returns F(r) (i.e., divided by total g).

  g2 is obtained from
  g = 0.5 [ g1 + g2 + (g1 - g2) K ].

  @note F1, F2 are the Volotka single-particle functions (see \ref VolotkaSP_F).

  @param mut  Total nuclear magnetic moment (in nuclear magnetons)
  @param It   Total nuclear spin
  @param mu1  Magnetic moment of nucleon 1
  @param I1   Spin of nucleon 1
  @param l1   Orbital angular momentum of nucleon 1
  @param gl1  Orbital g-factor of nucleon 1 (1 = proton, 0 = neutron)
  @param I2   Spin of nucleon 2
  @param l2   Orbital angular momentum of nucleon 2
  @param print Print model details

  @return Radial function F(r)
*/
RadialFunction doublyOddSP_F(double mut, double It, double mu1, double I1,
                             double l1, int gl1, double I2, double l2,
                             bool print = true);

//------------------------------------------------------------------------------
//! Converts reduced matrix element to A/B coeficients (takes k, 2J, 2J)
double convert_RME_to_HFSconstant_2J(int k, int tja, int tjb);

//! Converts reduced matrix element to A/B coeficients (takes k, kappa, kappa)
double convert_RME_to_HFSconstant(int k, int ka, int kb);

} // namespace Hyperfine

//==============================================================================
//==============================================================================
//==============================================================================

//! Generalised hyperfine-structure operator, including relevant nuclear moment
/*! @details

Implements the nuclear multipole hyperfine operator of rank @p k using a
specified finite-nucleus radial model.

By default, includes the nuclear moment, and relevant factors.
Input parameter @p GQ is the g-factor (k=1) or nuclear moment (k>1), 
in units of nuclear magnetons * barns^power.

That is, it calculates:

\f[
   GQ * t^k
\f]

with:

\f[
    t^k = 
    \frac{-1}{r^{k+1}}\,F(r)
    \begin{cases}
      C^k & \text{electric (even k)} \\
      \sqrt{\frac{k+1}{k}}
        \mathbf{\alpha}\!\cdot\!\mathbf{C}^{(0)}_k 
        & \text{magnetic (odd k)}
    \end{cases}
\f]

where F(r) is nuclear k-pole distribution (=1 for pointlike).

Convert to hyperfine constant by taking stretched state, and multiply by:

\f[
  M = 
  \begin{cases}
    1/J      &  k = 1 \\
    2        &  k = 2 \\
    -1       &  k = 3 \\
    1        &  k \geq 4 \\
  \end{cases}
\f]

This factor (including the steched 3j symbol) is returned by \ref Hyperfine::convert_RME_to_HFSconstant()

Units: 
- Assumes nuclear moments in units of:
- magnetic moments: @f$\mu_N\, b^{(k-1)/2}@f$
- electric moments: @f$b^{k/2}@f$
- Matrix element are in MHz by default, otherwise in atomic units.

Here \f$ \mu_N \f$ is the nuclear magneton and \f$ b \f$ is the barn.

See, e.g., Xiao _et al._, [Phys. Rev. A 102, 022810 (2020).](http://arxiv.org/abs/2007.06798)

The radial part is constructed from \ref Hyperfine::tk_radial().
*/
class hfs final : public TensorOperator {
  using RadialFunction = std::function<double(double, double)>;

public:
  //! @brief Constructs hyperfine operator of multipolarity k.
  /*! @details
    Initialises the hyperfine interaction operator.

    @param in_k     Tensor rank k (multipolarity) of the hyperfine interaction 
                    (e.g., 1,2,3,... for M1, E2, M3,... etc.)
    @param GQ       Nuclear g-factor for k=1, general moment for k>1
    @param rN_au    Nuclear radius (a.u.), used for finite-size magnetisation models.
    @param rgrid    Radial grid used to define the radial functions.
    @param hfs_F    Radial magnetisation distribution function
                    (see @ref DiracOperator::Hyperfine).
    @param MHzQ     If true outputs in MHz units; otherwise in atomic units.

    @see See also: \ref DiracOperator::Hyperfine
  */
  hfs(int in_k, double GQ, double rN_au, const Grid &rgrid,
      const RadialFunction &hfs_F = Hyperfine::pointlike_F(), bool MHzQ = true)
    : TensorOperator(in_k, Parity::even, GQ,
                     Hyperfine::tk_radial(in_k, rN_au, rgrid.r(), hfs_F)),
      k(in_k),
      magnetic(k % 2 != 0),
      cfg(magnetic ? 1.0 : 0.0),
      cff(magnetic ? 0.0 : 1.0),
      mMHzQ(MHzQ) {

    const auto power = magnetic ? (k - 1) / 2 : k / 2;

    // Assumes nuclear moment in muN*b^(k-1)/2 for magnetic,
    // and b^k/2 for electric
    const auto unit_au =
      magnetic ? PhysConst::muN_CGS * std::pow(PhysConst::barn_au, power) :
                 std::pow(PhysConst::barn_au, power);
    m_unit = mMHzQ ? unit_au * PhysConst::Hartree_MHz : unit_au;
  }

  std::string name() const override final {
    return "hfs" + std::to_string(k) + "";
  }
  std::string units() const override final { return mMHzQ ? "MHz" : "au"; }

  double angularF(const int ka, const int kb) const override final {
    return magnetic ? -double(ka + kb) / double(k) *
                        Angular::Ck_kk(k, ka, -kb) * m_unit :
                      -Angular::Ck_kk(k, ka, kb) * m_unit;
  }

  double angularCff(int, int) const override final { return cff; }
  double angularCgg(int, int) const override final { return cff; }
  double angularCfg(int, int) const override final { return cfg; }
  double angularCgf(int, int) const override final { return cfg; }

private:
  int k;
  bool magnetic;
  double cfg;
  double cff;
  bool mMHzQ;
  double m_unit{1.0};
};

//==============================================================================
//==============================================================================
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_hfs(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;

  input.check(
    {{"", "Most following will be taken from the default nucleus if "
          "not explicitely given"},
     {"mu", "Magnetic moment in mu_N"},
     {"Q", "Nuclear quadrupole moment, in barns. Also used as overall "
           "constant for any higher-order moments [1.0]"},
     {"k", "Multipolarity. 1=mag. dipole, 2=elec. quad, etc. [1]"},
     {"rrms",
      "nuclear (magnetic) rms radius, in Fermi (fm) (defult is charge rms)"},
     {"units", "Units for output (only for k=1,k=2). MHz or au [MHz]"},
     {"F", "F(r): Nuclear moment distribution: ball, point, shell, "
           "SingleParticle, or doublyOddSP [ball]"},
     {"F(r)", "Obselete; use 'F' from now - will be removed"},
     {"nuc_mag", "Obselete; use 'F' from now - will be removed"},
     {"printF", "Writes F(r) to a text file [false]"},
     {"print", "Write F(r) info to screen [true]"},
     {"", "The following are only for F=SingleParticle or doublyOddSP"},
     {"I", "Nuclear spin. Taken from nucleus"},
     {"parity", "Nulcear parity: +/-1"},
     {"l", "l for unpaired nucleon (automatically derived from I and "
           "parity; best to leave as default)"},
     {"gl", "=1 for proton, =0 for neutron"},
     {"", "The following are only used if F=doublyOddSP"},
     {"mu1", "mag moment of 'first' unpaired nucleon"},
     {"gl1", "gl of 'first' unpaired nucleon"},
     {"l1", "l of 'first' unpaired nucleon"},
     {"l2", "l of 'second' unpaired nucleon"},
     {"I1", "total spin (J) of 'first' unpaired nucleon"},
     {"I2", "total spin (J) of 'second' unpaired nucleon"},
     {"", "The following are only for u(r) function"},
     {"u", "u1 or u2 : u1(r) = (R-r)^n, u2(r) = r^n [u1]"},
     {"n", "n that appears above. Should be between 0 and 2 [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }

  const auto nuc = wf.nucleus();
  const auto isotope = Nuclear::findIsotopeData(nuc.z(), nuc.a());
  auto mu = input.get("mu", isotope.mu ? *isotope.mu : 1.0);
  auto I_nuc = input.get("I", isotope.I_N ? *isotope.I_N : 1.0);
  const auto print = input.get("print", true);
  const auto k = input.get("k", 1);

  const auto use_MHz =
    qip::ci_compare(input.get<std::string>("units", "MHz"), "MHz");

  if (k <= 0) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 246:\n");
    std::cout << "In hyperfine: invalid K=" << k << "! meaningless results\n";
  }
  if (I_nuc <= 0 && k == 1) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning 253:\n");
    std::cout << "In hyperfine: invalid I_nuc=" << I_nuc
              << "! meaningless results\nSetting I=1\n";
    I_nuc = 1;
  }
  if (mu == 0.0 && k == 1) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning 352:\n");
    std::cout << "Setting mu=1\n";
    mu = 1;
  }

  const auto g_or_Q = (k == 1) ? (mu / I_nuc) :
                      (k == 2) ? input.get("Q", isotope.q ? *isotope.q : 1.0) :
                                 input.get("Q", 1.0);

  enum class DistroType {
    point,
    ball,
    shell,
    SingleParticle,
    doublyOddSP,
    spu,
    Error
  };

  std::string default_distribution = "ball";

  // For compatability with old notation of 'F(r)' input option
  const auto Fr_str =
    input.has_option("F") ?
      input.get<std::string>("F", default_distribution) :
    input.has_option("nuc_mag") ?
      input.get<std::string>("nuc_mag", default_distribution) :
      input.get<std::string>("F(r)", default_distribution);

  const auto distro_type =
    (qip::ci_wc_compare(Fr_str, "point*") || qip::ci_compare(Fr_str, "1")) ?
      DistroType::point :
    qip::ci_compare(Fr_str, "ball")          ? DistroType::ball :
    qip::ci_compare(Fr_str, "shell")         ? DistroType::shell :
    qip::ci_wc_compare(Fr_str, "Single*")    ? DistroType::SingleParticle :
    qip::ci_wc_compare(Fr_str, "doublyOdd*") ? DistroType::doublyOddSP :
    qip::ci_compare(Fr_str, "spu")           ? DistroType::spu :
                                               DistroType::Error;
  if (distro_type == DistroType::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 271:\n");
    std::cout << "\nIn hyperfine. Unkown F(r) - " << Fr_str << "\n";
    std::cout << "Defaulting to pointlike!\n";
  }

  const auto r_rmsfm =
    distro_type == DistroType::point ? 0.0 : input.get("rrms", nuc.r_rms());
  const auto r_nucfm = std::sqrt(5.0 / 3) * r_rmsfm;
  const auto r_nucau = r_nucfm / PhysConst::aB_fm;

  if (print) {
    std::cout << "\nHyperfine structure: " << wf.atom() << "\n";
    std::cout << "K=" << k << " ("
              << (k == 1     ? "magnetic dipole" :
                  k == 2     ? "electric quadrupole" :
                  k % 2 == 0 ? "electric multipole" :
                               "magnetic multipole")
              << ")\n";
    std::cout << "Using " << Fr_str << " nuclear distro for F(r)\n"
              << "w/ r_N = " << r_nucfm << "fm = " << r_nucau
              << "au  (r_rms=" << r_rmsfm << "fm)\n";
    std::cout << "Points inside nucleus: " << wf.grid().getIndex(r_nucau)
              << "\n";
    if (k == 1) {
      std::cout << "mu = " << mu << ", I = " << I_nuc << ", g = " << g_or_Q
                << "\n";
    } else {
      std::cout << "Q = " << g_or_Q << "\n";
    }
  }

  // default is BALL:
  auto Fr = Hyperfine::sphericalBall_F(k);
  if (distro_type == DistroType::ball) {
    Fr = Hyperfine::sphericalBall_F(k);
  } else if (distro_type == DistroType::shell) {
    Fr = Hyperfine::sphericalShell_F();
  } else if (distro_type == DistroType::SingleParticle) {
    const auto pi = input.get("parity", isotope.parity ? *isotope.parity : 0);
    const auto l_tmp = int(I_nuc + 0.5 + 0.0001);
    auto l = ((l_tmp % 2 == 0) == (pi == 1)) ? l_tmp : l_tmp - 1;
    l = input.get("l", l); // can override derived 'l' (not recommended)
    const auto gl_default = wf.Znuc() % 2 == 0 ? 0 : 1; // unparied proton?
    const auto gl = input.get<int>("gl", gl_default);
    if (print) {
      std::cout << "Single-Particle (Volotka formula) for unpaired";
      if (gl == 1)
        std::cout << " proton ";
      else if (gl == 0)
        std::cout << " neturon ";
      else
        std::cout << " gl=" << gl << "??? program will run, but prob wrong!\n";
      std::cout << "with l=" << l << " (pi=" << pi << ")\n";
    }
    Fr = Hyperfine::VolotkaSP_F(mu, I_nuc, l, gl, print);
  } else if (distro_type == DistroType::spu) {
    const auto pi = input.get("parity", isotope.parity ? *isotope.parity : 0);
    const auto u_func = input.get("u", std::string{"u1"}); // u1=(R-r)^n, u2=r^n
    const bool u_option = u_func == std::string{"u1"};
    const auto n = input.get("n", 0.0); // u1=(R-r)^n, u2=r^n
    const auto l_tmp = int(I_nuc + 0.5 + 0.0001);
    auto l = ((l_tmp % 2 == 0) == (pi == 1)) ? l_tmp : l_tmp - 1;
    l = input.get("l", l); // can override derived 'l' (not recommended)
    const auto gl_default = wf.Znuc() % 2 == 0 ? 0 : 1; // unparied proton?
    const auto gl = input.get<int>("gl", gl_default);
    if (print) {
      std::cout << "Single-Particle (Volotka formula) with u(r) for unpaired";
      if (gl == 1)
        std::cout << " proton ";
      else if (gl == 0)
        std::cout << " neturon ";
      else
        std::cout << " gl=" << gl << "??? program will run, but prob wrong!\n";
      std::cout << "with l=" << l << " (pi=" << pi << ")\n";
    }
    Fr = Hyperfine::uSP(mu, I_nuc, l, gl, n, r_nucau, u_option, print);
  } else if (distro_type == DistroType::doublyOddSP) {
    const auto mu1 = input.get<double>("mu1", 1.0);
    const auto gl1 = input.get<int>("gl1", -1); // 1 or 0 (p or n)
    if (gl1 != 0 && gl1 != 1) {
      fmt2::styled_print(fg(fmt::color::red), "\nError 324:\n");
      std::cout << "In " << input.name() << " " << Fr_str
                << "; have gl1=" << gl1 << " but need 1 or 0\n";
      return std::make_unique<NullOperator>(NullOperator());
    }
    const auto l1 = input.get<double>("l1", -1.0);
    const auto l2 = input.get<double>("l2", -1.0);
    const auto I1 = input.get<double>("I1", -1.0);
    const auto I2 = input.get<double>("I2", -1.0);

    Fr = Hyperfine::doublyOddSP_F(mu, I_nuc, mu1, I1, l1, gl1, I2, l2, print);
  }

  // Optionally print F(r) function to file
  if (input.get<bool>("printF", false)) {
    std::ofstream of(wf.identity() + "_" + Fr_str + ".txt");
    of << "r/fm  F(r)\n";
    for (auto r : wf.grid()) {
      of << r * PhysConst::aB_fm << " "
         << Fr(r * PhysConst::aB_fm, r_nucau * PhysConst::aB_fm) << "\n";
    }
  }

  return std::make_unique<hfs>(k, g_or_Q, r_nucau, wf.grid(), Fr, use_MHz);
}

} // namespace DiracOperator
