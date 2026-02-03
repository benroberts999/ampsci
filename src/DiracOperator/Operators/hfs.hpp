#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "qip/Maths.hpp"
#include <functional>

namespace DiracOperator {

//==============================================================================
//! Functions for F(r) [nuclear magnetisation distribution] and similar
namespace Hyperfine {

using Func_R2_R = std::function<double(double, double)>; // save typing

//==============================================================================
//! Takes in F(r) and k, and forms hyperfine radial function: F(r,rN)/r^{k+1}
inline std::vector<double> RadialFunc(int k, double rN, const Grid &rgrid,
                                      const Func_R2_R &hfs_F) {
  std::vector<double> rfunc;
  rfunc.reserve(rgrid.num_points());
  for (const auto r : rgrid) {
    rfunc.push_back(hfs_F(r, rN) / std::pow(r, k + 1));
  }
  return rfunc;
}

//==============================================================================
//! Spherical ball F(r): (r/rN)^3 for r<rN, 1 for r>rN
inline auto sphericalBall_F(int k = 1) -> Func_R2_R {
  return [=](double r, double rN) {
    return (r > rN) ? 1.0 : std::pow(r / rN, 2 * k + 1);
  };
}

//! Spherical shell F(r): 0 for r<rN, 1 for r>rN
inline auto sphericalShell_F() -> Func_R2_R {
  return [=](double r, double rN) { return (r > rN) ? 1.0 : 0.0; };
}

//! Pointlike F(r): 1
inline auto pointlike_F() -> Func_R2_R {
  return [=](double, double) { return 1.0; };
}

//------------------------------------------------------------------------------
//! 'Volotka' single-particle model: see Phys. Rev. Lett. 125, 063002 (2020)
inline auto VolotkaSP_F(double mu, double I_nuc, double l_pn, int gl,
                        bool print = true) -> Func_R2_R
// Function that returns generates + returns F_BW Bohr-Weiskopf
// gl = 1 for proton, =0 for neutron. Double allowed for testing..
// mu is in units of Nuclear Magneton!
{
  const auto two_I = int(2 * I_nuc + 0.0001);
  const auto two_l = int(2 * l_pn + 0.0001);
  const auto g_l = double(gl); // just safety
  const auto gI = mu / I_nuc;

  const auto K = (l_pn * (l_pn + 1.0) - (3. / 4.)) / (I_nuc * (I_nuc + 1.0));
  const double g_s = (2.0 * gI - g_l * (K + 1.0)) / (1.0 - K);
  if (print)
    std::cout << "SingleParticle using: gl=" << g_l << ", gs=" << g_s
              << ", l=" << l_pn << ", gI=" << gI << " (I=" << I_nuc << ")\n";
  const double factor =
      (two_I == two_l + 1) ?
          g_s * (1 - two_I) / (4.0 * (two_I + 2)) + g_l * 0.5 * (two_I - 1) :
          g_s * (3 + two_I) / (4.0 * (two_I + 2)) +
              g_l * 0.5 * two_I * (two_I + 3) / (two_I + 2);
  if (two_I != two_l + 1 && two_I != two_l - 1) {
    std::cerr << "\nFAIL:59 in Hyperfine (VolotkaSP_F):\n "
                 "we must have I = l +/- 1/2, but we have: I,l = "
              << I_nuc << "," << l_pn << "\n";
    return [](double, double) { return 0.0; };
  }
  return [=](double r, double rN) {
    return (r > rN) ? 1.0 :
                      ((r * r * r) / (rN * rN * rN)) *
                          (1.0 - (3.0 / mu) * std::log(r / rN) * factor);
  };
}

//------------------------------------------------------------------------------
//! Elizarov et al., Optics and Spectroscopy (2006): u(r) = u0(R-r)^n
inline auto uSP(double mu, double I_nuc, double l_pn, int gl, double n,
                double R, bool u_option, bool print = true) -> Func_R2_R {
  const auto two_I = int(2 * I_nuc + 0.0001);
  const auto two_l = int(2 * l_pn + 0.0001);
  const auto g_l = double(gl); // just safety
  const auto gI = mu / I_nuc;

  const auto K = (l_pn * (l_pn + 1.0) - (3. / 4.)) / (I_nuc * (I_nuc + 1.0));
  const double g_s = (2.0 * gI - g_l * (K + 1.0)) / (1.0 - K);
  if (print) {
    std::cout << "uSP using: gl=" << g_l << ", gs=" << g_s << ", l=" << l_pn
              << ", gI=" << gI << " (I=" << I_nuc << ")\n";
    if (u_option) {
      std::cout << "u(r) = u1(r) = u0 (R-r)^n ";
    } else {
      std::cout << "u(r) = u2(r) = u0 r^n ";
    }
    std::cout << "with n = " << n << "\n";
  }
  if (two_I != two_l + 1 && two_I != two_l - 1) {
    std::cerr << "\nFAIL:59 in Hyperfine (uSP):\n "
                 "we must have I = l +/- 1/2, but we have: I,l = "
              << I_nuc << "," << l_pn << "\n";
    return [](double, double) { return 0.0; };
  }

  // nb: in theory, can read u(r) in from file!
  // r^2 * u(r)^2
  const auto r2u2 = [=](double r) {
    return u_option ? r * r * std::pow(R - r, 2.0 * n) :
                      r * r * std::pow(r, 2.0 * n);
  };

  // normalisation
  const auto u02 =
      1.0 / NumCalc::num_integrate(r2u2, 0.0, R, 100, NumCalc::linear);

  const auto F0r = [=](double r) {
    const auto f = [=](double x) { return u02 * r2u2(x); };
    return NumCalc::num_integrate(f, 0.0, r, 100, NumCalc::linear);
  };

  const auto FrR = [=](double r) {
    const auto f = [=](double x) {
      return u02 * r2u2(x) * r * r * r / x / x / x;
    };
    return NumCalc::num_integrate(f, r, R, 100, NumCalc::linear);
  };

  const auto two_Ip1 = 2.0 * (I_nuc + 1.0);
  const auto twoI_p3 = 2.0 * I_nuc + 3.0;
  const auto twoI_m1 = 2.0 * I_nuc - 1.0;

  double f1 = (two_I == two_l + 1) ?
                  0.5 * g_s + (I_nuc - 0.5) * g_l :
                  -I_nuc / two_Ip1 * g_s + I_nuc * twoI_p3 / two_Ip1 * g_l;

  double f2 =
      (two_I == two_l + 1) ?
          -twoI_m1 / (4.0 * two_Ip1) * g_s + (I_nuc - 0.5) * g_l :
          twoI_p3 / (4.0 * two_Ip1) * g_s + I_nuc * twoI_p3 / two_Ip1 * g_l;

  return [=](double r, double rN) {
    return (r > rN) ? 1.0 : (1.0 / mu) * (f1 * F0r(r) + f2 * FrR(r));
  };
}

//------------------------------------------------------------------------------
//! 'Volotka' SP model, for doubly-odd nuclei: Phys. Rev. Lett. 125, 063002 (2020)
inline auto doublyOddSP_F(double mut, double It, double mu1, double I1,
                          double l1, int gl1, double I2, double l2,
                          bool print = true) -> Func_R2_R
// F(r) * g = 0.5 [ g1F1 + g2F2 + (g1F1 - g2F2) * K]
// K = [I1(I1+1) - I2(I2+1)] / [I(I+1)]
// return F(r) [divide by g]
// VolotkaSP_F(mu, I, l, g_l); //gl is 1 or 0
// g2 : from: g = 0.5 [ g1 + g2 + (g1 - g2) * K]
{
  const auto K = (I1 * (I1 + 1.0) - I2 * (I2 + 1.0)) / (It * (It + 1.0));
  const auto gt = mut / It;
  const auto g1 = mu1 / I1;
  const auto g2 = (g1 * (K + 1.0) - 2.0 * gt) / (K - 1.0);
  const auto mu2 = g2 * I2;
  const auto gl2 = (gl1 == 0) ? 1 : 0;
  const auto F1 = VolotkaSP_F(mu1, I1, l1, gl1, print);
  const auto F2 = VolotkaSP_F(mu2, I2, l2, gl2, print);
  return [=](double r, double rN) {
    return (0.5 / gt) * (g1 * F1(r, rN) + g2 * F2(r, rN) +
                         K * (g1 * F1(r, rN) - g2 * F2(r, rN)));
  };
}

//------------------------------------------------------------------------------
//! Converts reduced matrix element to A/B coeficients (takes k, 2J, 2J)
inline double convert_RME_to_AB_2J(int k, int tja, int tjb) {
  // first, get stretched state:
  const auto tjz = std::min(tja, tjb); // arbitrary choice
  const auto s = Angular::neg1pow_2(tja - tjb);
  const auto tjs = Angular::threej_2(tja, 2 * k, tjb, -tjz, 0, tjz);
  // then: moment-specific factor
  // const auto f = k % 2 == 0 ? 2.0 : 1 / (0.5 * tjz);
  const auto f = [&]() {
    switch (k) {
    case 1: {
      // A: RME includes (μ/I)  ⇒  A = (1/J) <T1^e>_J (μ/I)
      const double J = 0.5 * tjz;
      return 1.0 / J;
    }
    case 2:
      // B: Q = 2 <T2^n>  ⇒  B = 2 Q <T2^e>_J
      return 2.0;
    case 3:
      // C: Ω = - <T3^n>  ⇒  C = -Ω <T3^e>_J
      return -1.0;
    case 4:
      // D: Π = <T4^n> (standard)  ⇒  D = Π <T4^e>_J
      return 1.0; // ???

    default:
      std::cout << "Warning: k=" << k
                << " has no standard A/B/C/D hyperfine-constant definition; "
                   "using f=1\n";
      return 1.0; // ???
    }
  }();
  return s * f * tjs;
}

//! Converts reduced matrix element to A/B coeficients (takes k, kappa, kappa)
inline double convert_RME_to_AB(int k, int ka, int kb) {
  const auto tja = Angular::twoj_k(ka);
  const auto tjb = Angular::twoj_k(kb);
  return convert_RME_to_AB_2J(k, tja, tjb);
}

inline double hfsA(const TensorOperator *h, const DiracSpinor &Fa) {
  auto Raa = h->radialIntegral(Fa, Fa); //?
  return Raa * Fa.kappa() / (Fa.jjp1()) * PhysConst::muN_CGS_MHz;
}

inline double hfsB(const TensorOperator *h, const DiracSpinor &Fa) {
  auto Raa = h->radialIntegral(Fa, Fa); //?
  return Raa * double(Fa.twoj() - 1) / double(Fa.twoj() + 2) *
         PhysConst::barn_MHz;
}

} // namespace Hyperfine

//==============================================================================
//==============================================================================
//==============================================================================
//! Units: Assumes nuclear moment in units of powers of nuclear magnetons and/or
//! barns - muN*b^(k-1)/2 for magnetic, and b^k/2 for electric
class hfs final : public TensorOperator {
  // see Xiao, ..., Derevianko, Phys. Rev. A 102, 022810 (2020).
  using Func_R2_R = std::function<double(double, double)>;

public:
  hfs(int in_k, double in_GQ, double rN_au, const Grid &rgrid,
      const Func_R2_R &hfs_F = Hyperfine::pointlike_F(), bool MHzQ = true)
      : TensorOperator(in_k, Parity::even, in_GQ,
                       Hyperfine::RadialFunc(in_k, rN_au, rgrid, hfs_F)),
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
                          Angular::Ck_kk(k, -ka, kb) * m_unit :
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
