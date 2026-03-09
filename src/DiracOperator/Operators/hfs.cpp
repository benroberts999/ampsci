#include "hfs.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "qip/Maths.hpp"
#include <functional>

namespace DiracOperator {
namespace Hyperfine {

//==============================================================================
// Takes in F(r) and k, and forms hyperfine radial function: F(r,rN)/r^{k+1}
std::vector<double> tk_radial(int k, double rN, const std::vector<double> &r,
                              const RadialFunction &hfs_F) {
  std::vector<double> rfunc;
  rfunc.reserve(r.size());
  for (const auto r_i : r) {
    rfunc.push_back(hfs_F(r_i, rN) / qip::pow(r_i, k + 1));
  }
  return rfunc;
}

//==============================================================================
// Spherical ball F(r): (r/rN)^3 for r<rN, 1 for r>rN
RadialFunction sphericalBall_F(int k) {
  return [=](double r, double rN) {
    return (r > rN) ? 1.0 : std::pow(r / rN, 2 * k + 1);
  };
}

// Spherical shell F(r): 0 for r<rN, 1 for r>rN
RadialFunction sphericalShell_F() {
  return [=](double r, double rN) { return (r > rN) ? 1.0 : 0.0; };
}

// Pointlike F(r): 1
RadialFunction pointlike_F() {
  return [=](double, double) { return 1.0; };
}

//------------------------------------------------------------------------------
// 'Volotka' single-particle model: see Phys. Rev. Lett. 125, 063002 (2020)
RadialFunction VolotkaSP_F(double mu, double I_nuc, double l_pn, int gl,
                           bool print)

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
// Elizarov et al., Optics and Spectroscopy (2006): u(r) = u0(R-r)^n
RadialFunction uSP(double mu, double I_nuc, double l_pn, int gl, double n,
                   double R, bool u_option, bool print) {
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
// 'Volotka' SP model, for doubly-odd nuclei: Phys. Rev. Lett. 125, 063002 (2020)
RadialFunction doublyOddSP_F(double mut, double It, double mu1, double I1,
                             double l1, int gl1, double I2, double l2,
                             bool print) {
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
// Converts reduced matrix element to A/B coeficients (takes k, 2J, 2J)
double convert_RME_to_HFSconstant_2J(int k, int tja, int tjb) {
  // first, get stretched state:
  const auto tjz = std::min(tja, tjb); // arbitrary choice
  const auto s = Angular::neg1pow_2(tja - tjb);
  const auto tjs = Angular::threej_2(tja, 2 * k, tjb, -tjz, 0, tjz);
  // then: moment-specific factor
  // const auto f = k % 2 == 0 ? 2.0 : 1 / (0.5 * tjz);
  const auto f = [&]() {
    switch (k) {
    case 1: {
      // A: RME includes (μ/I) ⇒ A = (1/J) <T1^e>_J (μ/I)
      const double J = 0.5 * tjz;
      return 1.0 / J;
    }
    case 2:
      // B: Q = 2 <T2^n> ⇒ B = 2 Q <T2^e>_J
      return 2.0;
    case 3:
      // C: Ω = - <T3^n> ⇒ C = -Ω <T3^e>_J
      return -1.0;
    case 4:
      // D: Π = <T4^n> (standard) ⇒ D = Π <T4^e>_J
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

// Converts reduced matrix element to A/B coeficients (takes k, kappa, kappa)
double convert_RME_to_HFSconstant(int k, int ka, int kb) {
  const auto tja = Angular::twoj_k(ka);
  const auto tjb = Angular::twoj_k(kb);
  return convert_RME_to_HFSconstant_2J(k, tja, tjb);
}

} // namespace Hyperfine
} // namespace DiracOperator