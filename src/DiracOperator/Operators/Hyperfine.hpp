#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "qip/Maths.hpp"
#include <functional>

namespace DiracOperator {

//******************************************************************************

//! Functions for F(r) [eg, nuclear magnetisation distribution] and similar
namespace Hyperfine {
using Func_R2_R = std::function<double(double, double)>; // save typing

inline auto sphericalBall_F() -> Func_R2_R {
  return [=](double r, double rN) {
    return (r > rN) ? 1.0 : (r * r * r) / (rN * rN * rN);
  };
}
inline auto sphericalShell_F() -> Func_R2_R {
  return [=](double r, double rN) { return (r > rN) ? 1.0 : 0.0; };
}
inline auto pointlike_F() -> Func_R2_R {
  return [=](double, double) { return 1.0; };
}

//! 'Volotka' single-particle model
inline auto volotkaBW_F(double mu, double I_nuc, double l_pn, int gl)
    -> Func_R2_R
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
  std::cout << "BW using: gl=" << g_l << ", gs=" << g_s << ", l=" << l_pn
            << ", gI=" << gI << " (I=" << I_nuc << ")\n";
  const double factor =
      (two_I == two_l + 1)
          ? g_s * (1 - two_I) / (4.0 * (two_I + 2)) + g_l * 0.5 * (two_I - 1)
          : g_s * (3 + two_I) / (4.0 * (two_I + 2)) +
                g_l * 0.5 * two_I * (two_I + 3) / (two_I + 2);
  if (two_I != two_l + 1 && two_I != two_l - 1) {
    std::cerr << "\nFAIL:59 in Hyperfine (volotkaBW_F):\n "
                 "we must have I = l +/- 1/2, but we have: I,l = "
              << I_nuc << "," << l_pn << "\n";
    return [](double, double) { return 0.0; };
  }
  return [=](double r, double rN) {
    return (r > rN) ? 1.0
                    : ((r * r * r) / (rN * rN * rN)) *
                          (1.0 - (3.0 / mu) * std::log(r / rN) * factor);
  };
}

//! 'Volotka' single-particle model, for doubly-odd nuclei
inline auto doublyOddBW_F(double mut, double It, double mu1, double I1,
                          double l1, int gl1, double I2, double l2) -> Func_R2_R
// F(r) * g = 0.5 [ g1F1 + g2F2 + (g1F1 - g2F2) * K]
// K = [I1(I1+1) - I2(I2+1)] / [I(I+1)]
// return F(r) [divide by g]
// volotkaBW_F(mu, I, l, g_l); //gl is 1 or 0
// g2 : from: g = 0.5 [ g1 + g2 + (g1 - g2) * K]
{
  const auto K = (I1 * (I1 + 1.0) - I2 * (I2 + 1.0)) / (It * (It + 1.0));
  const auto gt = mut / It;
  const auto g1 = mu1 / I1;
  const auto g2 = (g1 * (K + 1.0) - 2.0 * gt) / (K - 1.0);
  const auto mu2 = g2 * I2;
  const auto gl2 = (gl1 == 0) ? 1 : 0;
  const auto F1 = volotkaBW_F(mu1, I1, l1, gl1);
  const auto F2 = volotkaBW_F(mu2, I2, l2, gl2);
  return [=](double r, double rN) {
    return (0.5 / gt) * (g1 * F1(r, rN) + g2 * F2(r, rN) +
                         K * (g1 * F1(r, rN) - g2 * F2(r, rN)));
  };
}

inline std::vector<double> RadialFunc(int k, double rN, const Grid &rgrid,
                                      const Func_R2_R &hfs_F) {
  std::vector<double> rfunc;
  rfunc.reserve(rgrid.num_points());
  for (auto r : rgrid.r())
    rfunc.push_back(hfs_F(r, rN) / pow(r, k + 1));
  return rfunc;
}
} // namespace Hyperfine

//******************************************************************************
//******************************************************************************

// XXX Make special hfsA and hfsB operators, in terms of A and B coefs,
// and then general hfsK, just reduced matrix elements?

//******************************************************************************
//! @brief Magnetic hyperfine operator
//! @details Note: 'hfs_F' function (magnetization distribtuion) **includes**
//! the 1/r^2 (slightly different to definition in paper)
class HyperfineA final : public TensorOperator {

  using Func_R2_R = std::function<double(double, double)>; // save typing

public: // constructor
  HyperfineA(double muN, double IN, double rN, const Grid &rgrid,
             const Func_R2_R &hfs_F = Hyperfine::sphericalBall_F())
      : TensorOperator(1, Parity::even,
                       IN != 0.0 ? muN * PhysConst::muN_CGS_MHz / IN : 0.0,
                       Hyperfine::RadialFunc(1, rN, rgrid, hfs_F), 0)
  // , Inuc(IN)
  {
    if (IN == 0.0) {
      std::cout << "\nWarning: I=0 in Hyperfine operator; Setting gI to zero\n";
    }
  }
  std::string name() const override final { return "hfs"; }
  std::string units() const override final { return "MHz"; }

  double angularF(const int ka, const int kb) const override final {
    return (ka + kb) * Angular::Ck_kk(1, -ka, kb);
  }

  static double convertRMEtoA(const DiracSpinor &Fa, const DiracSpinor &Fb) {
    return 0.5 / Fa.jjp1() / Angular::Ck_kk(1, -Fa.k, Fb.k);
    // Correct for diag. Off diag? Prob not defined?
  }

  double hfsA(const DiracSpinor &Fa) const {
    auto Raa = radialIntegral(Fa, Fa);
    return Raa * Fa.k / (Fa.jjp1());
    // nb: in MHz
  }

  static double hfsA(const TensorOperator *h, const DiracSpinor &Fa) {
    auto Raa = h->radialIntegral(Fa, Fa);
    return Raa * Fa.k / (Fa.jjp1());
  }

  // XXX Make this a helper "conversion" function
  // double de_F(const DiracSpinor &Fa, double jF) const {
  //   auto Ahfs = hfsA(Fa); // nb: in MHz
  //   return 0.5 * Ahfs * (jF * (jF + 1.0) - Fa.jjp1() - Inuc * (Inuc + 1.0));
  // }

  double angularCff(int, int) const override final { return 0; }
  double angularCgg(int, int) const override final { return 0; }
  double angularCfg(int, int) const override final { return 1.0; }
  double angularCgf(int, int) const override final { return 1.0; }

private:
  // double Inuc;
};

//******************************************************************************
//! Units: Assumes g in nuc. magneton units (magnetic), and Q in barns
//! (electric)
class HyperfineK final : public TensorOperator {
  // see Xiao, ..., Derevianko, Phys. Rev. A 102, 022810 (2020).
  using Func_R2_R = std::function<double(double, double)>;

public:
  HyperfineK(int in_k, double in_GQ, double rN, const Grid &rgrid,
             const Func_R2_R &hfs_F = Hyperfine::sphericalBall_F())
      : TensorOperator(in_k, Parity::even, in_GQ,
                       Hyperfine::RadialFunc(in_k, rN, rgrid, hfs_F), 0),
        k(in_k),
        magnetic(k % 2 != 0),
        // gQ(in_GQ),
        cfg(magnetic ? 1.0 : 0.0),
        cff(magnetic ? 0.0 : 1.0) {}

  std::string name() const override final {
    return "hfs(" + std::to_string(k) + ")";
  }
  std::string units() const override final { return "MHz"; }

  double angularF(const int ka, const int kb) const override final {
    // inludes unit: Assumes g in nuc. magneton units, and/or Q in barns
    return magnetic
               ? (ka + kb) * Angular::Ck_kk(k, -ka, kb) * PhysConst::muN_CGS_MHz
               : -Angular::Ck_kk(k, ka, kb) * PhysConst::barn_MHz;

    // // For 'A' and 'B' constants:
    // auto j = Angular::j_k(kb);
    // return magnetic ? 0.5 * (ka + kb) / (j * (j + 1.0)) *
    // PhysConst::muN_CGS_MHz
    //                 : (2 * j - 1.0) / (2 * j + 2.0) * PhysConst::barn_MHz;
  }

  double angularCff(int, int) const override final { return cff; }
  double angularCgg(int, int) const override final { return cff; }
  double angularCfg(int, int) const override final { return cfg; }
  double angularCgf(int, int) const override final { return cfg; }

private:
  int k;
  bool magnetic;
  // double gQ;
  double cfg;
  double cff;
};

} // namespace DiracOperator
