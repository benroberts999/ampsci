#pragma once
#include "Dirac/DiracOperator.hpp"
#include "Maths/Grid.hpp"
// #include "Maths/SphericalBessel.hpp"
#include "Physics/Nuclear.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/Wigner_369j.hpp"
#include <cmath>
#include <functional>
#include <vector>

// XXX move impl into cpp !! XXX

// Classes (inherit from DriacOperator)

//******************************************************************************
class RadialFuncOperator : public ScalarOperator {
  // = some function of r
  // Pass only grid to just get r, or
  // either pass a lambda/function [f(r)], or a number, n, (for r^n)
  // Explicitely even with rank 0, so ...
public:
  std::vector<double> fillVec(const Grid &gr,
                              const std::function<double(double)> &f) {
    std::vector<double> f_r;
    f_r.reserve(gr.ngp);
    for (auto r : gr.r)
      f_r.push_back(f(r));
    return f_r;
  }

public:
  RadialFuncOperator(const Grid &rgrid, const std::function<double(double)> &f)
      : ScalarOperator(OperatorParity::even, 1.0, fillVec(rgrid, f)) {}
  RadialFuncOperator(const Grid &rgrid, const double n)
      : ScalarOperator(OperatorParity::even, 1.0, fillVec(rgrid, [n](double r) {
                         return std::pow(r, n);
                       })) {}
};

//******************************************************************************
class E1Operator : public DiracOperator
// v{d} = -e v{r}   [e=|e|=1]
// <a||d||b> = R (-1)^{ja+1/2} Sqrt([ja][jb]) tjs(ja,jb,1, -1/2,1/2,0)
//           = R <ja||C^k||jb>
// R = -e Int[ r(fafb + gagb) ]dr
{
public:
  E1Operator(const Grid &gr)
      : DiracOperator(1, OperatorParity::odd, -1.0, gr.r, 0) {}

  double reducedME(const DiracSpinor &Fa,
                   const DiracSpinor &Fb) const override {
    auto Rab = radialIntegral(Fa, Fb);
    // move to lookup table?
    auto Cab = Wigner::Ck_2j2j(1, Fa.twoj(), Fb.twoj());
    return Rab * Cab;
  }
};

//******************************************************************************
class E1Operator_vform : public DiracOperator
// d_v = (ie/w alpha) v{alpha}   [v{a} = g0v{g}]
// <a||dv||b> = -2e/(w alpha) Int[ fagb <ka||s||-kb> - gafb <-ka||s||kb>]
{
public:
  E1Operator_vform(const Grid &gr, const double alpha = PhysConst::alpha)
      : DiracOperator(1, OperatorParity::odd, -1.0,
                      std::vector<double>(gr.ngp, 1.0), 0),
        m_c(1.0 / alpha) {}

  double reducedME(const DiracSpinor &Fa,
                   const DiracSpinor &Fb) const override {
    if (Fa.k == Fb.k)
      return 0;
    auto Rab = radialIntegral(Fa, Fb);
    auto omega = Fa.en - Fb.en;
    return 2.0 * m_c * Rab / omega; // includes var-alpha
  }

private:
  virtual double angularCff(int, int) const override { return 0; }
  virtual double angularCgg(int, int) const override { return 0; }
  virtual double angularCfg(int ka, int kb) const override {
    return Wigner::S_kk(ka, -kb);
  }
  virtual double angularCgf(int ka, int kb) const override {
    return -Wigner::S_kk(-ka, kb);
  }

private:
  double m_c; // speed of light (including var-alpha)
};

// //******************************************************************************
// class M1Operator : public DiracOperator
// // M1 = mu . B
// // <a||M1||b> = 3 (ka + kb) <-ka||C^1||kb> Int [(fagb+gafb)j_1(kr)]
// // k = w/c, j1 is first spherical bessel
// // Two options; a) different operator for each w
// // b) generate jL each call
// {
// public:
//   M1Operator(const Grid &gr, double omega)
//       : DiracOperator(1, OperatorParity::even, 3.0,
//                       SphericalBessel::fillBesselVec(1, set_kr(omega, gr)),
//                       0) {
//   }
//
//   double reducedME(const DiracSpinor &Fa,
//                    const DiracSpinor &Fb) const override {
//     auto Rab = radialIntegral(Fa, Fb);
//     auto Cab = (Fa.k + Fb.k) * Wigner::Ck_kk(1, -Fa.k, Fb.k);
//     return Rab * Cab / PhysConst::alpha / PhysConst::muB_CGS; //???
//   }
//
// private:
//   virtual double angularCff(int, int) const override { return 0; }
//   virtual double angularCgg(int, int) const override { return 0; }
//   virtual double angularCfg(int, int) const override { return 1.0; }
//   virtual double angularCgf(int, int) const override { return 1.0; }
//
// private:
//   // std::vector<double> kr;
//   std::vector<double> set_kr(double omega, const Grid &gr) {
//     std::vector<double> kr;
//     kr.reserve(gr.ngp);
//     for (const auto &r : gr.r) {
//       kr.push_back(r * omega * PhysConst::alpha); // this? or var-alpha?
//     }
//     return kr;
//   }
// };

//******************************************************************************
class HyperfineOperator : public DiracOperator {

  using Func_R2_R = std::function<double(double, double)>; // save typing

public: // F(r) functions. NOTE: includes 1/r^2 !
  static inline auto sphericalBall_F() -> Func_R2_R {
    return [=](double r, double rN) {
      return (r > rN) ? 1. / (r * r) : r / (rN * rN * rN);
    };
  }
  static inline auto sphericalShell_F() -> Func_R2_R {
    return [=](double r, double rN) { return (r > rN) ? 1. / (r * r) : 0.0; };
  }
  static inline auto pointlike_F() -> Func_R2_R {
    return [=](double r, double) { return 1. / (r * r); };
  }

  static inline auto generate_F_BW(double mu, double I_nuc, double l_pn, int gl)
      -> Func_R2_R
  // Function that returns generates + returns F_BW Bohr-Weiskopf
  // gl = 1 for proton, =0 for neutron. Double allowed for testing..
  // mu is in units of Nuclear Magneton!
  {
    const auto two_I = int(2 * I_nuc + 0.0001);
    const auto two_l = int(2 * l_pn + 0.0001);
    auto g_l = double(gl); // just safety
    auto gI = mu / I_nuc;
    auto K = (l_pn * (l_pn + 1.0) - (3. / 4.)) / (I_nuc * (I_nuc + 1.0));
    const double g_s = (2.0 * gI - g_l * (K + 1.0)) / (1.0 - K);
    std::cout << "BW using: gl=" << g_l << ", gs=" << g_s << ", l=" << l_pn
              << ", gI=" << gI << " (I=" << I_nuc << ")\n";
    const double factor =
        (two_I == two_l + 1)
            ? g_s * (1 - two_I) / (4.0 * (two_I + 2)) + g_l * 0.5 * (two_I - 1)
            : g_s * (3 + two_I) / (4.0 * (two_I + 2)) +
                  g_l * 0.5 * two_I * (two_I + 3) / (two_I + 2);
    if (two_I != two_l + 1 && two_I != two_l - 1) {
      std::cerr << "\nFAIL:59 in HyperfineOperator (generate_F_BW):\n "
                   "we must have I = l +/- 1/2, but we have: I,l = "
                << I_nuc << "," << l_pn << "\n";
      return [](double, double) { return 0.0; };
    }
    return [=](double r, double rN) {
      return (r > rN) ? 1. / (r * r)
                      : (r / (rN * rN * rN)) *
                            (1.0 - (3.0 / mu) * std::log(r / rN) * factor);
    };
  }

  static inline auto generate_F_BW_doublyOdd(double mut, double It, double mu1,
                                             double I1, double l1, int gl1,
                                             double I2, double l2) -> Func_R2_R
  // F(r) * g = 0.5 [ g1F1 + g2F2 + (g1F1 - g2F2) * K]
  // K = [I1(I1+1) - I2(I2+1)] / [I(I+1)]
  // return F(r) [divide by g]
  // generate_F_BW(mu, I, l, g_l); //gl is 1 or 0
  // g2 : from: g = 0.5 [ g1 + g2 + (g1 - g2) * K]
  {
    auto K = (I1 * (I1 + 1.0) - I2 * (I2 + 1.0)) / (It * (It + 1.0));
    auto gt = mut / It;
    auto g1 = mu1 / I1;
    auto g2 = (g1 * (K + 1.0) - 2.0 * gt) / (K - 1.0);
    auto mu2 = g2 * I2;
    auto gl2 = (gl1 == 0) ? 1 : 0;
    auto F1 = generate_F_BW(mu1, I1, l1, gl1);
    auto F2 = generate_F_BW(mu2, I2, l2, gl2);
    return [=](double r, double rN) {
      return (0.5 / gt) * (g1 * F1(r, rN) + g2 * F2(r, rN) +
                           K * (g1 * F1(r, rN) - g2 * F2(r, rN)));
    };
  }

private: // helper
  static inline std::vector<double> RadialFunc(double rN, const Grid &rgrid,
                                               const Func_R2_R &hfs_F) {
    std::vector<double> rfunc;
    rfunc.reserve(rgrid.ngp);
    for (auto r : rgrid.r)
      rfunc.push_back(hfs_F(r, rN));
    return rfunc;
  }

public: // constructor
  HyperfineOperator(double muN, double IN, double rN, const Grid &rgrid,
                    const Func_R2_R &hfs_F = sphericalBall_F())
      : DiracOperator(1, OperatorParity::even, muN * PhysConst::muN_CGS / IN,
                      RadialFunc(rN, rgrid, hfs_F), 0),
        Inuc(IN) {}

  double reducedME(const DiracSpinor &Fa,
                   const DiracSpinor &Fb) const override {
    auto Rab = radialIntegral(Fa, Fb);
    return Rab * (Fa.k + Fb.k) * Wigner::Ck_kk(1, -Fa.k, Fb.k);
  }

  double hfsA(const DiracSpinor &Fa) {
    auto Raa = radialIntegral(Fa, Fa);
    return PhysConst::Hartree_MHz * Raa * Fa.k / (Fa.jjp1());
    // nb: in MHz
  }
  static double hfsA(const DiracOperator *h, const DiracSpinor &Fa) {
    auto Raa = h->radialIntegral(Fa, Fa);
    return PhysConst::Hartree_MHz * Raa * Fa.k / (Fa.jjp1());
  }

  double de_F(const DiracSpinor &Fa, double jF) {
    auto Ahfs = hfsA(Fa); // nb: in MHz
    return 0.5 * Ahfs * (jF * (jF + 1.0) - Fa.jjp1() - Inuc * (Inuc + 1.0));
  }

private:
  double Inuc;

private:
  virtual double angularCff(int, int) const override { return 0; }
  virtual double angularCgg(int, int) const override { return 0; }
  virtual double angularCfg(int, int) const override { return 1.0; }
  virtual double angularCgf(int, int) const override { return 1.0; }
};

//******************************************************************************
class DirectHamiltonian : public ScalarOperator {
  // Direct part of Radial Hamiltonian operator.
public:
  DirectHamiltonian(const std::vector<double> &vn,
                    const std::vector<double> &vd, const double alpha)
      : ScalarOperator(OperatorParity::even, 1.0, {}), cl(1.0 / alpha),
        vnuc(vn), vdir(vd), tmp_f(std::vector<double>(vn.size(), 0.0)) {}

public:
  void updateVdir(const std::vector<double> &in_vdir) { vdir = in_vdir; }
  const std::vector<double> &getVdir() const { return vdir; }
  const std::vector<double> &getVnuc() const { return vnuc; }

public:
  double reducedME(const DiracSpinor &Fa,
                   const DiracSpinor &Fb) const override {
    auto inv_threej = std::sqrt(Fb.twoj() + 1.0);
    return inv_threej * matrixEl(Fa, Fb);
  }
  double matrixEl(const DiracSpinor &Fa, const DiracSpinor &Fb) const override {
    if (Fa.k != Fb.k)
      return 0.0;
    const auto &drdu = Fa.p_rgrid->drdu;
    tmp_f = NumCalc::derivative(Fa.f, drdu, Fa.p_rgrid->du, 1);
    const auto max = Fa.pinf;
    for (std::size_t i = 0; i < max; i++) {
      tmp_f[i] += (Fa.k * Fa.f[i] / Fa.p_rgrid->r[i]) - cl * Fa.g[i];
    }
    auto Hz = 2.0 * cl * NumCalc::integrate(Fa.g, tmp_f, drdu, 1.0, 0, max);
    auto Hw = NumCalc::integrate(Fa.f, Fa.f, vnuc, drdu, 1.0, 0, max) +
              NumCalc::integrate(Fa.g, Fa.g, vnuc, drdu, 1.0, 0, max) +
              NumCalc::integrate(Fa.f, Fa.f, vdir, drdu, 1.0, 0, max) +
              NumCalc::integrate(Fa.g, Fa.g, vdir, drdu, 1.0, 0, max);
    return (Hw + Hz) * Fa.p_rgrid->du;
  }

private:
  const double cl;
  const std::vector<double> vnuc;
  std::vector<double> vdir;
  mutable std::vector<double> tmp_f;
};

//******************************************************************************
class PNCnsiOperator : public ScalarOperator {
  // PNCnsi = (Gf*10^11/2sqrt[2]) * rho(r) * gamma5
  // Ouput is in units of (Qw * 1.e-11.)
  // To get (Qw/-N), multiply by (-N) [can go into optional 'factor']
public:
  PNCnsiOperator(double c, double t, const Grid &rgrid, double factor = 1)
      : ScalarOperator(OperatorParity::odd,
                       factor * PhysConst::GFe11 / std::sqrt(8.),
                       Nuclear::fermiNuclearDensity_tcN(t, c, 1, rgrid),
                       GammaMatrix::ig5) {}
};
