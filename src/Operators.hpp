#pragma once
#include "DiracOperator.hpp"
#include "Grid.hpp"
#include "Physics/Nuclear.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>
#include <vector>

// Classes (inherit from DriacOperator)

// Various contructors: (all after v are optional)
// DiracOperator(C, v, DiracMatrix(a, b, c, d), diff_order, imag?)
// DiracOperator(C, DiracMatrix(a, b, c, d), diff_order, imag?)
// DiracOperator(v, DiracMatrix(a, b, c, d), diff_order, imag?)
// DiracOperator(DiracMatrix(a, b, c, d), diff_order, imag?)
// DiracOperator(diff_order, imag?)
// DiracOperator(DiracMatrix(a, b, c, d), imag?)

//******************************************************************************
static inline double hfs_SphericalBall_F(double r, double rN) {
  return (r > rN) ? 1. / (r * r) : r / (rN * rN * rN);
}

static inline std::vector<double> hfs_v(double rN, const Grid &rgrid,
                                        double (*hfs_F)(double, double)) {
  std::vector<double> invr2;
  invr2.reserve(rgrid.ngp);
  for (auto r : rgrid.r) {
    invr2.push_back(hfs_F(r, rN));
  }
  return invr2;
}

class HyperfineOperator : public DiracOperator {
public:
  HyperfineOperator(double muN, double IN, double rN, const Grid &rgrid,
                    double (*hfs_F)(double, double) = &hfs_SphericalBall_F)
      : DiracOperator(-0.5 * (muN / IN) * PhysConst::alpha / PhysConst::m_p,
                      hfs_v(rN, rgrid, hfs_F), DiracMatrix(0, 1, -1, 0), 0,
                      true) {}
};

//******************************************************************************
class PrOperator : public DiracOperator {
  // Pr = -i (d/dr)
public:
  PrOperator() : DiracOperator(-1, GammaMatrix::ident, 1, true) {}
};

//******************************************************************************
class DerivativeOperator : public DiracOperator {
  // = (d^n/dr^n), for any n>=0
public:
  DerivativeOperator(int d_order)
      : DiracOperator(GammaMatrix::ident, d_order) {}
};

//******************************************************************************
class E1Operator : public DiracOperator {
  // d_E1 = -er = -|e|r = -r
public:
  E1Operator(const Grid &rgrid) : DiracOperator(-1, rgrid.r) {}
};

//******************************************************************************
class PNCnsiOperator : public DiracOperator {
  // PNCnsi = (Gf*10^11/2sqrt[2]) * rho(r) * gamma5
  // Ouput is in units of (Qw * 1.e-11.)
  // To get (Qw/-N), multiply by (-N) [can go into optional 'factor']
public:
  PNCnsiOperator(double c, double t, const Grid &rgrid, double factor = 1)
      : DiracOperator(factor * PhysConst::GFe11 / sqrt(8.),
                      Nuclear::fermiNuclearDensity_tcN(t, c, 1, rgrid),
                      GammaMatrix::g5) {}
};

// double t = 2.3;
// double c = Nuclear::approximate_c_hdr(wf.Anuc());
// auto rho = Nuclear::fermiNuclearDensity_tcN(t, c, 1, wf.rgrid);
// double Gf = PhysConst::GFe11;
// double Cc = (Gf / sqrt(8.)) * (-wf.Nnuc()); // Qw/(-N)
// DiracOperator hpnc(Cc, rho, GammaMatrix::g5);

//******************************************************************************
class RadialOperator : public DiracOperator {
  // = some function of r
  // Pass only grid to just get r, or
  // either pass a lambda/function [f(r)], or a number, n, (for r^n)
public:
  RadialOperator(const Grid &rgrid, double (*f)(double r))
      : DiracOperator([&]() {
          std::vector<double> f_r;
          f_r.reserve(rgrid.ngp);
          for (auto r : rgrid.r) {
            f_r.push_back(f(r));
          }
          return f_r;
        }()) {}

  RadialOperator(const Grid &rgrid)
      : DiracOperator([&]() {
          std::vector<double> f_r;
          f_r.reserve(rgrid.ngp);
          for (auto r : rgrid.r) {
            f_r.push_back(r);
          }
          return f_r;
        }()) {}

  RadialOperator(const Grid &rgrid, int n)
      : DiracOperator([&]() {
          std::vector<double> f_r;
          f_r.reserve(rgrid.ngp);
          for (auto r : rgrid.r) {
            f_r.push_back(pow(r, n));
          }
          return f_r;
        }()) {}

  RadialOperator(const Grid &rgrid, double x)
      : DiracOperator([&]() {
          std::vector<double> f_r;
          f_r.reserve(rgrid.ngp);
          for (auto r : rgrid.r) {
            f_r.push_back(pow(r, x));
          }
          return f_r;
        }()) {}
};
