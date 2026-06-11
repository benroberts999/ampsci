#include "DiracHydrogen.hpp"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_result.h>

namespace DiracHydrogen {

//==============================================================================
namespace Hidden {
// Some helper functions:

double H1f1(double a, double b, double x) { //
  // return gsl_sf_hyperg_1F1(a, b, x);
  gsl_sf_result gsl_res;
  gsl_set_error_handler_off();
  auto status = gsl_sf_hyperg_1F1_e(a, b, x, &gsl_res);
  // return gsl_res.val;
  return status == 0 ? gsl_res.val : 0.0; //?
}
double Gamma(double x) { return std::tgamma(x); }

double nn(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  const auto ak = double(std::abs(k.v));
  return std::sqrt((n * n).v - 2.0 * (n.v - ak) * (ak - gamma(k, z, a)));
}

// N_{nk}: normalization prefactor [Methods Eq. (233)].
// Mass-independent: c1 ~ sqrt(m) and x = 2*lambda_m*r already account for the
// full m-scaling, f_m(r) = sqrt(m) f_1(m*r).
double Norm(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  const auto denom = nn(n, k, z, a) * Gamma(2.0 * gamma(k, z, a) + 1);
  const auto argGamma = n.v + double(1 - std::abs(k.v));
  const auto arg1 = z.v * Gamma(2.0 * gamma(k, z, a) + argGamma);
  const auto arg2 = 2.0 * Gamma(argGamma) * (nn(n, k, z, a) - double(k.v));
  return std::sqrt(arg1 / arg2) / denom;
}

double lambda(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m) {
  const auto e = enk(n, k, z, a, m);
  const auto a2 = (a * a).v;
  return std::sqrt(-e * (2.0 * m + a2 * e));
}

double x(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m) {
  return r.v * 2.0 * lambda(n, k, z, a, m);
}
} // namespace Hidden

//==============================================================================
// Stable form for enk (avoids cancellation between O(c^2) terms):
//   nbar = gamma + n - |k|
//   s    = sqrt(nbar^2 + (alpha*Z)^2)
//   enk  = -m*Z^2 / (s * (nbar + s))
double enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m) {
  const auto nbar = gamma(k, z, a) + n.v - double(std::abs(k.v));
  const auto az = a.v * z.v;
  const auto s = std::sqrt(nbar * nbar + az * az);
  return -m * z.v * z.v / (s * (nbar + s));
}

double Enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m) {
  return m / (a * a).v + enk(n, k, z, a, m);
}

double gamma(DiracQN k, Zeff z, AlphaFS a) {
  return std::sqrt(double((k * k).v) - (a * a).v * (z * z).v); //
}

double f(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m) {
  using namespace Hidden;
  const auto xr = x(r, n, k, z, a, m);
  const auto gam = gamma(k, z, a);
  const auto kmn = double(std::abs(k.v)) - n.v;
  const auto a2 = (a * a).v;
  const auto en = enk(n, k, z, a, m);
  const auto c1 = std::sqrt(2.0 * m + a2 * en);
  const auto c2 = Norm(n, k, z, a) * std::exp(-0.5 * xr) * std::pow(xr, gam);
  const auto d1 =
    (nn(n, k, z, a) - double(k.v)) * H1f1(kmn, 2.0 * gam + 1.0, xr);
  const auto d2 = kmn * H1f1(kmn + 1.0, 2.0 * gam + 1.0, xr);
  const auto sk = k.v < 0 ? 1.0 : -1.0;
  return sk * c1 * c2 * (d1 + d2);
}

double g(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m) {
  using namespace Hidden;
  const auto xr = x(r, n, k, z, a, m);
  const auto gam = gamma(k, z, a);
  const auto kmn = double(std::abs(k.v)) - n.v;
  const auto a2 = (a * a).v;
  const auto en = enk(n, k, z, a, m);
  const auto c1 = std::sqrt(-a2 * en);
  const auto c2 = Norm(n, k, z, a) * std::exp(-0.5 * xr) * std::pow(xr, gam);
  const auto d1 =
    (nn(n, k, z, a) - double(k.v)) * H1f1(kmn, 2.0 * gam + 1.0, xr);
  const auto d2 = kmn * H1f1(kmn + 1.0, 2.0 * gam + 1.0, xr);
  const auto sk = k.v < 0 ? 1.0 : -1.0;
  return -sk * c1 * c2 * (d1 - d2);
}

double gfratio(double r, int k, double z, double a, double e, double m) {
  using namespace Hidden;

  const auto gam = std::sqrt(double(k * k) - (a * a * z * z));
  const auto a2 = a * a;
  const auto absk = (double)std::abs(k);

  const auto n =
    ((z * (m + e * a2)) / std::sqrt(-(e * (2.0 * m + e * a2)))) - gam + absk;

  const auto xr = r * 2.0 * std::sqrt(-e * (2.0 * m + e * a2));
  const auto kmn = double(absk) - n;
  const auto c1_f = std::sqrt(2.0 * m + a2 * e);
  const auto nn = std::sqrt((n * n) - 2.0 * (n - absk) * (absk - gam));
  const auto d1 = (nn - double(k)) * H1f1(kmn, 2.0 * gam + 1.0, xr);
  const auto d2 = kmn * H1f1(kmn + 1.0, 2.0 * gam + 1.0, xr);
  const auto ff = c1_f * (d1 + d2);

  const auto c1_g = std::sqrt(-a2 * e);
  // sign factor s_kappa cancels in the ratio g/f
  const auto gg = c1_g * (d2 - d1);
  return gg / ff;
}

} // namespace DiracHydrogen
