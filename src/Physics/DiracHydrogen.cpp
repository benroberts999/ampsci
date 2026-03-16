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

double Norm(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  const auto denom = nn(n, k, z, a) * Gamma(2.0 * gamma(k, z, a) + 1);
  const auto argGamma = n.v + double(1 - std::abs(k.v));
  const auto arg1 = z.v * Gamma(2.0 * gamma(k, z, a) + argGamma);
  const auto arg2 = 2.0 * Gamma(argGamma) * (nn(n, k, z, a) - double(k.v));
  return std::sqrt(arg1 / arg2) / denom;
}

double lambda(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  // const auto a2 = (a * a).v;
  // const auto en = Enk(n, k, z, a);
  // return std::sqrt(1.0 / a2 - a2 * en * en);
  const auto e = enk(n, k, z, a);
  const auto c2 = 1.0 / (a * a).v;
  return std::sqrt(-e * (2.0 + e / c2));
}

double x(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  return r.v * 2.0 * lambda(n, k, z, a);
}
} // namespace Hidden

//==============================================================================
double enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  using namespace Hidden;
  const auto a2 = (a * a).v;
  const auto c2 = 1.0 / a2;
  const auto arg = gamma(k, z, a) + n.v - double(std::abs(k.v));
  const auto w2 = (z * z).v / (arg * arg);
  const auto d = 1.0 + a2 * w2;
  return -w2 / (2.0 * d) -
         (0.5 * a2 * w2 + 1.0 - std::sqrt(1.0 + a2 * w2)) * (c2 / d);
}

double Enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  const auto a2 = (a * a).v;
  return 1.0 / a2 + enk(n, k, z, a);
}

double gamma(DiracQN k, Zeff z, AlphaFS a) {
  return std::sqrt(double((k * k).v) - (a * a).v * (z * z).v); //
}

double f(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  using namespace Hidden;
  const auto xr = x(r, n, k, z, a);
  const auto g = gamma(k, z, a);
  const auto kmn = double(std::abs(k.v)) - n.v;
  const auto c1 = std::sqrt(1.0 + (a * a).v * Enk(n, k, z, a));
  const auto c2 = Norm(n, k, z, a) * std::exp(-0.5 * xr) * std::pow(xr, g);
  const auto d1 = (nn(n, k, z, a) - double(k.v)) * H1f1(kmn, 2.0 * g + 1.0, xr);
  const auto d2 = kmn * H1f1(kmn + 1.0, 2.0 * g + 1.0, xr);
  return c1 * c2 * (d1 + d2);
}

double g(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a) {
  using namespace Hidden;
  const auto xr = x(r, n, k, z, a);
  const auto g = gamma(k, z, a);
  const auto kmn = double(std::abs(k.v)) - n.v;
  const auto c1 = std::sqrt(1.0 - (a * a).v * Enk(n, k, z, a));
  const auto c2 = Norm(n, k, z, a) * std::exp(-0.5 * xr) * std::pow(xr, g);
  const auto d1 = (nn(n, k, z, a) - double(k.v)) * H1f1(kmn, 2.0 * g + 1.0, xr);
  const auto d2 = kmn * H1f1(kmn + 1.0, 2.0 * g + 1.0, xr);
  return -c1 * c2 * (d1 - d2);
}

} // namespace DiracHydrogen
