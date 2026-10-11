#include "Hypergeometric.hpp"
#include <cassert>
#include <cmath>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_result.h>
#include <limits>
#include <type_traits>
#include <utility>

#if defined(AMPSCI_USE_FLINT3) || defined(AMPSCI_USE_FLINT2)
// FLINT >= 3.0 merged Arb in: headers live under flint/ and link with -lflint
// (AMPSCI_USE_FLINT3). FLINT 2.x keeps acb_hypgeom in the separate Arb library:
// top-level headers, link with -lflint-arb (Debian/Ubuntu) or -larb (upstream)
// (AMPSCI_USE_FLINT2); see configure.sh. The acb_*/arb_*/arf_* API used below is
// identical either way.
#ifdef AMPSCI_USE_FLINT2
#include <acb.h>
#include <acb_hypgeom.h>
#else
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#endif
#endif

// Manuall switch off (for testing only) (doesn't change has_flint)
// #undef AMPSCI_USE_FLINT3
// #undef AMPSCI_USE_FLINT2

namespace Hypergeometric {

// Testing toggle (cpp-only; does not change has_flint). When true, use arb
// ball arithmetic and escalate the working precision until the result is good
// to double precision. When false, do a single fixed-precision call at
// double-ish precision ("call normally"). Escalation exists to handle the
// large cancellations for continuum (large imaginary) args.
[[maybe_unused]] constexpr bool use_arb = true;
// used only when use_arb == false
[[maybe_unused]] constexpr long fixed_prec = 64;

// Testing toggle (cpp-only). When true, complex H1f1 first attempts a fast
// double-precision evaluation (Maclaurin series or asymptotic expansion,
// each with a running error estimate), and only falls back to FLINT ball
// arithmetic when the estimated error is too large. This covers the vast
// majority of continuum-state evaluations at ~100x the speed of ball
// arithmetic. When false, always use FLINT (old behaviour).
constexpr bool use_fast_path = true;

// Estimated relative error above which the fast path defers to FLINT.
// The running estimates can be slightly optimistic (a few x), are relative
// to the modulus (a component can be several x less accurate), and the
// result should be good to ~1e-11: so sit well below that:
constexpr double fast_path_tol = 1.0e-13;

//------------------------------------------------------------------------------
// ln Gamma(z) for complex z, via GSL. The imaginary part (phase) is
// returned mod 2*pi, which is fine inside exp(). Returns NaN on error
// (e.g., at poles), which safely poisons the fast path (falls back).
std::complex<double> lngamma_complex(std::complex<double> z) {
  // static: set once only (thread safe)
  [[maybe_unused]] static const auto hndl = gsl_set_error_handler_off();
  gsl_sf_result lnr, arg;
  const auto status = gsl_sf_lngamma_complex_e(z.real(), z.imag(), &lnr, &arg);
  if (status != 0) {
    return {std::nan(""), std::nan("")};
  }
  return {lnr.val, arg.val};
}

//------------------------------------------------------------------------------
// Maclaurin series for e^{s} * 1F1(a, b, z): sum_n (a)_n z^n / ((b)_n n!).
// Returns {value, estimated relative error}; the estimate tracks float
// cancellation (eps * max|term| / |sum|). Non-convergence gives err = inf.
// The scale e^{s} is folded in up-front: the cancellation estimate is
// unchanged, but intermediate sums cannot overflow when e^{s} is small and
// compensating a large 1F1 (the usual continuum case).
std::pair<std::complex<double>, double>
series_1f1(std::complex<double> a, double b, std::complex<double> z, double s) {
  std::complex<double> term = std::exp(s);
  std::complex<double> sum = term;
  double max_term = std::abs(term);
  constexpr int max_n = 1500;
  bool converged = false;
  for (int n = 0; n < max_n; ++n) {
    term *= (a + double(n)) * z / ((b + double(n)) * double(n + 1));
    sum += term;
    max_term = std::max(max_term, std::abs(term));
    if (n > 3 && std::abs(term) < 1.0e-17 * std::abs(sum)) {
      converged = true;
      break;
    }
  }
  if (!converged) {
    return {sum, std::numeric_limits<double>::infinity()};
  }
  return {sum, 1.0e-16 * max_term / std::abs(sum)};
}

//------------------------------------------------------------------------------
// Divergent asymptotic series sum_n (p)_n (q)_n w^n / n!, truncated at the
// smallest term (optimal truncation). Returns {sum, estimated relative
// error}; the error is of order the smallest included term.
std::pair<std::complex<double>, double> asymptotic_sum(std::complex<double> p,
                                                       std::complex<double> q,
                                                       std::complex<double> w) {
  std::complex<double> term{1.0, 0.0};
  std::complex<double> sum = term;
  double err = 1.0;
  constexpr int max_n = 128;
  for (int n = 0; n < max_n; ++n) {
    const auto next =
      term * (p + double(n)) * (q + double(n)) * w / double(n + 1);
    if (n > 2 && std::abs(next) > std::abs(term)) {
      // terms started growing: stop at the smallest (previous) term
      break;
    }
    term = next;
    sum += term;
    err = std::abs(term);
    if (err < 1.0e-17 * std::abs(sum)) {
      break;
    }
  }
  return {sum, err / std::abs(sum)};
}

//------------------------------------------------------------------------------
// Large-|z| asymptotic expansion of e^{s} * 1F1(a, b, z) (DLMF 13.7.2):
//   1F1(a,b,z) ~ Gamma(b) * [ e^{sg*i*pi*a} z^{-a} / Gamma(b-a) * S1
//                           + e^{z} z^{a-b} / Gamma(a) * S2 ],
// with S1 = sum_n (a)_n (a-b+1)_n / (n! (-z)^n),
//      S2 = sum_n (b-a)_n (1-a)_n / (n! z^n),
// and sg = +1 for Im(z) >= 0, -1 otherwise. The scale s and all Gamma
// factors are combined inside the exponents, so there is no intermediate
// over/underflow. Returns {value, estimated relative error}.
std::pair<std::complex<double>, double> asymptotic_1f1(std::complex<double> a,
                                                       double b,
                                                       std::complex<double> z,
                                                       double s) {
  const std::complex<double> I{0.0, 1.0};
  const auto [S1, e1] = asymptotic_sum(a, a - b + 1.0, -1.0 / z);
  const auto [S2, e2] = asymptotic_sum(b - a, 1.0 - a, 1.0 / z);
  const auto lnGb = lngamma_complex({b, 0.0});
  const auto lnz = std::log(z);
  const double sg = z.imag() >= 0.0 ? 1.0 : -1.0;
  const auto t1 =
    std::exp(s + lnGb - lngamma_complex(b - a) + sg * I * M_PI * a - a * lnz) *
    S1;
  const auto t2 =
    std::exp(s + lnGb - lngamma_complex(a) + z + (a - b) * lnz) * S2;
  const auto sum = t1 + t2;
  const auto mag = std::abs(t1) + std::abs(t2);
  // truncation errors of each series, plus float cancellation between the
  // two terms:
  const auto err =
    (e1 * std::abs(t1) + e2 * std::abs(t2) + 1.0e-16 * mag) / std::abs(sum);
  return {sum, err};
}

//------------------------------------------------------------------------------
// Fast (double precision) evaluation of e^{s} * 1F1(a, b, z), returning
// {value, estimated relative error}. Chooses between the Maclaurin series
// and the large-|z| asymptotic expansion. If the returned estimate exceeds
// fast_path_tol, the value must not be used (fall back to FLINT).
std::pair<std::complex<double>, double>
fast_1f1(std::complex<double> a, double b, std::complex<double> z, double s) {
  const auto inf = std::numeric_limits<double>::infinity();
  // any non-finite value or estimate means "failed":
  const auto sanitise = [=](std::pair<std::complex<double>, double> v) {
    if (!std::isfinite(v.first.real()) || !std::isfinite(v.first.imag()) ||
        !std::isfinite(v.second)) {
      v.second = inf;
    }
    return v;
  };
  // Beyond this, the asymptotic expansion typically reaches full double
  // precision (it needs |z| somewhat larger than |a|^2 to converge); below
  // it, the series is usually both accurate and cheaper:
  const bool big_z = std::abs(z) > 25.0 + 0.1 * std::norm(a);
  // The series needs ~|z| terms, and its float cancellation grows
  // exponentially with |z|: hopeless beyond ~700 (overflows double):
  const bool series_viable = std::abs(z) < 700.0;

  auto best = std::pair<std::complex<double>, double>{{0.0, 0.0}, inf};
  if (big_z) {
    best = sanitise(asymptotic_1f1(a, b, z, s));
    if (best.second < fast_path_tol) {
      return best;
    }
  }
  if (series_viable) {
    const auto ser = sanitise(series_1f1(a, b, z, s));
    if (ser.second < best.second) {
      best = ser;
    }
    if (best.second < fast_path_tol) {
      return best;
    }
  }
  if (!big_z) {
    const auto asym = sanitise(asymptotic_1f1(a, b, z, s));
    if (asym.second < best.second) {
      best = asym;
    }
  }
  return best;
}

//==============================================================================

template <typename T>
T H1f1(T a, double b, T z, double s) {
  static_assert(std::is_same_v<T, double> ||
                  std::is_same_v<T, std::complex<double>>,
                "H1f1: T must be double or std::complex<double>");

  if constexpr (std::is_same_v<T, double>) {
    // Real case: GSL, in double precision
    // static: set once only (thread safe)
    [[maybe_unused]] static const auto hndl = gsl_set_error_handler_off();
    // GSL recurses without end (stack overflow) when both a and b are nan
    if (std::isnan(a) || std::isnan(b) || std::isnan(z)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    gsl_sf_result res;
    const auto status = gsl_sf_hyperg_1F1_e(a, b, z, &res);
    const auto val = status == 0 ? res.val : 0.0;
    return s == 0.0 ? val : std::exp(s) * val;
  } else {
    // Complex case. First, attempt the fast double-precision evaluation
    // (Maclaurin series or asymptotic expansion, with error estimates);
    // this covers the vast majority of continuum-state calls. Fall back to
    // FLINT where the estimated error is too large (strong cancellation:
    // large Im(a) with moderate |z|):
    if constexpr (use_fast_path) {
      const auto [val, err] = fast_1f1(a, b, z, s);
      if (err < fast_path_tol) {
        return val;
      }
    }
#if defined(AMPSCI_USE_FLINT3) || defined(AMPSCI_USE_FLINT2)
    // Complex case: FLINT ball arithmetic; working precision is increased
    // until the result is accurate to (at least) full double precision
    // https://flintlib.org/doc/acb_hypgeom.html
    acb_t res, aa, bb, zz;
    acb_init(res);
    acb_init(aa);
    acb_init(bb);
    acb_init(zz);
    acb_set_d_d(aa, a.real(), a.imag());
    acb_set_d(bb, b);
    acb_set_d_d(zz, z.real(), z.imag());
    // Following example ("arb"itrary precission):
    // https://arblib.org/using.html
    // Starting precision: the cancellation (in bits) scales with Im(a)
    // (fast path fails, and lands here, when Im(a) is large with moderate
    // |z|), so skip the doomed low-precision attempts:
    slong prec = 64;
    const double bits_est = 53.0 + 3.0 * std::abs(a.imag());
    while (double(prec) < bits_est && prec < 16384) {
      prec *= 2;
    }
    if constexpr (use_arb) {
      for (; prec <= 16384; prec *= 2) {
        acb_hypgeom_m(res, aa, bb, zz, 0, prec);
        if (acb_rel_accuracy_bits(res) >= 53) {
          break;
        }
      }
    } else {
      // Single fixed-precision call ("call normally"):
      prec = fixed_prec;
      acb_hypgeom_m(res, aa, bb, zz, 0, prec);
    }
    if (s != 0.0) {
      // apply scale e^{s} at matching precision (no accuracy loss)
      arb_t scale;
      arb_init(scale);
      arb_set_d(scale, s);
      arb_exp(scale, scale, prec);
      acb_mul_arb(res, res, scale, prec);
      arb_clear(scale);
    }
    const auto re = arf_get_d(arb_midref(acb_realref(res)), ARF_RND_NEAR);
    const auto im = arf_get_d(arb_midref(acb_imagref(res)), ARF_RND_NEAR);
    acb_clear(res);
    acb_clear(aa);
    acb_clear(bb);
    acb_clear(zz);
    return {re, im};
#else
    // Stub, when compiled without FLINT:
    // Run-time failure" if ever called:
    assert(false && "Requires FLINT to call H1f1 with complex args. See docs");
    (void)a;
    (void)b;
    (void)z;
    (void)s;
    return {std::nan(""), std::nan("")};
#endif
  }
}

// explicit instantiations:
template double H1f1(double, double, double, double);
template std::complex<double> H1f1(std::complex<double>, double,
                                   std::complex<double>, double);

} // namespace Hypergeometric
