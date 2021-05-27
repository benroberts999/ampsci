#include "RadiativePotential.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <iostream>
#include <string>
#include <vector>

static inline double ExpInt1(double x) {
  // E1(x) = -std::expint(-x)
  // return -std::expint(-x);
  return gsl_sf_expint_E1(x);
}

namespace RadiativePotential {
using namespace Helper;

std::vector<double> form_Hel(const std::vector<double> &r, double x_simple,
                             double x_Ueh, double x_SEe_h, double x_SEe_l,
                             double r_rms_Fermi, double z, double alpha,
                             double rcut) {

  const auto rN_rad = r_rms_Fermi * std::sqrt(5.0 / 3.0) / PhysConst::aB_fm;

  if (rcut == 0.0)
    rcut = r.back();
  const auto imax = r.size();
  std::vector<double> Vel_tmp(r.size());
  std::cout << "\n";

  if (std::abs(x_simple) > 0) {
    std::cout << "Forming simple exponential potential (Ak=" << x_simple
              << ")\n";
    for (std::size_t i = 0; i < imax; ++i) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      Vel_tmp[i] -=
          x_simple * RadiativePotential::vSimpleExp(ri, rN_rad, z, alpha);
    }
  }

  if (x_Ueh > 0) {
    std::cout << "Forming Uehling potential (scale=" << x_Ueh << ")\n";
    // #pragma omp parallel for
    for (std::size_t i = 0; i < imax; ++i) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      const auto v_Ueh = RadiativePotential::vUehling(ri, rN_rad, z, alpha);
      Vel_tmp[i] -= x_Ueh * v_Ueh;
    }
  }

  if (x_SEe_h > 0) {
    std::cout << "Forming Self-Energy (electric high-f) potential (scale="
              << x_SEe_h << ")\n";
    // The SE high-f part is very slow. We calculate it every few points only,
    // and then interpolate the intermediate points
    // const std::size_t stride = 6; // r.size()
    const auto stride = std::max(r.size() / 1000, 1ul);
    const auto i_max_rcut = [&]() {
      for (std::size_t i = 0; i < imax; ++i) {
        if (r[i] > rcut)
          return i;
      }
      return imax;
    }();
    const auto tmp_max = i_max_rcut / stride;
    std::vector<double> x(tmp_max), y(tmp_max);
#pragma omp parallel for
    for (std::size_t i = 0; i < tmp_max; ++i) {
      x[i] = r[i * stride];
      y[i] = RadiativePotential::vSEh(x[i], rN_rad, z, alpha);
    }
    const auto vec_SE_h = Interpolator::interpolate(x, y, r);
    for (std::size_t i = 0; i < imax; i++) {
      if (r[i] > rcut)
        break;
      Vel_tmp[i] -= x_SEe_h * vec_SE_h[i];
    }
  }

  if (x_SEe_l > 0) {
    std::cout << "Forming Self-Energy (electric low-f) potential (scale="
              << x_SEe_l << ")\n";
    // #pragma omp parallel for
    for (std::size_t i = 0; i < imax; i++) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      const auto v_SE_l = RadiativePotential::vSEl(ri, rN_rad, z, alpha);
      Vel_tmp[i] -= x_SEe_l * v_SE_l;
    }
  }

  return Vel_tmp;
}
//------------------------------------------------------------------------------
std::vector<double> form_Hmag(const std::vector<double> &r, double x_SEm,
                              double r_rms_Fermi, double z, double alpha,
                              double rcut) {

  const auto rN_rad = r_rms_Fermi * std::sqrt(5.0 / 3.0) / PhysConst::aB_fm;

  if (rcut == 0)
    rcut = r.back();
  const auto imax = r.size();
  std::vector<double> Hmag_tmp(r.size());

  if (x_SEm > 0) {
    std::cout << "Forming Self-Energy (magnetic) potential "
              << "(scale=" << x_SEm << ")\n";
    // #pragma omp parallel for
    for (std::size_t i = 0; i < imax; i++) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      const auto Hr = RadiativePotential::vSE_Hmag(ri, rN_rad, z, alpha);
      Hmag_tmp[i] += x_SEm * Hr; // nb: plus!
    }
  }
  return Hmag_tmp;
}

//******************************************************************************
bool read_write_qed(const std::vector<double> &r, std::vector<double> &Hel,
                    std::vector<double> &Hmag, const std::string &fname,
                    IO::FRW::RoW rw) {
  //
  const auto readQ = rw == IO::FRW::read;

  if (readQ && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str = !readQ ? "Writing to " : "Reading from ";
  std::cout << rw_str << "QED rad. pot. file: " << fname << " ... "
            << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  // read-write grid:
  std::vector<double> r_in;
  auto sizer = r.size();
  rw_binary(iofs, rw, sizer);
  if (readQ) {
    r_in.resize(sizer);
    for (auto &ri : r_in)
      rw_binary(iofs, rw, ri);
  } else {
    for (auto ri : r)
      rw_binary(iofs, rw, ri);
  }

  // check if grid the same
  bool grid_ok = r.size() == r_in.size();
  if (grid_ok) {
    for (auto i = 0ul; i < r.size(); i += 10) {
      if (std::abs(r[i] / r_in[i] - 1.0) > 0.001) {
        grid_ok = false;
        break;
      }
    }
  }

  auto sizeE = Hel.size();
  rw_binary(iofs, rw, sizeE);
  if (readQ) {
    Hel.resize(sizeE);
  }
  for (auto &hel : Hel) {
    rw_binary(iofs, rw, hel);
  }

  auto sizeM = Hmag.size();
  rw_binary(iofs, rw, sizeM);
  if (readQ) {
    Hmag.resize(sizeM);
  }
  for (auto &hmag : Hmag) {
    rw_binary(iofs, rw, hmag);
  }

  if (readQ && !grid_ok) {
    std::cout << "\nInterpolating QED rad-pot onto current grid ... ";
    auto tmp_copy_Hrad = Hel;
    Hel = Interpolator::interpolate(r_in, tmp_copy_Hrad, r);
    tmp_copy_Hrad = Hmag;
    Hmag = Interpolator::interpolate(r_in, tmp_copy_Hrad, r);
  }

  std::cout << "done.\n";
  return true;
}

//******************************************************************************
//******************************************************************************
double vSimpleExp(double r, double, double z, double alpha) {
  return -z * z * alpha * std::exp(-r / alpha);
}

//******************************************************************************
double vUehling(double r, double rN, double z, double alpha) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);

  // Routines return the first approximation which has an absolute error
  // smaller than abs_err_lim or a relative error smaller than rel_err_lim.
  // Note that this is an either-or constraint, not simultaneous. To compute to
  // a specified absolute error, set epsrel to zero (etc.)
  static constexpr double abs_err_lim = 0.0;
  static constexpr double rel_err_lim = 1.0e-6;
  // max_num_subintvls < size(gsl_int_wrk)
  static constexpr unsigned long max_num_subintvls = 750; //?

  gsl_set_error_handler_off(); //?
  if (rN <= 0) {
    rN = 1.0e-7;
  }

  RadPot_params params = {r, rN, z, alpha};

  gsl_function f_gsl;
  f_gsl.function = (r < rN) ? &gslfunc_Ueh_smallr : &gslfunc_Ueh_larger;
  f_gsl.params = &params;

  // This workspace handles the memory for the subinterval ranges, results and
  // error estimates. max_num_subintvls < size(gsl_int_wrk)
  // Allocates a workspace sufficient to hold n double precision
  // intervals, their integration results and error estimates.
  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);
  // nb: i allocate + destroy wrk EACH r... doesn't much matter though...

  double int_result, abs_err;
  gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                        max_num_subintvls, gsl_int_wrk, &int_result, &abs_err);
  gsl_integration_workspace_free(gsl_int_wrk);

  auto pre_factor = z * (alpha / M_PI) / r;

  return pre_factor * int_result;
}

//******************************************************************************
double Q_MLVP(double r, double rN) {
  // magnetic loop vacuum polarisation (VP vertex)
  // Performs t integral. This multiplies regular operator

  // variables needed for gsl
  constexpr double abs_err_lim = 0.0;
  constexpr double rel_err_lim = 1.0e-6;
  constexpr unsigned long max_num_subintvls = 750;
  gsl_set_error_handler_off();

  // compute the integral at each radial grid point

  // intialise gsl
  gsl_function f_gsl;
  RadiativePotential::MLVP_params params{r, rN};
  f_gsl.function = &Helper::Qt_MLVP;
  f_gsl.params = &params;

  // integrate using gsl
  double result{0.0};
  double abs_err{0.0};
  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);
  gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                        max_num_subintvls, gsl_int_wrk, &result, &abs_err);
  gsl_integration_workspace_free(gsl_int_wrk);

  return result;
}

//------------------------------------------------------------------------------
double Helper::Qt_MLVP(double t, void *p) {
  // This in Qt, defined via Q(r) = Int[ Qt(t,r) , {t,1,infty}]

  // Main:
  // A. V. Volotka, D. A. Glazov, I. I. Tupitsyn, N. S. Oreshkina, G.
  // Plunien, and V. M. Shabaev, Phys. Rev. A 78, 062507 (2008) - Equation (18)

  // also:
  // P. Sunnergren, H. Persson, S. Salomonson, S. M. Schneider, I. Lindgren,
  // and G. Soff, Phys. Rev. A 58, 1055 (1998) -- Equation (58)
  // * note: This has typo compared to above

  // obtain the radial grid parameters and convert to relativistic units
  const auto [r, r_N] = *(static_cast<MLVP_params *>(p));

  const auto chi = std::min(r, r_N) / PhysConst::alpha;
  const auto eta = std::max(r, r_N) / PhysConst::alpha;

  constexpr double a0 = PhysConst::alpha / (3.0 * M_PI);
  const auto t2 = t * t;
  const auto top = std::sqrt(t2 - 1.0) * (2 * t2 + 1.0) * (2 * t * eta + 1.0);
  const auto exp_t4 = std::exp(-2 * t * eta) / (t2 * t2);
  const auto rn_tmp = r_N; // to avoid clang error re: capturing binding ref
  const auto finite = [=]() {
    if (rn_tmp > 1.0e-6) {
      const auto x = 2.0 * chi * t;
      return (3.0 / (x * x * x)) * xCoshxMinusSinhx(x);
    }
    return 1.0;
  }(); // immediately invoked lambda

  return a0 * top * exp_t4 * finite;
}

//******************************************************************************
static Fit_AB fit_AB;
double vSEh(double r, double rN, double z, double alpha) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);

  static constexpr double abs_err_lim = 1.0e-6;
  static constexpr double rel_err_lim = 1.0e-3;
  static constexpr unsigned long max_num_subintvls = 750; //?

  gsl_set_error_handler_off(); //?
  if (rN <= 0) {
    rN = 1.0e-7;
  }

  RadPot_params params = {r, rN, z, alpha};

  gsl_function f_gsl;
  f_gsl.function = (r < rN) ? &gslfunc_SEh_smallr : &gslfunc_SEh_larger;
  f_gsl.params = &params;

  double int_result, abs_err;
  auto rel_err_targ = rel_err_lim;
  int i = 0;
  while (rel_err_targ < 1.0) {
    gsl_integration_workspace *gsl_int_wrk =
        gsl_integration_workspace_alloc(max_num_subintvls + 1);
    gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_targ,
                          max_num_subintvls, gsl_int_wrk, &int_result,
                          &abs_err);
    gsl_integration_workspace_free(gsl_int_wrk);
    ++i;
    if (std::abs(abs_err / int_result) < rel_err_targ ||
        abs_err <= abs_err_lim) {
      break;
    } else {
      int_result = 0.0;
      abs_err = 0.0;
      rel_err_targ *= 10.0;
    }
  }

  auto l = 0; // XXX
  auto al = fit_AB.Al(l, z * alpha);

  auto pre_factor = -1.5 * al * (alpha / M_PI) * z / r;

  return pre_factor * int_result;
}

//------------------------------------------------------------------------------
double vSEl(double r, double rN, double z, double alpha) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  //
  auto l = 0; // XXX
  auto bl = fit_AB.Bl(l, z * alpha);

  if (rN <= 0) {
    rN = 1.0e-7;
  }

  auto pre_factor =
      -1.5 * bl * z * z * alpha * alpha * alpha / (rN * rN * rN * r);

  auto f = [=](double x) {
    auto arg_neg = z * std::abs(r - x);
    auto arg_pos = z * (r + x);
    auto a = (arg_neg + 1.0) * std::exp(-arg_neg);
    auto b = (arg_pos + 1.0) * std::exp(-arg_pos);
    return x * (a - b);
  };

  return pre_factor * NumCalc::num_integrate(f, 0.0, rN, 1000, NumCalc::linear);

  // This is the "point-nucleus" version. Gives same result!
  // auto pre_factor2 = -bl * z * z * z * z * alpha * alpha * alpha;
  // return pre_factor2 * std::exp(-z * r);
}

//******************************************************************************
double vSE_Hmag(double r, double rN, double z, double alpha) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);

  static constexpr double abs_err_lim = 0.0;
  static constexpr double rel_err_lim = 1.0e-6;
  static constexpr unsigned long max_num_subintvls = 750; //?

  gsl_set_error_handler_off(); //?
  if (rN <= 0) {
    rN = 1.0e-7;
  }

  RadPot_params params = {r, rN, z, alpha};

  gsl_function f_gsl;
  f_gsl.function = &gslfunc_SEmag;
  f_gsl.params = &params;

  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);

  double int_result, abs_err;
  gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                        max_num_subintvls, gsl_int_wrk, &int_result, &abs_err);
  gsl_integration_workspace_free(gsl_int_wrk);

  const auto pre_factor = 3.0 * z * alpha / M_PI;
  return pre_factor * int_result;
}

//******************************************************************************
//******************************************************************************

namespace Helper {

//------------------------------------------------------------------------------
double vUehcommon(double t, double chi) {
  auto t2 = t * t;
  auto chi3 = chi * chi * chi;
  auto a = std::sqrt(t2 - 1.0);
  auto b = (1.0 / t2 + 1.0 / (2.0 * t2 * t2));
  auto c = 2.0 / chi3;
  return a * b * c;
}

//------------------------------------------------------------------------------
double vUehf_smallr(double r, double rN, double chi) {
  const auto a = (r / rN) * chi;
  const auto b = std::exp(-chi) * (1.0 + chi);
  const auto c = std::sinh(a);
  return a - b * c;
}
double vUehf_larger(double r, double rN, double chi) {
  const auto tmp = (r / rN) * chi;
  const auto a = std::exp(-tmp);
  const auto b = xCoshxMinusSinhx(chi);
  return a * b;
}

//------------------------------------------------------------------------------
double gslfunc_Ueh_smallr(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto g = vUehcommon(t, chi);
  const auto f = vUehf_smallr(r, rN, chi);
  return f * g;
}
double gslfunc_Ueh_larger(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto g = vUehcommon(t, chi);
  const auto f = vUehf_larger(r, rN, chi);
  return f * g;
}

//******************************************************************************
double gb_I1(double t, double z, double alpha) {
  const auto t2 = t * t;
  const auto za = z * alpha;
  const double a = 1.0 / std::sqrt(t2 - 1.0);
  const double b = 1.0 - 1.0 / (2.0 * t2);
  const double c = std::log(t2 - 1.0) + 4.0 * std::log(1.0 / za + 0.5);
  const double d = -1.5 + 1.0 / t2;
  return a * (b * c + d);
}

double gb_I2(double t, double r, double rN, double z, double alpha) {
  const auto za = z * alpha;
  const auto ttoa = 2.0 * t / alpha;
  const double rA = 0.07 * za * za * alpha;
  const auto expr = std::exp(rA * ttoa);

  const auto f = [=](double x) {
    const auto arg1 = (std::abs(r - x) + rA) * ttoa;
    const auto arg2 = (r + x + rA) * ttoa;
    const auto a = ExpInt1(arg1) - ExpInt1(arg2);
    return x * a;
  };

  return expr * NumCalc::num_integrate(f, 0.0, rN, 1000, NumCalc::linear);
}

//------------------------------------------------------------------------------
double gb_GSEh_larger(double r, double rN, double chi) {
  // Essentially duplicate of Uehling, but swapped large/small r ?
  const auto chi3 = chi * chi * chi;
  const auto a = std::exp(-chi * r / rN) * 2.0 / chi3;
  const auto b = xCoshxMinusSinhx(chi);
  return a * b;
}

double gb_GSEh_smallr(double r, double rN, double chi) {
  // Essentially duplicate of Uehling, but swapped large/small r ?
  const auto chi3 = chi * chi * chi;
  const auto a = 2.0 / chi3;
  const auto b = (chi * r / rN);
  const auto c = std::exp(-chi) * (1.0 + chi) * std::sinh(chi * r / rN);
  return a * (b - c);
}

//------------------------------------------------------------------------------
double gslfunc_SEh_smallr(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  const auto chi = 2.0 * t * rN / alpha;
  const double rA = 0.07 * z * z * alpha * alpha * alpha;
  return gb_I1(t, z, alpha) *
         (gb_GSEh_smallr(r, rN, chi) -
          (rA / (rN * rN * rN)) * gb_I2(t, r, rN, z, alpha));
}
double gslfunc_SEh_larger(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  const auto chi = 2.0 * t * rN / alpha;
  const double rA = 0.07 * z * z * alpha * alpha * alpha;
  return gb_I1(t, z, alpha) *
         (gb_GSEh_larger(r, rN, chi) -
          (rA / (rN * rN * rN)) * gb_I2(t, r, rN, z, alpha));
}

//******************************************************************************
double seJmag(double x, double y) {
  const auto y3 = y * y * y;
  const auto sihncosh =
      (y < 0.5) ? -y3 / 3.0 - (y3 * y * y) / 30.0 - (y3 * y3 * y) / 840.0
                : sinh(y) - y * std::cosh(y);
  return std::exp(-x) * (1.0 + x) * sihncosh + y3 / 3.0;
}

double gslfunc_SEmag(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto eta = 2.0 * t * r / alpha;
  const auto factor =
      1.0 / (chi * chi * chi * eta * eta * std::sqrt(t * t - 1.0));
  return (r <= rN) ? factor * seJmag(chi, eta) : factor * seJmag(eta, chi);
}

} // namespace Helper
} // namespace RadiativePotential
