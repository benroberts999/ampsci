#pragma once
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/Operators/Vee.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <gsl/gsl_sf.h>

namespace BSM_Vee {

DiracSpinor V_Fv(const std::vector<DiracSpinor> &core, const bool eN,
                 const DiracSpinor &Fv, const std::string type,
                 const int kappa_n, const double y, const bool contact,
                 const double mu) {
  if ((type != "sp") && (type != "va")) {
    return 0.0 * Fv;
  } else if (type == "va" && contact == false) {
    return 0.0 * Fv;
  }

  DiracSpinor VFv(Fv.n(), Fv.kappa(), Fv.grid_sptr());
  const double relative_kappa = ((type == "ss") || (type == "vv")) ? 1.0 : -1.0;

  if (kappa_n != relative_kappa * Fv.kappa()) {
    return 0.0 * Fv;
  }

  // For eN (testing sign)
  if (eN) {

    // Make sure this sign is same as what we end up doing for ee
    const auto VFv_radial = g0(g5(Fv));
    double VFv_angular = 1.0;

    // Only for Caesium 133!
    const auto c = Nuclear::c_hdr_formula_rrms_t(Nuclear::find_rrms(55, 133));
    const auto t = Nuclear::default_t;
    const auto rho = Nuclear::fermiNuclearDensity_tcN(t, c, 1.0, Fv.grid());

    return rho * VFv_radial * VFv_angular;
  }

  // For ee
  for (auto Fa : core) {
    const auto Fdir = calc_Fdir(type, Fa, Fv, contact, mu);

    if (type == "va") {
      VFv += 2 * Fdir;
    } else {
      const auto Fexch = calc_Fexch(type, kappa_n, Fa, Fv, contact, mu);
      VFv += Fdir - Fexch;
    }
  }

  double mu_factor;

  if (type == "va" || contact == true) {
    mu_factor = 1.0 / (mu * mu);
  } else if (mu == 0) {
    mu_factor = 1.0;
  } else {
    mu_factor = mu;
  }

  return mu_factor * VFv;
}

DiracSpinor calc_Fdir(const std::string type, const DiracSpinor &Fa,
                      const DiracSpinor &Fv, const bool contact,
                      const double mu) {
  return (1.0 / (4 * M_PI)) * Fa.twojp1() *
         Bk_ab_v(0, contact, mu, true, type, Fa, Fa, Fv);
}

DiracSpinor calc_Fexch(const std::string type, const int kappa_n,
                       const DiracSpinor &Fa, const DiracSpinor &Fv,
                       const bool contact, const double mu) {

  DiracSpinor Fexch(Fa.n(), Fa.kappa(), Fa.grid_sptr());

  const auto phase =
      Angular::neg1pow_2(Fv.twoj() - Fa.twoj()) * (1.0 / (Fv.twojp1()));

  for (int k = 0; k <= 0.5 * (Fa.twoj() + Fv.twoj()); ++k) {

    const auto twok = 2.0 * k;
    const auto twokp1 = twok + 1;

    const auto A_naav = Angular::Ck_kk(k, kappa_n, Fa.kappa()) *
                        Angular::Ck_kk(k, Fa.kappa(), -Fv.kappa());
    const auto Fexch_AB =
        A_naav * Bk_ab_v(k, contact, mu, false, type, Fa, Fv, Fa);

    const auto A_anva = Angular::Ck_kk(k, Fa.kappa(), Fv.kappa()) *
                        Angular::Ck_kk(k, kappa_n, -Fa.kappa());
    const auto Fexch_BA =
        A_anva * Bk_ab_v(k, contact, mu, true, type, Fa, Fv, Fa);

    Fexch += phase * twokp1 * (Fexch_AB + Fexch_BA);
  }

  return (1.0 / (4 * M_PI)) * Fexch;
}

DiracSpinor Bk_ab_v(const int k, const bool contact, const double mu,
                    const bool betaalpha, const std::string type,
                    const DiracSpinor &Fa, const DiracSpinor &Fb,
                    const DiracSpinor &Fv) {
  const auto &gr = Fa.grid();
  const auto &r = gr.r();

  auto mod_Fb = Fb;
  auto mod_Fv = Fv;

  if (betaalpha) {
    if (type == "sp") {
      mod_Fb = g0(Fb);
      mod_Fv = -1.0 * g5(Fv);
    } else if (type == "va") {
      mod_Fv = g0(g5(Fv));
    } else if (type == "ss") {
      mod_Fb = g0(Fb);
      mod_Fv = g0(Fv);
    }
  } else {
    if (type == "sp") {
      mod_Fb = -1.0 * g5(Fb);
      mod_Fv = g0(Fv);
    } else if (type == "va") {
      mod_Fb = g0(g5(Fb));
    } else if (type == "ss") {
      mod_Fb = g0(Fb);
      mod_Fv = g0(Fv);
    }
  }

  if (contact == true) {
    std::vector<double> FaFb_rr(gr.size());
    for (int i_gr = 0; i_gr < gr.size(); ++i_gr) {
      FaFb_rr[i_gr] =
          (Fa.f(i_gr) * mod_Fb.f(i_gr) + Fa.g(i_gr) * mod_Fb.g(i_gr)) /
          (gr.r(i_gr) * gr.r(i_gr));
    }
    return FaFb_rr * mod_Fv;
  } else if (mu == 0) {
    return Coulomb::yk_ab(k, Fa, mod_Fb) * mod_Fv * (1.0 / (2.0 * k + 1));
  } else {
    // Modified spherical Bessel functions
    std::vector<double> i_k(gr.size());
    std::vector<double> k_k(gr.size());

    for (int i_gr = 0; i_gr < gr.size(); ++i_gr) {
      const auto x = mu * r[i_gr];
      i_k[i_gr] = mod_sph_bessel_i(k, x);
      k_k[i_gr] = mod_sph_bessel_k(k, x);
    }

    // Integrate
    std::vector<double> result(gr.size());

    for (int i_mid = 0; i_mid < gr.size(); ++i_mid) {
      double lower_ff =
          NumCalc::integrate(1.0, 0, i_mid, i_k, Fa.f(), mod_Fb.f(), gr.drdu());

      double lower_gg =
          NumCalc::integrate(1.0, 0, i_mid, i_k, Fa.g(), mod_Fb.g(), gr.drdu());

      // For r0 point
      if (i_mid == 0) {
        lower_ff = 0;
        lower_gg = 0;
      }

      const double upper_ff = NumCalc::integrate(1.0, i_mid, gr.size(), k_k,
                                                 Fa.f(), mod_Fb.f(), gr.drdu());

      const double upper_gg = NumCalc::integrate(1.0, i_mid, gr.size(), k_k,
                                                 Fa.g(), mod_Fb.g(), gr.drdu());

      result[i_mid] = (k_k[i_mid] * (lower_ff + lower_gg) +
                       i_k[i_mid] * (upper_ff + upper_gg)) *
                      gr.du();
    }

    return result * mod_Fv;
  }
}

DiracSpinor g0(const DiracSpinor &Fa) {
  DiracSpinor Fb(Fa);
  Fb.g() = (-1.0 * Fa).g();
  return Fb;
}

DiracSpinor g5(const DiracSpinor &Fa) {
  DiracSpinor Fb(Fa);
  Fb.f() = Fa.g();
  Fb.g() = Fa.f();
  return Fb;
}

// Contact limit only, includes VA and AV (zero)
double u_anav_contact(const Wavefunction &wf, const DiracSpinor &Fn,
                      const DiracSpinor &Fa, const DiracSpinor &Fv) {

  const auto gr = wf.grid();
  // Radial part
  std::vector<double> B_aa(gr.size());
  for (int i = 0; i < B_aa.size(); ++i) {
    const double inv_r2 = 1.0 / (gr.r(i) * gr.r(i));
    const double ff_gg = Fa.f(i) * Fa.f(i) + Fa.g(i) * Fa.g(i);
    B_aa[i] = inv_r2 * ff_gg;
  }

  DiracSpinor modFv(Fv);
  modFv.f() = Fv.g();
  modFv.g() = (-1.0 * Fv).f();

  // const double R_anav = Fn * (B_aa * modFv);
  const double R_anav = R_AV_abcd(Fn, Fa, Fv, Fa);

  // Angular part
  const double A_factor =
      (1.0 / (4 * M_PI)) * sqrt(Fa.twojp1() * (1.0 / Fv.twojp1()));
  const double A_cc = Angular::Ck_kk(0, Fa.kappa(), Fa.kappa()) *
                      Angular::Ck_kk(0, Fn.kappa(), -Fv.kappa());

  // Alternatively, calculate full Ck_kk MEs
  double full_cc_VA = 0;
  double full_cc_AV = 0;
  const int tmv = std::min(Fv.twoj(), Fn.twoj());

  for (int k = 0; k <= 0.5 * (Fa.twoj() + Fv.twoj()); ++k) {

    for (int tq = -2 * k; tq <= 2 * k; tq += 2) {
      for (int tma = -Fa.twoj(); tma <= Fa.twoj(); tma += 2) {
        const auto factor = Angular::neg1pow_2(tq) * (2 * k + 1) / (4.0 * M_PI);

        const auto cc =
            factor *
            Angular::Ck_kk_mmq(k, Fn.kappa(), -Fv.kappa(), tmv, tmv, tq) *
            Angular::Ck_kk_mmq(k, Fa.kappa(), Fa.kappa(), tma, tma, -tq);

        full_cc_VA += cc;
      }
    }
  }

  const double spin_factor = 1.0;
  full_cc_VA = full_cc_VA > 1e-10 ? full_cc_VA : 0.0;
  // return spin_factor * A_factor * A_cc * R_anav;
  return full_cc_VA * R_anav;
}

// Contact limit only, includes VA and AV
double u_anva_contact(const Wavefunction &wf, const DiracSpinor &Fn,
                      const DiracSpinor &Fa, const DiracSpinor &Fv) {
  const auto gr = wf.grid();

  // in VA: g5 on 'a' - big contribution
  // In AV: g5 on 'v' - small contribution

  // Radial part
  std::vector<double> B_av_VA(gr.size());
  std::vector<double> B_na_AV(gr.size());
  for (int i = 0; i < gr.size(); ++i) {
    const double inv_r2 = 1.0 / (gr.r(i) * gr.r(i));

    const double ff_gg_VA = Fa.f(i) * Fv.f(i) + Fa.g(i) * Fv.g(i);
    const double ff_gg_AV = Fn.f(i) * Fa.f(i) + Fn.g(i) * Fa.g(i);

    B_av_VA[i] = inv_r2 * ff_gg_VA;
    B_na_AV[i] = inv_r2 * ff_gg_AV;
  }

  DiracSpinor modFa_VA(Fa);
  modFa_VA.f() = Fa.g();
  modFa_VA.g() = (-1.0 * Fa).f();

  DiracSpinor modFv_AV(Fv);
  modFv_AV.f() = Fv.g();
  modFv_AV.g() = (-1.0 * Fv).f();

  // const double R_anva_VA = Fn * (B_av_VA * modFa_VA);
  // const double R_anva_AV = Fa * (B_na_AV * modFv_AV);

  const double R_anva_VA = R_AV_abcd(Fn, Fa, Fa, Fv);
  const double R_anva_AV = R_AV_abcd(Fa, Fn, Fv, Fa);

  // if (Fn.kappa() == -Fv.kappa()) {
  //   std::cout << R_AV_abcd(Fn, Fa, Fa, Fv) << " " << R_AV_abcd(Fa, Fn, Fv, Fa)
  //             << " " << R_AV_abcd(Fn, Fa, Fv, Fa) << "\n";
  // }

  // Angular part
  const double delta = Fn.kappa() == -Fv.kappa() ? 1.0 : 0.0;

  const double A_factor = delta * Angular::neg1pow_2(Fv.twoj() - Fa.twoj()) /
                          (4.0 * M_PI * Fv.twojp1());

  double kA_cc_VA = 0;
  double kA_cc_AV = 0;

  double full_cc_VA = 0;
  double full_cc_AV = 0;
  const int tmv = std::min(Fv.twoj(), Fn.twoj());

  for (int k = 0; k <= Fa.twoj() + Fv.twoj(); ++k) {
    const double twokp1 = 2.0 * k + 1.0;
    kA_cc_VA += twokp1 * Angular::Ck_kk(k, Fn.kappa(), -Fa.kappa()) *
                Angular::Ck_kk(k, Fa.kappa(), Fv.kappa());

    kA_cc_AV += twokp1 * Angular::Ck_kk(k, Fn.kappa(), Fa.kappa()) *
                Angular::Ck_kk(k, Fa.kappa(), -Fv.kappa());

    // Alternatively, calculate full Ck_kk MEs
    for (int tq = -2 * k; tq <= 2 * k; tq += 2) {
      for (int tma = -Fa.twoj(); tma <= Fa.twoj(); tma += 2) {
        const auto factor =
            Angular::neg1pow_2(tq) * (2.0 * k + 1.0) / (4.0 * M_PI);

        const auto cc_VA =
            factor *
            Angular::Ck_kk_mmq(k, Fn.kappa(), -Fa.kappa(), tmv, tma, tq) *
            Angular::Ck_kk_mmq(k, Fa.kappa(), Fv.kappa(), tma, tmv, -tq);

        full_cc_VA += cc_VA;

        const auto cc_AV =
            factor *
            Angular::Ck_kk_mmq(k, Fn.kappa(), Fa.kappa(), tmv, tma, tq) *
            Angular::Ck_kk_mmq(k, Fa.kappa(), -Fv.kappa(), tma, tmv, -tq);

        full_cc_AV += cc_AV;
      }
    }
  }

  full_cc_AV = full_cc_AV > 1e-10 ? full_cc_AV : 0.0;
  full_cc_VA = full_cc_VA > 1e-10 ? full_cc_VA : 0.0;

  // return A_factor * (kA_cc_VA * R_anva_VA + kA_cc_AV * R_anva_AV);
  return full_cc_AV * R_anva_AV + full_cc_VA * R_anva_VA;
}

double R_AV_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
                 const DiracSpinor &Fc, const DiracSpinor &Fd) {
  const auto gr = Fa.grid();
  std::vector<double> neg_fg_gf_inv_r2;
  for (int i = 0; i < gr.size(); ++i) {
    const auto inv_r2 = 1.0 / (gr.r(i) * gr.r(i));
    neg_fg_gf_inv_r2.push_back(1.0 * (Fa.f(i) * Fc.g(i) - Fa.g(i) * Fc.f(i)) *
                               inv_r2);
  }
  return Fb * (neg_fg_gf_inv_r2 * Fd);
}

DiracSpinor old_ig5(const DiracSpinor &Fa) {
  DiracSpinor Fb(Fa);
  Fb.f() = (-1.0 * Fa).g();
  Fb.g() = Fa.f();
  return Fb;
}

DiracSpinor old_i_g0_g5(const DiracSpinor &Fa) {
  // Fb = i*γ5*γ0*Fa

  DiracSpinor Fb(Fa);
  Fb.f() = (-1.0 * Fa).g();
  Fb.g() = (-1.0 * Fa).f();
  return Fb;
}

double mod_sph_bessel_i(double n, double x) {
  gsl_set_error_handler_off();
  gsl_sf_result i_k;
  const int gsl_status = gsl_sf_bessel_Inu_e(n + 0.5, x, &i_k);

  if (gsl_status == GSL_SUCCESS) {
    /*if (i_k.err / i_k.val > 0.01) {
      std::cout << "\nWARNING: error in i_k greater than 1\%\ni_k = " << i_k.val
                << " \u00b1 " << i_k.err << "\n";
    }*/
    return std::sqrt(M_PI / (2.0 * x)) * i_k.val;
  } else if (gsl_status == GSL_EOVRFLW) {
    return 0.0;
  }
  std::cout << "Need GSL_SUCCESS = " << GSL_SUCCESS
            << " or GSL_EOVRFLW = " << GSL_EOVRFLW;
  std::cout << "\n\ngsl_status = " << gsl_status << "\n";

  throw std::bad_function_call();
}

double mod_sph_bessel_k(double n, double x) {
  gsl_set_error_handler_off();
  gsl_sf_result k_k;
  // std::cout << "test..." << std::flush;
  const int gsl_status = gsl_sf_bessel_Knu_e(n + 0.5, x, &k_k);
  // std::cout << "...success";

  if (gsl_status == GSL_SUCCESS) {
    /*
    if (k_k.err / k_k.val > 0.01) {
      std::cout << "\nWARNING: error in i_k greater than 1\%\n " << k_k.val
                << " \u00b1 " << k_k.err << "\n";
    }*/
    return std::sqrt(2.0 / (M_PI * x)) * k_k.val;
  } else if (gsl_status == GSL_EUNDRFLW) {
    return 0.0;
  }
  std::cout << "Need GSL_SUCCESS = " << GSL_SUCCESS
            << " or GSL_EUNDRFLW = " << GSL_EUNDRFLW;
  std::cout << "\n\ngsl_status = " << gsl_status << "\n";

  throw std::bad_function_call();
}

} // namespace BSM_Vee