#include "StructureRad.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "HF/ExternalField.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <numeric>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
constexpr bool use_omp = true;
#else
constexpr bool use_omp = true;
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

namespace MBPT {

//******************************************************************************
StructureRad::StructureRad(const std::vector<DiracSpinor> &basis,
                           double en_core, std::pair<int, int> nminmax)
    : mBasis(basis), // store a local copy
      mY(basis.front().rgrid, &mBasis) {

  // nb: en_core defined such that: Fa.en < en_core ==> core state!

  // Sort into core/excited; store pointers
  // nb: this makes it faster..
  const auto [n_min, n_max] = nminmax;
  for (const auto &Fn : mBasis) {
    if (Fn.en < en_core && Fn.n >= n_min) {
      mCore.push_back(Fn);
    } else if (Fn.en > en_core && Fn.n <= n_max) {
      mExcited.push_back(Fn);
    }
  }
}

//******************************************************************************
std::pair<double, double>
StructureRad::srTB(const DiracOperator::TensorOperator *const h,
                   const DiracSpinor &w, const DiracSpinor &v, double omega,
                   const HF::ExternalField *const dV) const {

  if (h->isZero(w.k, v.k))
    return {0.0, 0.0};

  const auto k = h->rank();

  const std::size_t num_para_threads =
      use_omp ? std::size_t(2 * omp_get_max_threads()) : 1;

  std::vector<double> sr(num_para_threads);
  std::vector<double> sr_dv(num_para_threads);

#pragma omp parallel for num_threads(num_para_threads)
  for (auto ir = 0ul; ir < mExcited.size(); ++ir) {
    const auto &r = mExcited[ir];
    const auto tid = std::size_t(omp_get_thread_num());
    for (const auto &c : mCore) {

      if (h->isZero(c.k, r.k))
        continue;

      const auto inv_erc_pw = 1.0 / (r.en - c.en + omega);
      const auto inv_erc_mw = 1.0 / (r.en - c.en - omega);

      const auto t_cr = h->reducedME(c, r);
      const auto t_rc = h->reducedME(r, c); // blah

      const auto T_wrvc = t1(k, w, r, v, c) + t2(k, w, r, v, c) +
                          t3(k, w, r, v, c) + t4(k, w, r, v, c);

      // Is this correct?? (i.e., T=B if v=w), up to energy denom..
      const auto B_wcvr = v == w ? T_wrvc
                                 : b1(k, w, c, v, r) + b2(k, w, c, v, r) +
                                       b3(k, w, c, v, r) + b4(k, w, c, v, r);

      sr[tid] += (t_cr * T_wrvc * inv_erc_pw) + (t_rc * B_wcvr * inv_erc_mw);

      if (dV) {
        const auto tdv_cr = t_cr + dV->dV(c, r);
        const auto tdv_rc = t_rc + dV->dV(r, c); // blah
        sr_dv[tid] +=
            (tdv_cr * T_wrvc * inv_erc_pw) + (tdv_rc * B_wcvr * inv_erc_mw);
      }
    }
  }

  // nb: this just a test.. return (t+b) in final version
  const auto tb = std::accumulate(cbegin(sr), cend(sr), 0.0);
  const auto dv = std::accumulate(cbegin(sr_dv), cend(sr_dv), 0.0);

  return {tb, dv};
}

//******************************************************************************
std::pair<double, double>
StructureRad::srC(const DiracOperator::TensorOperator *const h,
                  const DiracSpinor &w, const DiracSpinor &v,
                  const HF::ExternalField *const dV) const {

  if (h->isZero(w.k, v.k))
    return {0.0, 0.0};

  const auto k = h->rank();

  // For parallelisation:
  const std::size_t num_para_threads =
      use_omp ? std::size_t(2 * omp_get_max_threads()) : 1;

  std::vector<double> src(num_para_threads);
  std::vector<double> src_dv(num_para_threads);

#pragma omp parallel for num_threads(num_para_threads) collapse(2)
  for (auto ia = 0ul; ia < mCore.size(); ++ia) {
    for (auto ic = 0ul; ic < mCore.size(); ++ic) {
      const auto &a = mCore[ia];
      const auto &c = mCore[ic];

      if (h->isZero(c.k, a.k))
        continue;

      const auto t_ca = h->reducedME(c, a);
      const auto C_wavc = c1(k, w, a, v, c) + c2(k, w, a, v, c);

      // nb: -ve
      const auto tid = std::size_t(omp_get_thread_num());
      src[tid] -= t_ca * C_wavc;

      if (dV) {
        const auto tdv_ca = t_ca + dV->dV(c, a);
        src_dv[tid] -= tdv_ca * C_wavc;
      }
    }
  }

// collapse(2)
#pragma omp parallel for num_threads(num_para_threads)
  for (auto im = 0ul; im < mExcited.size(); ++im) {
    const auto &m = mExcited[im];
    const auto tid = std::size_t(omp_get_thread_num());
    for (const auto &r : mExcited) {

      if (h->isZero(m.k, r.k))
        continue;

      const auto t_mr = h->reducedME(m, r);
      const auto C_wrvm = c3(k, w, r, v, m) + c4(k, w, r, v, m);

      // nb: -ve
      src[tid] -= t_mr * C_wrvm;

      if (dV) {
        const auto tdv_mr = t_mr + dV->dV(m, r);
        src_dv[tid] -= tdv_mr * C_wrvm;
      }
    }
  }

  const auto t = std::accumulate(cbegin(src), cend(src), 0.0);
  const auto dv = std::accumulate(cbegin(src_dv), cend(src_dv), 0.0);

  return {t, dv};
}

//******************************************************************************
std::pair<double, double>
StructureRad::norm(const DiracOperator::TensorOperator *const h,
                   const DiracSpinor &w, const DiracSpinor &v,
                   const HF::ExternalField *const dV) const {

  if (h->isZero(w.k, v.k))
    return {0.0, 0.0};

  const auto t_wv = h->reducedME(w, v);
  const auto tdv_wv = dV ? t_wv + dV->dV(w, v) : 0.0;

  const auto nv = n1(v) + n2(v);
  const auto nw = w == v ? nv : n1(w) + n2(w);

  return {-0.5 * t_wv * (nv + nw), -0.5 * tdv_wv * (nv + nw)};
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
double StructureRad::t1(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  // (-1)^{w-v+K}
  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto inv_e_rwab = 1.0 / (r.en + w.en - a.en - b.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_rwab) < 1.0 / e_thresh))
      //   continue;

      const auto [minL, maxL] = mY.k_minmax(a, c);

      for (int u = 0; u <= max_u; ++u) {

        const auto z = mY.Zk(u, w, r, b, a);
        if (Angular::zeroQ(z))
          continue;

        for (int l = minL; l <= maxL; ++l) {
          const auto sj1 = mY.sixj(w, v, K, l, u, b);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, a);
          if (Angular::zeroQ(sj2))
            continue;

          const auto x = mY.Xk(l, b, a, v, c); // *
          if (Angular::zeroQ(x))
            continue;

          t += sj1 * sj2 * z * x * inv_e_rwab;
        }
      }
    }
  }
  return s * t;

  // * note: It would seem that it should be faster to calculate x
  // and store in array first. Turns out, this is not the case.
  // Probably, because we take advantage of selection rules, and only calc x
  // when the other terms are all non-zero!
}

//******************************************************************************
double StructureRad::t2(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_rnav = 1.0 / (r.en + n.en - a.en - v.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_rnav) < 1.0 / e_thresh))
      //   continue;

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, r, c, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2 = Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u);
        const auto f = (2 * u + 1);

        const auto z1 = mY.Zk(u, a, w, n, c);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, n, r, a, v);
        if (Angular::zeroQ(z2))
          continue;

        t += s2 * sj * z1 * z2 * inv_e_rnav / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::t3(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_nwac = 1.0 / (n.en + w.en - a.en - c.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_nwac) < 1.0 / e_thresh))
      //   continue;

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, r, c, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2 = Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u);
        const auto f = (2 * u + 1);

        const auto z1 = mY.Zk(u, n, w, a, c);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, a, r, n, v);
        if (Angular::zeroQ(z2))
          continue;

        t += s2 * sj * z1 * z2 * inv_e_nwac / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::t4(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_l = mY.Ck().max_k();

  for (const auto &n : mExcited) {
    for (const auto &m : mExcited) {

      const auto inv_e_nmcv = 1.0 / (n.en + m.en - c.en - v.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_nmcv) < 1.0 / e_thresh))
      //   continue;

      const auto [minU, maxU] = mY.k_minmax(r, m);

      for (int u = minU; u <= maxU; ++u) {

        const auto x = mY.Xk(u, w, r, n, m);
        if (Angular::zeroQ(x))
          continue;

        for (int l = 0; l <= max_l; ++l) {

          const auto z = mY.Zk(l, n, m, v, c);
          if (Angular::zeroQ(z))
            continue;

          const auto sj1 = mY.sixj(w, v, K, l, u, n);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          t += sj1 * sj2 * x * z * inv_e_nmcv;
        }
      }
    }
  }
  return s * t;
}

//******************************************************************************
//******************************************************************************
double StructureRad::b1(const int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto inv_e_rvab = 1.0 / (r.en + v.en - a.en - b.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_rvab) < 1.0 / e_thresh))
      //   continue;

      const auto [minL, maxL] = mY.k_minmax(c, a);

      for (int u = 0; u <= max_u; ++u) {

        const auto z = mY.Zk(u, b, a, v, r);
        if (Angular::zeroQ(z))
          continue;

        for (int l = minL; l <= maxL; ++l) {
          const auto sj1 = mY.sixj(w, v, K, u, l, b);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, a);
          if (Angular::zeroQ(sj2))
            continue;

          const auto x = mY.Xk(l, w, c, b, a);
          if (Angular::zeroQ(x))
            continue;

          t += sj1 * sj2 * z * x * inv_e_rvab;
        }
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b2(const int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_rnaw = 1.0 / (r.en + n.en - a.en - w.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_rnaw) < 1.0 / e_thresh))
      //   continue;

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, c, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2 = Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u);
        const auto f = (2 * u + 1);

        const auto z1 = mY.Zk(u, n, c, a, v);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, a, w, n, r);
        if (Angular::zeroQ(z2))
          continue;

        t += s2 * sj * z1 * z2 * inv_e_rnaw / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b3(const int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_nvac = 1.0 / (n.en + v.en - a.en - c.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_nvac) < 1.0 / e_thresh))
      //   continue;

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, c, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2 = Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u);
        const auto f = (2 * u + 1);

        const auto z1 = mY.Zk(u, a, c, n, v);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, n, w, a, r);
        if (Angular::zeroQ(z2))
          continue;

        t += s2 * sj * z1 * z2 * inv_e_nvac / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b4(const int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_l = mY.Ck().max_k();

  for (const auto &n : mExcited) {
    for (const auto &m : mExcited) {

      const auto inv_e_nmcw = 1.0 / (n.en + m.en - c.en - w.en);

      // if (e_thresh > 0.0 && (std::abs(inv_e_nmcw) < 1.0 / e_thresh))
      //   continue;

      const auto [minU, maxU] = mY.k_minmax(m, r);

      for (int u = minU; u <= maxU; ++u) {
        const auto x = mY.Xk(u, n, m, v, r);
        if (Angular::zeroQ(x))
          continue;

        for (int l = 0; l <= max_l; ++l) {

          const auto z = mY.Zk(l, w, c, n, m);
          if (Angular::zeroQ(z))
            continue;

          const auto sj1 = mY.sixj(w, v, K, u, l, n);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          t += sj1 * sj2 * x * z * inv_e_nmcw;
        }
      }
    }
  }
  return s * t;
}

//******************************************************************************
//******************************************************************************
double StructureRad::c1(const int K, const DiracSpinor &w, const DiracSpinor &a,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto max_u = mY.Ck().max_k();

  for (const auto &b : mCore) {
    for (const auto &n : mExcited) {

      const auto e_nvbc = n.en + v.en - b.en - c.en;
      const auto e_nwab = n.en + w.en - a.en - b.en;
      // if (e_thresh > 0.0 &&
      //     (std::abs(e_nvbc) > e_thresh || std::abs(e_nwab) > e_thresh))
      //   continue;

      const auto invde = 1.0 / (e_nvbc * e_nwab);

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, c, a, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto z1 = mY.Zk(u, w, n, a, b);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, b, c, n, v);
        if (Angular::zeroQ(z2))
          continue;

        const auto s = Angular::neg1pow_2(b.twoj() - n.twoj() + 2 * (K + u));
        const auto f = 2 * u + 1;

        t += s * sj * z1 * z2 * invde / f;
      }
    }
  }
  return t;
}

//******************************************************************************
double StructureRad::c2(const int K, const DiracSpinor &w, const DiracSpinor &a,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_l = mY.Ck().max_k();

  for (const auto &m : mExcited) {
    for (const auto &n : mExcited) {

      const auto e_mnaw = m.en + n.en - a.en - w.en;
      const auto e_nmcv = n.en + m.en - c.en - v.en;
      // if (e_thresh > 0.0 &&
      //     (std::abs(e_mnaw) > e_thresh || std::abs(e_nmcv) > e_thresh))
      //   continue;

      const auto invde = 1.0 / (e_mnaw * e_nmcv);

      const auto [minU, maxU] = mY.k_minmax(n, v);
      for (int u = minU; u <= maxU; ++u) {
        for (int l = 0; l <= max_l; ++l) {

          const auto sj1 = mY.sixj(w, v, K, u, l, n);
          if (Angular::zeroQ(sj1))
            continue;

          const auto sj2 = mY.sixj(c, a, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          const auto x = mY.Xk(u, m, n, c, v);
          if (Angular::zeroQ(x))
            continue;

          const auto z = mY.Zk(l, a, w, m, n);
          if (Angular::zeroQ(z))
            continue;

          t += sj1 * sj2 * x * z * invde;
        }
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::c3(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &m) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_l = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto e_mvab = m.en + v.en - a.en - b.en;
      const auto e_rwba = r.en + w.en - b.en - a.en;
      // if (e_thresh > 0.0 &&
      //     (std::abs(e_mvab) > e_thresh || std::abs(e_rwba) > e_thresh))
      //   continue;

      const auto invde = 1.0 / (e_mvab * e_rwba);

      const auto [minU, maxU] = mY.k_minmax(b, m);
      for (int u = minU; u <= maxU; ++u) {
        for (int l = 0; l <= max_l; ++l) {

          const auto sj1 = mY.sixj(w, v, K, u, l, a);
          if (Angular::zeroQ(sj1))
            continue;

          const auto sj2 = mY.sixj(m, r, K, l, u, b);
          if (Angular::zeroQ(sj2))
            continue;

          const auto x = mY.Xk(u, a, b, v, m);
          if (Angular::zeroQ(x))
            continue;

          const auto z = mY.Zk(l, r, w, b, a);
          if (Angular::zeroQ(z))
            continue;

          t += sj1 * sj2 * x * z * invde;
        }
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::c4(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &m) const {

  double t = 0.0;

  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto e_mnav = m.en + n.en - a.en - v.en;
      const auto e_rnaw = r.en + n.en - a.en - w.en;
      // if (e_thresh > 0.0 &&
      //     (std::abs(e_mnav) > e_thresh || std::abs(e_rnaw) > e_thresh))
      //   continue;

      const auto invde = 1.0 / (e_mnav * e_rnaw);

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, m, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto z1 = mY.Zk(u, a, w, n, r);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, m, n, v, a);
        if (Angular::zeroQ(z2))
          continue;

        const auto s = Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * (K + u));
        const auto f = 2 * u + 1;

        t += s * sj * z1 * z2 * invde / f;
      }
    }
  }
  return t;
}

//******************************************************************************
double StructureRad::n1(const DiracSpinor &v) const {
  //
  const auto tjvp1 = v.twoj() + 1;

  double t = 0.0;

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {
      for (const auto &m : mExcited) {

        const auto e_vmab = (v.en + m.en - a.en - b.en);
        const auto ide2 = 1.0 / (e_vmab * e_vmab);

        const auto [minU, maxU] = mY.k_minmax(m, b);
        for (int u = minU; u <= maxU; ++u) {

          const auto tup1 = 2 * u + 1;

          const auto x = mY.Xk(u, v, m, a, b);
          if (Angular::zeroQ(x))
            continue;

          // XXX Note: Break Z into Q+P, saves time!!!
          const auto z = mY.Zk(u, v, m, a, b);

          t += x * z * ide2 / tup1;
        }
      }
    }
  }
  return t / tjvp1;
}

//******************************************************************************
double StructureRad::n2(const DiracSpinor &v) const {
  //
  const auto tjvp1 = v.twoj() + 1;

  double t = 0.0;

  for (const auto &a : mCore) {
    for (const auto &m : mExcited) {
      for (const auto &n : mExcited) {

        const auto e_nmav = (n.en + m.en - a.en - v.en);
        const auto ide2 = 1.0 / (e_nmav * e_nmav);

        const auto [minU, maxU] = mY.k_minmax(a, m);
        for (int u = minU; u <= maxU; ++u) {

          const auto tup1 = 2 * u + 1;

          const auto x = mY.Xk(u, v, a, n, m);
          if (Angular::zeroQ(x))
            continue;

          // XXX Note: Break Z into Q+P, saves time!!!
          const auto z = mY.Zk(u, v, a, n, m);

          t += x * z * ide2 / tup1;
        }
      }
    }
  }
  return t / tjvp1;
}
} // namespace MBPT
