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

  const auto s = Angular::neg1pow_2(v.twoj() + r.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto sab = Angular::neg1pow_2(b.twoj() + a.twoj());
      const auto inv_e_rwab = 1.0 / (r.en + w.en - a.en - b.en);

      const auto [minU, maxU] = mY.k_minmax_W(w, r, b, a);
      const auto [minL, maxL] = mY.k_minmax_Q(v, a, b, c);
      for (int u = minU; u <= maxU; ++u) {

        const auto wu = mY.Wk(u, w, r, b, a);
        if (Angular::zeroQ(wu))
          continue;

        for (int l = minL; l <= maxL; ++l) {
          const auto sj1 = mY.sixj(w, v, K, l, u, b);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, a);
          if (Angular::zeroQ(sj2))
            continue;

          const auto ql = mY.Qk(l, v, a, b, c); // *

          t += sab * sj1 * sj2 * wu * ql * inv_e_rwab;
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

  const auto s = Angular::neg1pow_2(w.twoj() - c.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_rnav = 1.0 / (r.en + n.en - a.en - v.en);

      const auto [minU1, maxU1] = mY.k_minmax_W(w, a, c, n);
      const auto [minU2, maxU2] = mY.k_minmax_W(v, a, r, n);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.sixj(w, v, K, r, c, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto su = Angular::neg1pow(u);
        const auto f = (2 * u + 1);

        const auto wu1 = mY.Wk(u, w, a, c, n);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = mY.Wk(u, v, a, r, n);

        t += su * sj * wu1 * wu2 * inv_e_rnav / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::t3(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - c.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_nwac = 1.0 / (n.en + w.en - a.en - c.en);

      const auto [minU1, maxU1] = mY.k_minmax_W(w, n, c, a);
      const auto [minU2, maxU2] = mY.k_minmax_W(v, n, r, a);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.sixj(w, v, K, r, c, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto su = Angular::neg1pow(u);
        const auto f = (2 * u + 1);

        const auto wu1 = mY.Wk(u, w, n, c, a);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = mY.Wk(u, v, n, r, a);

        t += su * sj * wu1 * wu2 * inv_e_nwac / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::t4(const int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(v.twoj() + r.twoj() + 2 * K);

  for (const auto &n : mExcited) {
    for (const auto &m : mExcited) {

      const auto inv_e_nmcv = 1.0 / (n.en + m.en - c.en - v.en);

      const auto [minU, maxU] = mY.k_minmax_Q(w, r, n, m);
      const auto [minL, maxL] = mY.k_minmax_W(n, m, v, c);
      for (int u = minU; u <= maxU; ++u) {

        const auto snm = Angular::neg1pow_2(n.twoj() + m.twoj());

        const auto qu = mY.Qk(u, w, r, n, m);
        if (Angular::zeroQ(qu))
          continue;

        for (int l = minL; l <= maxL; ++l) {

          const auto sj1 = mY.sixj(w, v, K, l, u, n);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          const auto wl = mY.Wk(l, n, m, v, c);

          t += snm * sj1 * sj2 * qu * wl * inv_e_nmcv;
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

  const auto s = Angular::neg1pow_2(v.twoj() + c.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto inv_e_rvab = 1.0 / (r.en + v.en - a.en - b.en);

      const auto [minU, maxU] = mY.k_minmax_W(v, r, b, a);
      const auto [minL, maxL] = mY.k_minmax_Q(w, c, b, a);
      for (int u = minU; u <= maxU; ++u) {

        const auto sab = Angular::neg1pow_2(a.twoj() + b.twoj());

        const auto wu = mY.Wk(u, v, r, b, a);
        if (Angular::zeroQ(wu))
          continue;

        for (int l = minL; l <= maxL; ++l) {
          const auto sj1 = mY.sixj(w, v, K, u, l, b);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, a);
          if (Angular::zeroQ(sj2))
            continue;

          const auto ql = mY.Qk(l, w, c, b, a);

          t += sab * sj1 * sj2 * wu * ql * inv_e_rvab;
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

  const auto s = Angular::neg1pow_2(w.twoj() - r.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_rnaw = 1.0 / (r.en + n.en - a.en - w.en);

      const auto [minU1, maxU1] = mY.k_minmax_W(v, a, c, n);
      const auto [minU2, maxU2] = mY.k_minmax_W(w, a, r, n);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.sixj(w, v, K, c, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto su = Angular::neg1pow(u);
        const auto f = (2 * u + 1);

        const auto wu1 = mY.Wk(u, v, a, c, n);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = mY.Wk(u, w, a, r, n);

        t += su * sj * wu1 * wu2 * inv_e_rnaw / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b3(const int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - r.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_nvac = 1.0 / (n.en + v.en - a.en - c.en);

      const auto [minU1, maxU1] = mY.k_minmax_W(v, n, c, a);
      const auto [minU2, maxU2] = mY.k_minmax_W(w, n, r, a);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.sixj(w, v, K, c, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto su = Angular::neg1pow(u);
        const auto f = (2 * u + 1);

        const auto wu1 = mY.Wk(u, v, n, c, a);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = mY.Wk(u, w, n, r, a);

        t += su * sj * wu1 * wu2 * inv_e_nvac / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b4(const int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(v.twoj() + c.twoj() + 2 * K);

  for (const auto &n : mExcited) {
    for (const auto &m : mExcited) {

      const auto inv_e_nmcw = 1.0 / (n.en + m.en - c.en - w.en);

      const auto [minU, maxU] = mY.k_minmax_Q(v, m, n, r);
      const auto [minL, maxL] = mY.k_minmax_W(w, c, n, m);
      for (int u = minU; u <= maxU; ++u) {

        const auto snm = Angular::neg1pow_2(n.twoj() + m.twoj());

        const auto qu = mY.Qk(u, v, m, n, r);
        if (Angular::zeroQ(qu))
          continue;

        for (int l = minL; l <= maxL; ++l) {

          const auto sj1 = mY.sixj(w, v, K, u, l, n);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          const auto wl = mY.Wk(l, w, c, n, m);

          t += snm * sj1 * sj2 * qu * wl * inv_e_nmcw;
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

  const auto s = Angular::neg1pow_2(w.twoj() - c.twoj() + 2 * K);

  for (const auto &b : mCore) {
    for (const auto &n : mExcited) {

      const auto e_nvbc = n.en + v.en - b.en - c.en;
      const auto e_nwab = n.en + w.en - a.en - b.en;
      const auto invde = 1.0 / (e_nvbc * e_nwab);

      const auto [minU1, maxU1] = mY.k_minmax_W(w, n, a, b);
      const auto [minU2, maxU2] = mY.k_minmax_W(v, n, c, b);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto su = Angular::neg1pow(u);
        const auto f = 2 * u + 1;

        const auto sj = mY.sixj(w, v, K, c, a, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto wu1 = mY.Wk(u, w, n, a, b);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = mY.Wk(u, v, n, c, b);

        t += su * sj * wu1 * wu2 * invde / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::c2(const int K, const DiracSpinor &w, const DiracSpinor &a,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(v.twoj() + a.twoj() + 2 * K);

  for (const auto &m : mExcited) {
    for (const auto &n : mExcited) {

      const auto e_mnaw = m.en + n.en - a.en - w.en;
      const auto e_nmcv = n.en + m.en - c.en - v.en;
      const auto invde = 1.0 / (e_mnaw * e_nmcv);
      const auto smn = Angular::neg1pow_2(m.twoj() + n.twoj());

      const auto [minU, maxU] = mY.k_minmax_Q(v, m, n, c);
      const auto [minL, maxL] = mY.k_minmax_W(w, a, n, m);
      for (int u = minU; u <= maxU; ++u) {

        const auto qu = mY.Qk(u, v, m, n, c);
        if (Angular::zeroQ(qu))
          continue;

        for (int l = minL; l <= maxL; ++l) {

          const auto sj1 = mY.sixj(w, v, K, u, l, n);
          if (Angular::zeroQ(sj1))
            continue;

          const auto sj2 = mY.sixj(c, a, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          const auto wl = mY.Wk(l, w, a, n, m);

          t += smn * sj1 * sj2 * qu * wl * invde;
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

  const auto s = Angular::neg1pow_2(v.twoj() + r.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto sab = Angular::neg1pow_2(a.twoj() + b.twoj());

      const auto e_mvab = m.en + v.en - a.en - b.en;
      const auto e_rwba = r.en + w.en - b.en - a.en;
      const auto invde = 1.0 / (e_mvab * e_rwba);

      const auto [minU, maxU] = mY.k_minmax_Q(v, b, a, m);
      const auto [minL, maxL] = mY.k_minmax_W(w, r, a, b);
      for (int u = minU; u <= maxU; ++u) {

        const auto qu = mY.Qk(u, v, b, a, m);
        if (Angular::zeroQ(qu))
          continue;

        for (int l = minL; l <= maxL; ++l) {

          const auto sj1 = mY.sixj(w, v, K, u, l, a);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(m, r, K, l, u, b);
          if (Angular::zeroQ(sj2))
            continue;

          const auto wl = mY.Wk(l, w, r, a, b);

          t += sab * sj1 * sj2 * qu * wl * invde;
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

  const auto s = Angular::neg1pow_2(w.twoj() - m.twoj() + 2 * K);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto e_mnav = m.en + n.en - a.en - v.en;
      const auto e_rnaw = r.en + n.en - a.en - w.en;
      const auto invde = 1.0 / (e_mnav * e_rnaw);

      const auto [minU1, maxU1] = mY.k_minmax_W(w, a, r, n);
      const auto [minU2, maxU2] = mY.k_minmax_W(v, a, m, n);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto su = Angular::neg1pow(u);
        const auto f = 2 * u + 1;

        const auto sj = mY.sixj(w, v, K, m, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto wu1 = mY.Wk(u, w, a, r, n);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = mY.Wk(u, v, a, m, n);

        t += su * sj * wu1 * wu2 * invde / f;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::n1(const DiracSpinor &v) const {

  double t = 0.0;

  const auto tjvp1 = v.twoj() + 1;

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {
      for (const auto &m : mExcited) {

        const auto e_vmab = (v.en + m.en - a.en - b.en);
        const auto ide2 = 1.0 / (e_vmab * e_vmab);

        const auto [minU, maxU] = mY.k_minmax_W(v, m, a, b);
        for (int u = minU; u <= maxU; ++u) {

          const auto tup1 = 2 * u + 1;

          const auto q = mY.Qk(u, v, m, a, b);
          if (Angular::zeroQ(q))
            continue;

          const auto w = q + mY.Pk(u, v, m, a, b);

          t += q * w * ide2 / tup1;
        }
      }
    }
  }
  return t / tjvp1;
}

//******************************************************************************
double StructureRad::n2(const DiracSpinor &v) const {

  double t = 0.0;

  const auto tjvp1 = v.twoj() + 1;

  for (const auto &a : mCore) {
    for (const auto &m : mExcited) {
      for (const auto &n : mExcited) {

        const auto e_nmav = (n.en + m.en - a.en - v.en);
        const auto ide2 = 1.0 / (e_nmav * e_nmav);

        const auto [minU, maxU] = mY.k_minmax_W(v, a, n, m);
        for (int u = minU; u <= maxU; ++u) {

          const auto tup1 = 2 * u + 1;

          const auto q = mY.Qk(u, v, a, n, m);
          if (Angular::zeroQ(q))
            continue;

          const auto w = q + mY.Pk(u, v, a, n, m);

          t += q * w * ide2 / tup1;
        }
      }
    }
  }
  return t / tjvp1;
}
} // namespace MBPT
