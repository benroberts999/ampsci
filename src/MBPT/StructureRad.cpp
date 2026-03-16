#include "StructureRad.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/omp.hpp"
#include <numeric>
#include <vector>

namespace MBPT {

//==============================================================================
StructureRad::StructureRad(const std::vector<DiracSpinor> &basis,
                           double en_core, std::pair<int, int> nminmax,
                           const std::string &Qk_fname)
    : m_use_Qk(!Qk_fname.empty()) {
  // nb: en_core defined such that: Fa.en() < en_core ==> core state!

  // nb: this makes it faster..
  const auto [n_min, n_max] = nminmax;
  for (const auto &Fn : basis) {
    if (Fn.en() < en_core && Fn.n() >= n_min) {
      mCore.push_back(Fn);
    } else if (Fn.en() > en_core && Fn.n() <= n_max) {
      mExcited.push_back(Fn);
    }
  }
  auto both = mCore;
  both.insert(both.end(), mExcited.begin(), mExcited.end());
  mY.calculate(both);

  // nb: don't need mY if using QkTable.
  // However, require SixJTable!

  // std::string Qfname;
  if (m_use_Qk) {
    std::cout << "\nFill Qk table:\n";
    mQ = Coulomb::QkTable{};
    const auto ok = mQ->read(Qk_fname);
    if (!ok) {
      mQ->fill(both, mY);
      mQ->write(Qk_fname);
    }
  }

  // set index lists, for faster
}

//==============================================================================
std::pair<double, double>
StructureRad::srTB(const DiracOperator::TensorOperator *const h,
                   const DiracSpinor &w, const DiracSpinor &v, double omega,
                   const ExternalField::CorePolarisation *const dV) const {
  if (h->isZero(w.kappa(), v.kappa()))
    return {0.0, 0.0};

  const auto k = h->rank();

  std::vector<std::pair<std::size_t, std::size_t>> index_ra;
  for (auto ir = 0ul; ir < mExcited.size(); ++ir) {
    const auto &r = mExcited[ir];
    for (auto ia = 0ul; ia < mCore.size(); ++ia) {
      const auto &a = mCore[ia];
      if (h->isZero(a.kappa(), r.kappa()))
        continue;
      index_ra.emplace_back(ir, ia);
    }
  }

  double tb{0.0}, dv{0.0};
#pragma omp parallel for reduction(+ : tb) reduction(+ : dv)
  for (std::size_t i = 0; i < index_ra.size(); ++i) {
    const auto [ir, ia] = index_ra[i];
    const auto &r = mExcited[ir];
    const auto &a = mCore[ia];

    const auto inv_era_pw = 1.0 / (r.en() - a.en() + omega);
    const auto inv_era_mw = 1.0 / (r.en() - a.en() - omega);

    const auto t_ar = h->reducedME(a, r);
    const auto t_ra = h->symm_sign(a, r) * t_ar;

    const auto T_wrva = t1234(k, w, r, v, a);
    const auto B_wavr = v == w ? T_wrva : b1234(k, w, a, v, r);

    tb += (t_ar * T_wrva * inv_era_pw) + (t_ra * B_wavr * inv_era_mw);

    if (dV) {
      const auto dVar = dV->dV(a, r);
      const auto tdv_ar = t_ar + dVar;
      const auto tdv_ra = h->symm_sign(a, r) * tdv_ar;
      dv += (tdv_ar * T_wrva * inv_era_pw) + (tdv_ra * B_wavr * inv_era_mw);
    }
  }

  return {tb, dv};
}

//==============================================================================
std::pair<double, double>
StructureRad::srC(const DiracOperator::TensorOperator *const h,
                  const DiracSpinor &w, const DiracSpinor &v,
                  const ExternalField::CorePolarisation *const dV) const {
  if (h->isZero(w.kappa(), v.kappa()))
    return {0.0, 0.0};

  const auto k = h->rank();

  std::vector<std::pair<std::size_t, std::size_t>> index_ab;
  for (auto ia = 0ul; ia < mCore.size(); ++ia) {
    const auto &a = mCore[ia];
    for (auto ib = 0ul; ib < mCore.size(); ++ib) {
      const auto &b = mCore[ib];
      if (h->isZero(a.kappa(), b.kappa()))
        continue;
      index_ab.emplace_back(ia, ib);
    }
  }

  std::vector<std::pair<std::size_t, std::size_t>> index_mr;
  for (auto im = 0ul; im < mExcited.size(); ++im) {
    const auto &m = mExcited[im];
    for (auto ir = 0ul; ir < mExcited.size(); ++ir) {
      const auto &r = mExcited[ir];
      if (h->isZero(m.kappa(), r.kappa()))
        continue;
      index_mr.emplace_back(im, ir);
    }
  }

  double c{0.0}, dv{0.0};
#pragma omp parallel for reduction(- : c) reduction(- : dv)
  for (std::size_t i = 0; i < index_ab.size(); ++i) {
    const auto [ia, ib] = index_ab[i];
    const auto &a = mCore[ia];
    const auto &b = mCore[ib];

    const auto t_ba = h->reducedME(b, a);
    const auto C_wavb = c1(k, w, a, v, b) + c2(k, w, a, v, b);

    // nb: -ve
    c -= t_ba * C_wavb;

    if (dV) {
      const auto tdv_ba = t_ba + dV->dV(b, a);
      dv -= tdv_ba * C_wavb;
    }
  }

#pragma omp parallel for reduction(- : c) reduction(- : dv)
  for (std::size_t i = 0; i < index_mr.size(); ++i) {
    const auto [im, ir] = index_mr[i];
    const auto &m = mExcited[im];
    const auto &r = mExcited[ir];

    const auto t_mr = h->reducedME(m, r);
    const auto C_wrvm = d2(k, w, r, v, m) + d1(k, w, r, v, m);

    // nb: -ve
    c -= t_mr * C_wrvm;

    if (dV) {
      const auto tdv_mr = t_mr + dV->dV(m, r);
      dv -= tdv_mr * C_wrvm;
    }
  }

  return {c, dv};
}

//==============================================================================
std::pair<double, double>
StructureRad::norm(const DiracOperator::TensorOperator *const h,
                   const DiracSpinor &w, const DiracSpinor &v,
                   const ExternalField::CorePolarisation *const dV) const {
  if (h->isZero(w.kappa(), v.kappa()))
    return {0.0, 0.0};

  const auto t_wv = h->reducedME(w, v);
  const auto tdv_wv = dV ? t_wv + dV->dV(w, v) : 0.0;

  const auto nv = n1(v) + n2(v);
  const auto nw = w == v ? nv : n1(w) + n2(w);

  return {-0.5 * t_wv * (nv + nw), -0.5 * tdv_wv * (nv + nw)};
}

//==============================================================================
Coulomb::meTable<std::pair<double, double>>
StructureRad::srn_table(const DiracOperator::TensorOperator *const h,
                        const std::vector<DiracSpinor> &as,
                        const std::vector<DiracSpinor> &tbs, double omega,
                        const ExternalField::CorePolarisation *const dV) const {
  const auto &bs = tbs.empty() ? as : tbs;

  Coulomb::meTable<std::pair<double, double>> tab;

  // First, fill with zeros. Must be done in serial
  for (const auto &a : as) {
    for (const auto &b : bs) {
      if (h->isZero(a, b))
        continue;
      tab.add(a, b, {0.0, 0.0});
      tab.add(b, a, {0.0, 0.0});
    }
  }

#pragma omp parallel for collapse(2)
  for (std::size_t ia = 0; ia < as.size(); ++ia) {
    for (std::size_t ib = 0; ib < bs.size(); ++ib) {
      const auto &a = as[ia];
      const auto &b = bs[ib];
      if (h->isZero(a, b))
        continue; //?
      auto el_ab = tab.get(a, b);
      assert(el_ab != nullptr);
      if (el_ab->first == 0.0) {
        *el_ab = srn(h, a, b, omega, dV);
        const auto s = h->symm_sign(a, b);
        *tab.get(b, a) = {s * el_ab->first, s * el_ab->second};
      }
    }
  }

  return tab;
}

//==============================================================================
//==============================================================================
double StructureRad::t1234(int k, const DiracSpinor &w, const DiracSpinor &r,
                           const DiracSpinor &v, const DiracSpinor &c) const {
  return t1(k, w, r, v, c) + t2(k, w, r, v, c) + t3(k, w, r, v, c) +
         t4(k, w, r, v, c);
}

double StructureRad::b1234(int k, const DiracSpinor &w, const DiracSpinor &c,
                           const DiracSpinor &v, const DiracSpinor &r) const {
  const auto st = Angular::neg1pow_2(v.twoj() - w.twoj() + c.twoj() - r.twoj());
  return st * t1234(k, v, r, w, c);
}

//==============================================================================
double StructureRad::t1(const int k, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(v.twoj() + r.twoj() + 2 * k);

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto sab = Angular::neg1pow_2(b.twoj() + a.twoj());
      const auto inv_e_rwab = 1.0 / (r.en() + w.en() - a.en() - b.en());

      const auto [minU, maxU] = Coulomb::k_minmax_W(w, r, b, a);
      const auto [minL, maxL] = Coulomb::k_minmax_Q(v, a, b, c);
      if (minU > maxU)
        continue;
      for (int l = minL; l <= maxL; l += 2) {

        const auto ql = Q(l, v, a, b, c);
        if (Angular::zeroQ(ql))
          continue;

        for (int u = minU; u <= maxU; ++u) {

          // const auto sj1 = mY.SixJ().get_2(w, v, k, l, u, b);
          const auto sj1 = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, 2 * l,
                                           2 * u, b.twoj());
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.SixJ().get_2(r.twoj(), c.twoj(), 2 * k, 2 * l,
                                           2 * u, a.twoj());
          if (Angular::zeroQ(sj2))
            continue;

          const auto wu = W(u, w, r, b, a); // *

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

//==============================================================================
double StructureRad::t2(const int k, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - c.twoj() + 2 * k);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_rnav = 1.0 / (r.en() + n.en() - a.en() - v.en());

      const auto [minU1, maxU1] = Coulomb::k_minmax_W(w, a, c, n);
      const auto [minU2, maxU2] = Coulomb::k_minmax_W(v, a, r, n);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, r.twoj(),
                                        c.twoj(), 2 * u);
        if (Angular::zeroQ(sj))
          continue;

        const auto su = Angular::neg1pow(u);
        const auto f = (2 * u + 1);

        const auto wu1 = W(u, w, a, c, n);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = W(u, v, a, r, n);

        t += su * sj * wu1 * wu2 * inv_e_rnav / f;
      }
    }
  }
  return s * t;
}

//==============================================================================
double StructureRad::t3(const int k, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - c.twoj() + 2 * k);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_nwac = 1.0 / (n.en() + w.en() - a.en() - c.en());

      const auto [minU1, maxU1] = Coulomb::k_minmax_W(w, n, c, a);
      const auto [minU2, maxU2] = Coulomb::k_minmax_W(v, n, r, a);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, r.twoj(),
                                        c.twoj(), 2 * u);
        if (Angular::zeroQ(sj))
          continue;

        const auto su = Angular::neg1pow(u);
        const auto f = (2 * u + 1);

        const auto wu1 = W(u, w, n, c, a);
        if (Angular::zeroQ(wu1))
          continue;

        const auto wu2 = W(u, v, n, r, a);

        t += su * sj * wu1 * wu2 * inv_e_nwac / f;
      }
    }
  }
  return s * t;
}

//==============================================================================
double StructureRad::t4(const int k, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(v.twoj() + r.twoj() + 2 * k);

  for (const auto &n : mExcited) {
    for (const auto &m : mExcited) {

      const auto snm = Angular::neg1pow_2(n.twoj() + m.twoj());
      const auto inv_e_nmcv = 1.0 / (n.en() + m.en() - c.en() - v.en());

      const auto [minU, maxU] = Coulomb::k_minmax_Q(w, r, n, m);
      const auto [minL, maxL] = Coulomb::k_minmax_W(n, m, v, c);
      if (minL > maxL)
        continue;
      for (int u = minU; u <= maxU; u += 2) {

        const auto qu = Q(u, w, r, n, m);
        if (Angular::zeroQ(qu))
          continue;

        for (int l = minL; l <= maxL; ++l) {

          const auto sj1 = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, 2 * l,
                                           2 * u, n.twoj());
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.SixJ().get_2(r.twoj(), c.twoj(), 2 * k, 2 * l,
                                           2 * u, m.twoj());
          if (Angular::zeroQ(sj2))
            continue;

          const auto wl = W(l, v, c, n, m);

          t += snm * sj1 * sj2 * qu * wl * inv_e_nmcv;
        }
      }
    }
  }
  return s * t;
}

//==============================================================================
//==============================================================================
double StructureRad::c1(const int k, const DiracSpinor &w, const DiracSpinor &a,
                        const DiracSpinor &v, const DiracSpinor &c) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - c.twoj() + 2 * k);

  for (const auto &b : mCore) {
    for (const auto &n : mExcited) {

      if (!Angular::sixjTriads(w.twoj(), v.twoj(), 2 * k, c.twoj(), a.twoj(),
                               {}))
        continue;

      const auto e_nvbc = n.en() + v.en() - b.en() - c.en();
      const auto e_nwab = n.en() + w.en() - a.en() - b.en();
      const auto invde = 1.0 / (e_nvbc * e_nwab);

      const auto [minU1, maxU1] = Coulomb::k_minmax_W(w, n, a, b);
      const auto [minU2, maxU2] = Coulomb::k_minmax_W(v, n, c, b);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, c.twoj(),
                                        a.twoj(), 2 * u);
        if (sj == 0.0)
          continue;

        const auto su = Angular::neg1pow(u);
        const auto f = 2.0 * u + 1.0;

        const auto wu1 = W(u, w, n, a, b);
        const auto wu2 = W(u, v, n, c, b);

        t += su * sj * wu1 * wu2 * invde / f;
      }
    }
  }
  return s * t;
}

//==============================================================================
double StructureRad::c2(const int k, const DiracSpinor &w, const DiracSpinor &a,
                        const DiracSpinor &v, const DiracSpinor &c) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(v.twoj() + a.twoj() + 2 * k);

  for (const auto &m : mExcited) {
    for (const auto &n : mExcited) {

      if (!Angular::sixjTriads(w.twoj(), v.twoj(), 2 * k, {}, {}, n.twoj()))
        continue;
      if (!Angular::sixjTriads(c.twoj(), a.twoj(), 2 * k, {}, {}, m.twoj()))
        continue;

      const auto e_mnaw = m.en() + n.en() - a.en() - w.en();
      const auto e_nmcv = n.en() + m.en() - c.en() - v.en();
      const auto invde = 1.0 / (e_mnaw * e_nmcv);
      const auto smn = Angular::neg1pow_2(m.twoj() + n.twoj());

      const auto [minU, maxU] = Coulomb::k_minmax_Q(v, m, n, c);
      const auto [minL, maxL] = Coulomb::k_minmax_W(w, a, n, m);
      if (minL > maxL)
        continue;
      for (int u = minU; u <= maxU; u += 2) {

        const auto s_qu_invde = smn * Q(u, v, m, n, c) * invde;

        for (int l = minL; l <= maxL; ++l) {

          const auto sj1 = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, 2 * u,
                                           2 * l, n.twoj());
          if (sj1 == 0.0)
            continue;

          const auto sj2 = mY.SixJ().get_2(c.twoj(), a.twoj(), 2 * k, 2 * l,
                                           2 * u, m.twoj());
          if (sj2 == 0.0)
            continue;

          const auto wl = W(l, w, a, n, m);

          t += sj1 * sj2 * s_qu_invde * wl;
        }
      }
    }
  }
  return s * t;
}

//==============================================================================
double StructureRad::d1(const int k, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &m) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - m.twoj() + 2 * k);

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      if (!Angular::sixjTriads(w.twoj(), v.twoj(), 2 * k, m.twoj(), r.twoj(),
                               {}))
        continue;

      const auto e_mnav = m.en() + n.en() - a.en() - v.en();
      const auto e_rnaw = r.en() + n.en() - a.en() - w.en();
      const auto invde = 1.0 / (e_mnav * e_rnaw);

      const auto [minU1, maxU1] = Coulomb::k_minmax_W(w, a, r, n);
      const auto [minU2, maxU2] = Coulomb::k_minmax_W(v, a, m, n);
      const auto minU = std::max(minU1, minU2);
      const auto maxU = std::min(maxU1, maxU2);
      for (int u = minU; u <= maxU; ++u) {

        const auto sj = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, m.twoj(),
                                        r.twoj(), 2 * u);
        if (sj == 0.0)
          continue;

        const auto wu1 = W(u, w, a, r, n);
        const auto wu2 = W(u, v, a, m, n);
        const auto fsu = Angular::neg1pow(u) / (2.0 * u + 1.0);

        t += fsu * sj * wu1 * wu2 * invde;
      }
    }
  }
  return s * t;
}

//==============================================================================
double StructureRad::d2(const int k, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &m) const {
  double t = 0.0;

  const auto s = Angular::neg1pow_2(v.twoj() + r.twoj() + 2 * k);

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      if (!Angular::sixjTriads(w.twoj(), v.twoj(), 2 * k, {}, {}, a.twoj()))
        continue;
      if (!Angular::sixjTriads(m.twoj(), r.twoj(), 2 * k, {}, {}, b.twoj()))
        continue;

      const auto sab = Angular::neg1pow_2(a.twoj() + b.twoj());

      const auto e_mvab = m.en() + v.en() - a.en() - b.en();
      const auto e_rwba = r.en() + w.en() - b.en() - a.en();
      const auto invde = 1.0 / (e_mvab * e_rwba);

      const auto [minU, maxU] = Coulomb::k_minmax_Q(v, b, a, m);
      const auto [minL, maxL] = Coulomb::k_minmax_W(w, r, a, b);
      if (minL > maxL)
        continue;
      for (int u = minU; u <= maxU; u += 2) {

        const auto s_qu_invde = sab * Q(u, v, b, a, m) * invde;

        for (int l = minL; l <= maxL; ++l) {

          const auto sj1 = mY.SixJ().get_2(w.twoj(), v.twoj(), 2 * k, 2 * u,
                                           2 * l, a.twoj());
          if (sj1 == 0.0)
            continue;
          const auto sj2 = mY.SixJ().get_2(m.twoj(), r.twoj(), 2 * k, 2 * l,
                                           2 * u, b.twoj());
          if (sj2 == 0.0)
            continue;

          const auto wl = W(l, w, r, a, b);

          t += sj1 * sj2 * s_qu_invde * wl;
        }
      }
    }
  }
  return s * t;
}

//==============================================================================
double StructureRad::dSigma_dE(const DiracSpinor &v, const DiracSpinor &i,
                               const DiracSpinor &j,
                               const DiracSpinor &k) const {
  double t = 0.0;

  const auto e_vijk = (v.en() + i.en() - j.en() - k.en());
  const auto ide2 = 1.0 / (e_vijk * e_vijk);
  const auto tjvp1 = v.twoj() + 1;

  const auto [minU, maxU] = Coulomb::k_minmax_W(v, i, j, k);
  for (int u = minU; u <= maxU; ++u) {

    const auto tup1 = (2 * u + 1);

    const auto q = Q(u, v, i, j, k);
    if (Angular::zeroQ(q))
      continue;

    const auto w = q + P(u, v, i, j, k);

    t += q * w * ide2 / tup1;
  }
  return t / tjvp1;
}
//==============================================================================
double StructureRad::n1(const DiracSpinor &v) const {
  double t = 0.0;
  for (const auto &a : mCore) {
    for (const auto &b : mCore) {
      for (const auto &n : mExcited) {
        t += dSigma_dE(v, n, a, b);
      }
    }
  }
  return t;
}

//==============================================================================
double StructureRad::n2(const DiracSpinor &v) const {
  double t = 0.0;
  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {
      for (const auto &m : mExcited) {
        t += dSigma_dE(v, a, m, n);
      }
    }
  }
  return t;
}

//==============================================================================
std::pair<double, double>
StructureRad::z_bo(const DiracOperator::TensorOperator *const h,
                   const DiracSpinor &w, const DiracSpinor &v, bool transpose,
                   const ExternalField::CorePolarisation *const dV) const {
  //

  double zBO = 0.0, zBO_dv = 0.0;

  std::vector<const DiracSpinor *> i_ptr;
  for (const auto &t_basis : {&mCore, &mExcited}) {
    for (const auto &i : *t_basis) {
      if (i.kappa() != v.kappa() || i == v)
        continue;
      i_ptr.push_back(&i);
    }
  }

#pragma omp parallel for reduction(+ : zBO) reduction(+ : zBO_dv)
  for (std::size_t ii = 0; ii < i_ptr.size(); ++ii) {
    const auto &i = *i_ptr[ii];

    const auto ftr = transpose ? h->symm_sign(w, i) : 1.0;
    const auto h_wi = ftr * h->reducedME(w, i);
    const auto dV_wi = ftr * (dV ? dV->dV(w, i) : 0.0);

    const auto de_vi = v.en() - i.en();
    // const auto sign = Angular::neg1pow_2(v.twoj() - i.twoj());

    const auto f_iv = h_wi / v.twojp1() / de_vi;
    const auto dv_iv = (h_wi + dV_wi) / v.twojp1() / de_vi;

    double zBO_i = 0.0;

    for (const auto &m : mExcited) {
      for (const auto &a : mCore) {
        for (const auto &b : mCore) {
          const auto [k0, km] = Coulomb::k_minmax_Q(i, m, a, b);
          for (auto k = k0; k <= km; k += 2) {
            const auto e_vmab = v.en() + m.en() - a.en() - b.en();
            zBO_i += Q(k, i, m, a, b) * W(k, v, m, a, b) / e_vmab / (2 * k + 1);
          }
        }
        for (const auto &n : mExcited) {
          const auto [k0, km] = Coulomb::k_minmax_Q(i, a, n, m);
          for (auto k = k0; k <= km; k += 2) {
            const auto e_nmav = v.en() + a.en() - n.en() - m.en();
            zBO_i += Q(k, i, a, n, m) * W(k, v, a, n, m) / e_nmav / (2 * k + 1);
          }
        }
      }
    }

    zBO += zBO_i * f_iv;
    zBO_dv += zBO_i * dv_iv;
  }

  return {zBO, zBO_dv};
}

//==============================================================================
std::pair<double, double>
StructureRad::BO(const DiracOperator::TensorOperator *const h,
                 const DiracSpinor &w, const DiracSpinor &v,
                 const ExternalField::CorePolarisation *const dV, double fw,
                 double fv) const {

  auto z1 = z_bo(h, w, v, false, dV);
  const auto z2 = v == w ? z1 : z_bo(h, v, w, true, dV);

  z1.first *= fw;
  z1.second *= fw;

  z1.first += fv * z2.first;
  z1.second += fv * z2.second;
  return z1;
}

} // namespace MBPT
