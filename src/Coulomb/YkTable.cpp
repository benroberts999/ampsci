#include "YkTable.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace Coulomb {

//******************************************************************************
void YkTable::calculate(const std::vector<DiracSpinor> &a_orbs,
                        const std::vector<DiracSpinor> &b_orbs) {
  // note: make use of the symmetries: force a<=b
  // Note: we RE-calculate all integrals, and calculate any new ones.
  // do not use this to 'extend' the Yab table

  const auto max_2j =
      std::max(DiracSpinor::max_tj(a_orbs), DiracSpinor::max_tj(b_orbs));

  m_Ck.fill(max_2j); // XXX DiragramRPA test fail without +1???
  m_6j.fill(max_2j);

  allocate_space(a_orbs, b_orbs);

  const auto a_is_b = (&a_orbs == &b_orbs);

// for (const auto &a : a_orbs) {
#pragma omp parallel for
  for (auto ia = 0ul; ia < a_orbs.size(); ++ia) {
    const auto &a = a_orbs[ia];
    for (const auto &b : b_orbs) {
      if (a_is_b && b > a)
        continue;
      const auto [k0, kI] = k_minmax(a, b);
      for (auto k = k0; k <= kI; k += 2) {
        // auto &ykab = get_or_insert(std::size_t(k), ab_key(a, b));
        auto &ykab = get_ref(k, a, b);
        Coulomb::yk_ab(a, b, k, ykab);
      }
    }
  }
}

//******************************************************************************
void YkTable::extend_angular(int new_max_2j) {
  m_Ck.fill(new_max_2j);
  m_6j.fill(new_max_2j);
}

//******************************************************************************
void YkTable::allocate_space(const std::vector<DiracSpinor> &a_orbs,
                             const std::vector<DiracSpinor> &b_orbs) {

  const auto a_is_b = (&a_orbs == &b_orbs);

  for (const auto &a : a_orbs) {
    for (const auto &b : b_orbs) {
      if (a_is_b && b > a)
        continue;
      const auto [k0, kI] = k_minmax(a, b);
      for (auto k = k0; k <= kI; k += 2) {
        get_or_insert(std::size_t(k), ab_key(a, b));
      }
    }
  }
}

//******************************************************************************
const std::vector<double> *YkTable::get(const int k, const DiracSpinor &Fa,
                                        const DiracSpinor &Fb) const {
  const auto sk = static_cast<std::size_t>(k);
  if (sk >= m_Y.size())
    return nullptr;
  // const auto key = ;
  const auto it = m_Y[sk].find(ab_key(Fa, Fb));
  if (it == m_Y[sk].cend()) {
    return nullptr;
  }
  return &(it->second);
}

//******************************************************************************
uint32_t YkTable::ab_key(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
  // symmetric: define F1=min(Fa,Fb), F2=max(Fa,Fb)
  static_assert(sizeof(uint32_t) == 2 * sizeof(uint16_t), "32=2*16");
  const std::array<uint16_t, 2> key_array{
      static_cast<uint16_t>(std::min(Fa.nk_index(), Fb.nk_index())),
      static_cast<uint16_t>(std::max(Fa.nk_index(), Fb.nk_index())) //
  };
  uint32_t tkey;
  std::memcpy(&tkey, &key_array, sizeof(uint32_t));
  return tkey;
}

//******************************************************************************
std::vector<double> &YkTable::get_or_insert(std::size_t k, uint32_t key) {
  if (k >= m_Y.size()) {
    m_Y.resize(k + 1);
  }
  // Returns reference to vector at 'key'; if no such vector exists, creates
  // it first.
  const auto [new_it, ok] = m_Y[k].emplace(key, std::vector<double>{});
  return new_it->second;
}

//******************************************************************************
std::vector<double> &YkTable::get_ref(const int k, const DiracSpinor &Fa,
                                      const DiracSpinor &Fb) {
  // May only call this function is assured vector exists
  const auto sk = static_cast<std::size_t>(k);
  assert(sk < m_Y.size());
  const auto it = m_Y[sk].find(ab_key(Fa, Fb));
  assert(it != m_Y[sk].cend());
  return it->second;
}

//****************************************************************************
double YkTable::Q(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                  const DiracSpinor &Fc, const DiracSpinor &Fd) const {
  // // nb: b and d _MUST_ be in {a},{b} orbitals
  // // assert(m_aisb && "May only use Qk if Yk init with {a}={b}");
  // const auto ykbd = get(k, Fb, Fd);
  // return ykbd ? Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k, *ykbd, m_Ck) : 0.0;

  const auto tCac = m_Ck.get_tildeCkab(k, Fa.k, Fc.k);
  if (Angular::zeroQ(tCac))
    return 0.0;
  const auto tCbd = m_Ck.get_tildeCkab(k, Fb.k, Fd.k);
  if (Angular::zeroQ(tCbd))
    return 0.0;
  const auto ykbd = get(k, Fb, Fd);
  const auto Rkabcd = Coulomb::Rk_abcd(Fa, Fc, *ykbd);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return m1tk * tCac * tCbd * Rkabcd;
}

//******************************************************************************
double YkTable::P(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                  const DiracSpinor &Fc, const DiracSpinor &Fd) const {
  // Pk = [k] Sum_l {a,c,k,b,d,l} Q^l_abdc

  if (Angular::triangle(Fa.twoj(), Fc.twoj(), 2 * k) == 0)
    return 0.0;
  if (Angular::triangle(2 * k, Fb.twoj(), Fd.twoj()) == 0)
    return 0.0;

  double pk = 0.0;
  const auto [l0, lI] = Coulomb::k_minmax_Q(Fa, Fb, Fd, Fc);
  for (int l = l0; l <= lI; l += 2) {
    const auto ylbc = get(l, Fb, Fc);
    assert(ylbc != nullptr);

    const auto Ql = Q(l, Fa, Fb, Fd, Fc);
    const auto sj = m_6j.get(Fa, Fc, k, Fb, Fd, l);

    pk += Ql * sj;
  }
  return pk * (2 * k + 1);
}

//******************************************************************************
double YkTable::W(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                  const DiracSpinor &Fc, const DiracSpinor &Fd) const {
  return Q(k, Fa, Fb, Fc, Fd) + P(k, Fa, Fb, Fc, Fd);
}

//******************************************************************************
DiracSpinor YkTable::Qkv_bcd(int kappa, const DiracSpinor &Fb,
                             const DiracSpinor &Fc, const DiracSpinor &Fd,
                             const int k) const {
  DiracSpinor Qkv{0, kappa, Fb.rgrid};
  const auto tCac = m_Ck.get_tildeCkab(k, kappa, Fc.k);
  const auto tCbd = m_Ck.get_tildeCkab(k, Fb.k, Fd.k);
  const auto tCC = tCbd * tCac;
  if (tCC == 0.0) {
    // Qkv.scale(0.0);
    return Qkv;
  }
  // const auto ylbc = get(l, Fb, Fc);
  const auto ykbd = get(k, Fb, Fd);
  Coulomb::Rkv_bcd(&Qkv, Fc, *ykbd);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  Qkv.scale(m1tk * tCC);
  return Qkv;
}

//******************************************************************************
DiracSpinor YkTable::Pkv_bcd(int kappa, const DiracSpinor &Fb,
                             const DiracSpinor &Fc, const DiracSpinor &Fd,
                             const int k,
                             const std::vector<double> &f2k) const {

  DiracSpinor Pkv{0, kappa, Fb.rgrid};

  const auto fk = [&f2k](int l) {
    // nb: only screens l, k assumed done outside...
    if (l < int(f2k.size())) {
      return f2k[std::size_t(l)];
    }
    return 1.0;
  };

  const auto tkp1 = 2 * k + 1;
  const auto [l0, lI] = Coulomb::k_minmax_Q(Pkv, Fb, Fd, Fc);
  for (int l = l0; l <= lI; l += 2) {
    const auto ylbc = get(l, Fb, Fc);
    assert(ylbc != nullptr);

    const auto sj = fk(l) * m_6j.get_2(Fc.twoj(), Angular::twoj_k(kappa), 2 * k,
                                       Fd.twoj(), Fb.twoj(), 2 * l);

    if (Angular::zeroQ(sj))
      continue;
    Pkv += sj * Qkv_bcd(Pkv.k, Fb, Fd, Fc, l);
  }
  Pkv *= tkp1;
  return Pkv;
}

} // namespace Coulomb
