#include "DiagramRPA.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <string>
#include <vector>

namespace MBPT {

//******************************************************************************
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const DiagramRPA *const drpa)
    : m_k(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ()) {
  //
  if (m_k != drpa->m_k || m_pi != drpa->m_pi) {
    std::cerr << "\nFAIL21 in DiagramRPA: Cannot use 'eat' constructor for "
                 "different rank/parity operators!\n";
    std::abort();
  }

  // Set up basis:
  holes = drpa->holes;
  excited = drpa->excited;

  setup_ts(h);

  Wanmb = drpa->Wanmb;
  Wabmn = drpa->Wabmn;
  Wmnab = drpa->Wmnab;
  Wmban = drpa->Wmban;
}

//******************************************************************************
void DiagramRPA::setup_ts(const DiracOperator::TensorOperator *const h) {
  // Calc t0 (and setup t)
  for (const auto &Fa : holes) {
    std::vector<double> t0a_m;
    for (const auto &Fm : excited) {
      t0a_m.push_back(h->reducedME(Fa, Fm));
    }
    t0am.push_back(t0a_m);
  }
  for (const auto &Fm : excited) {
    std::vector<double> t0m_a;
    for (const auto &Fa : holes) {
      t0m_a.push_back(h->reducedME(Fm, Fa));
    }
    t0ma.push_back(t0m_a);
  }
  clear_tam();
}

//******************************************************************************
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const std::vector<DiracSpinor> &basis,
                       const std::vector<DiracSpinor> &core)
    : m_k(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ()) {

  // Set up basis:
  for (const auto &Fi : basis) {
    const bool inCore = std::find(cbegin(core), cend(core), Fi) != cend(core);
    if (inCore) {
      holes.push_back(Fi);
    } else {
      excited.push_back(Fi);
    }
  }

  // Calc t0 (and setup t) [RPA MEs for hole-excited]
  setup_ts(h);

  // RPA: store W Coulomb integrals (used only for Core RPA its)
  std::cout << "Filling RPA Diagram matrix .. " << std::flush;
  Wanmb.resize(holes.size());
  Wabmn.resize(holes.size());
#pragma omp parallel for
  for (std::size_t i = 0; i < holes.size(); i++) {
    const auto &Fa = holes[i];
    std::vector<std::vector<std::vector<double>>> Wa_nmb;
    std::vector<std::vector<std::vector<double>>> Wa_bmn;
    Wa_nmb.reserve(excited.size());
    Wa_bmn.reserve(excited.size());
    for (const auto &Fn : excited) {
      std::vector<std::vector<double>> Wan_mb;
      std::vector<std::vector<double>> Wab_mn;
      Wan_mb.reserve(excited.size());
      Wab_mn.reserve(excited.size());
      for (const auto &Fm : excited) {
        std::vector<double> Wanm_b;
        std::vector<double> Wabm_n;
        Wanm_b.reserve(holes.size());
        Wabm_n.reserve(holes.size());
        for (const auto &Fb : holes) {
          // XXX Use P,Q,Y
          const bool zero = h->isZero(Fb.k, Fn.k);
          auto x = zero ? 0.0 : Coulomb::Wk_abcd(Fa, Fn, Fm, Fb, m_k);
          auto y = zero ? 0.0 : Coulomb::Wk_abcd(Fa, Fb, Fm, Fn, m_k);
          Wanm_b.push_back(x);
          Wabm_n.push_back(y);
        }
        Wan_mb.push_back(Wanm_b);
        Wab_mn.push_back(Wabm_n);
      }
      Wa_nmb.push_back(Wan_mb);
      Wa_bmn.push_back(Wab_mn);
    }
    Wanmb[i] = Wa_nmb;
    Wabmn[i] = Wa_bmn;
  }
  Wmnab.resize(excited.size());
  Wmban.resize(excited.size());
  std::cout << "." << std::flush;
#pragma omp parallel for
  for (std::size_t i = 0; i < excited.size(); i++) {
    const auto &Fm = excited[i];
    std::vector<std::vector<std::vector<double>>> Wa_nmb;
    std::vector<std::vector<std::vector<double>>> Wa_bmn;
    for (const auto &Fn : excited) {
      std::vector<std::vector<double>> Wan_mb;
      std::vector<std::vector<double>> Wab_mn;
      for (const auto &Fa : holes) {
        std::vector<double> Wanm_b;
        std::vector<double> Wabm_n;
        for (const auto &Fb : holes) {
          // XXX Use P,Q,Y
          const bool zero = h->isZero(Fb.k, Fn.k);
          auto x = zero ? 0.0 : Coulomb::Wk_abcd(Fm, Fn, Fa, Fb, m_k);
          auto y = zero ? 0.0 : Coulomb::Wk_abcd(Fm, Fb, Fa, Fn, m_k);
          Wanm_b.push_back(x);
          Wabm_n.push_back(y);
        }
        Wan_mb.push_back(Wanm_b);
        Wab_mn.push_back(Wabm_n);
      }
      Wa_nmb.push_back(Wan_mb);
      Wa_bmn.push_back(Wab_mn);
    }
    Wmnab[i] = Wa_nmb;
    Wmban[i] = Wa_bmn;
  }
  std::cout << " done.\n" << std::flush;
}

//******************************************************************************
void DiagramRPA::clear_tam() {
  tam = t0am;
  tma = t0ma;
}

//******************************************************************************
double DiagramRPA::dV(const DiracSpinor &Fw, const DiracSpinor &Fv,
                      const bool first_order) const {

  const auto orderOK = Fv.en <= Fw.en;
  const auto &Fi = orderOK ? Fv : Fw;
  const auto &Ff = orderOK ? Fw : Fv;
  const auto f = (1.0 / (2 * m_k + 1));

  double sum = 0.0;
  for (std::size_t ia = 0; ia < holes.size(); ia++) {
    const auto &Fa = holes[ia];
    const auto s1 = ((Fa.twoj() - Ff.twoj() + 2 * m_k) % 4 == 0) ? 1 : -1;
    double sum_a = 0.0;
    for (std::size_t im = 0; im < excited.size(); im++) {
      const auto &Fm = excited[im];
      if (t0am[ia][im] == 0.0)
        continue;
      const auto s2 = ((Fa.twoj() - Fm.twoj()) % 4 == 0) ? 1 : -1;
      const auto Wwmva = Coulomb::Wk_abcd(Ff, Fm, Fi, Fa, m_k);
      const auto Wwavm = Coulomb::Wk_abcd(Ff, Fa, Fi, Fm, m_k);
      const auto ttam = first_order ? t0am[ia][im] : tam[ia][im];
      const auto ttma = first_order ? t0ma[im][ia] : tma[im][ia];
      const auto A = ttam * Wwmva / (Fa.en - Fm.en - m_omega);
      const auto B = Wwavm * ttma / (Fa.en - Fm.en + m_omega);
      sum_a += s1 * (A + s2 * B);
    }
    sum += sum_a;
  }

  return f * sum;
}

//******************************************************************************
void DiagramRPA::rpa_core(const double omega, const bool print) {

  m_omega = std::abs(omega);

  if (print) {
    printf("RPA(D) (w=%.3f): .. \r", m_omega);
    std::cout << std::flush;
  }
  int it = 1;
  auto eps = 0.0;
  const auto f = (1.0 / (2 * m_k + 1));
  for (; it <= max_its; it++) {
    eps = 0.0;
#pragma omp parallel for
    for (std::size_t ia = 0; ia < holes.size(); ia++) {
      const auto &Fa = holes[ia];
      for (std::size_t im = 0; im < excited.size(); im++) {
        const auto &Fm = excited[im];

        double sum_am = 0;
        double sum_ma = 0;

        for (std::size_t ib = 0; ib < holes.size(); ib++) {
          const auto &Fb = holes[ib];
          const auto s1 = ((Fb.twoj() - Fa.twoj() + 2 * m_k) % 4 == 0) ? 1 : -1;
          const auto s3 = ((Fb.twoj() - Fm.twoj() + 2 * m_k) % 4 == 0) ? 1 : -1;
          for (std::size_t in = 0; in < excited.size(); in++) {
            const auto &Fn = excited[in];

            const auto tbn = tam[ib][in];
            if (tbn == 0.0)
              continue;
            const auto tnb = tma[in][ib];
            const auto dem = Fb.en - Fn.en - m_omega;
            const auto dep = Fb.en - Fn.en + m_omega;
            const auto A = tbn * Wanmb[ia][in][im][ib] / dem;
            const auto B = tnb * Wabmn[ia][in][im][ib] / dep;
            const auto C = tbn * Wmnab[im][in][ia][ib] / dem;
            const auto D = tnb * Wmban[im][in][ia][ib] / dep;
            const auto s2 = ((Fb.twoj() - Fn.twoj()) % 4 == 0) ? 1 : -1;
            sum_am += s1 * (A + s2 * B);
            sum_ma += s3 * (C + s2 * D);
          }
        }
        const auto prev = tam[ia][im];
        tam[ia][im] = 0.5 * (tam[ia][im] + t0am[ia][im] + f * sum_am);
        tma[im][ia] = 0.5 * (tma[im][ia] + t0ma[im][ia] + f * sum_ma);
        const auto delta = std::abs((tam[ia][im] - prev) / tam[ia][im]);
#pragma omp critical(compare_eps)
        {
          if (delta > eps) {
            eps = std::abs(delta);
          }
        }
      }
    }
    if (eps < eps_targ)
      break;
    if (print && it % 15 == 0) {
      printf("RPA(D) (w=%.3f): %2i %.1e \r", m_omega, it, eps);
      std::cout << std::flush;
    }
  }
  if (print) {
    printf("RPA(D) (w=%.3f): %2i %.1e\n", m_omega, it, eps);
  }
  m_core_eps = eps;
}

} // namespace MBPT
