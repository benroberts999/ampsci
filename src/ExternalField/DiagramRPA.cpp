#include "DiagramRPA.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

namespace ExternalField {

//==============================================================================
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const std::vector<DiracSpinor> &basis,
                       const HF::HartreeFock *in_hf, const std::string &atom,
                       bool print)
  : CorePolarisation(h), p_hf(in_hf) {

  assert(p_hf != nullptr && "Hartree-Fock cannot be null in DiagramRPA");
  assert(h != nullptr && "Operator cannot be null in DiagramRPA");

  // Set up basis (split core-excited):
  const int n_min_core = 1;
  auto [t_holes, t_excited] =
    DiracSpinor::split_by_core(basis, p_hf->core(), n_min_core);
  m_holes = std::move(t_holes);
  m_excited = std::move(t_excited);

  // Setup faster Breit
  if (p_hf->vBreit() != nullptr) {
    m_Br = *p_hf->vBreit();

    // nb: This [fill_gb()] uses HUGE amount of memory, leads to ~2x speedup
    // Decided not worth the speedup
    // m_Br->fill_gb(basis);
  }

  // filename for W^k (Coulmb integrals) output
  const auto basis_string = DiracSpinor::state_config(basis);
  const auto fname = atom + "_" + std::to_string(m_rank) +
                     (m_pi == 1 ? "+" : "-") + "_" + basis_string + ".rpad.abf";
  const auto do_read_write = atom != "" && atom != "false";

  // Attempt to read W's from a file:
  const auto read_ok = do_read_write ? read_write(fname, IO::FRW::read) : false;
  if (!read_ok) {
    // If not, calc W's, and write to file
    fill_W_matrix(h, print);
    if (!m_holes.empty() && !m_excited.empty() && do_read_write)
      read_write(fname, IO::FRW::write);
  }
}

//==============================================================================
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const DiagramRPA *const drpa)
  : CorePolarisation(h), p_hf(drpa->p_hf) {

  assert(m_rank == drpa->m_rank &&
         "If constructing a DiagramRPA from another, must have same rank");
  assert(m_pi == drpa->m_pi &&
         "If constructing a DiagramRPA from another, must have same parity");

  // Set up basis:
  m_holes = drpa->m_holes;
  m_excited = drpa->m_excited;

  // "eat" W matrices from other rpa
  m_Wanmb = drpa->m_Wanmb;
  m_Wabmn = drpa->m_Wabmn;
  m_Wmnab = drpa->m_Wmnab;
  m_Wmban = drpa->m_Wmban;

  update_t0s(h);
}

//==============================================================================
bool DiagramRPA::read_write(const std::string &fname, IO::FRW::RoW rw) {
  // Note: only writes W (depends on k/pi, and basis). Do not write t's, since
  // they depend on operator. This makes it very fast when making small changes
  // to operator (don't need to re-calc W)

  const auto readQ = rw == IO::FRW::read;

  if (readQ && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str = !readQ ? "Writing to " : "Reading from ";
  std::cout << rw_str << "RPA(diagram) file: " << fname << " ("
            << DiracSpinor::state_config(m_holes) << "/"
            << DiracSpinor::state_config(m_excited) << ") ... " << std::flush;

  if (readQ)
    std::cout
      << "\nNote: still uses Basis for summation (only reads in W matrix)\n";

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  if (m_holes.empty() || m_excited.empty()) {
    return false;
  }

  // Note: Basis states must match exactly (since use their index across arrays)
  // Check if same. If not, print status and calc W from scratch

  std::size_t hs = m_holes.size(), es = m_excited.size();
  rw_binary(iofs, rw, hs, es);
  if (readQ) {
    if (hs != m_holes.size() || es != m_excited.size()) {
      std::cout << "\nCannot read from " << fname << ". Basis mis-match (read "
                << hs << "," << es << "; expected " << m_holes.size() << ","
                << m_excited.size() << ").\n"
                << "Will recalculate rpa_Diagram matrix, and overwrite file.\n";
      return false;
    }
  }

  for (const auto porbs : {&m_holes, &m_excited}) {
    for (const auto &Fn : *porbs) {
      int n = Fn.n();
      int k = Fn.kappa();
      rw_binary(iofs, rw, n, k);
      if (readQ) {
        if (Fn.n() != n || Fn.kappa() != k) {
          std::cout
            << "\nCannot read from " << fname << ". Basis mis-match (read " << n
            << "," << k << "; expected " << Fn.n() << "," << Fn.kappa()
            << ").\n"
            << "Will recalculate rpa_Diagram matrix, and overwrite file.\n";
          return false;
        }
      }
    }
  }

  // read/write Ws:
  rw_binary(iofs, rw, m_Wanmb, m_Wabmn, m_Wmnab, m_Wmban);
  std::cout << "done.\n";

  return true;
}

//==============================================================================
void DiagramRPA::fill_W_matrix(const DiracOperator::TensorOperator *const h,
                               bool print) {
  if (m_holes.empty() || m_excited.empty()) {
    std::cout << "\nWARNING 64 in DiagramRPA: no basis! RPA will be zero\n";
    return;
  }

  Coulomb::YkTable Yhe(m_holes, m_excited);
  Yhe.calculate(m_excited);
  Yhe.calculate(m_holes);

  // const auto Vbr = p_hf->vBreit(); // a pointer, may be null

  // RPA: store W Coulomb integrals (used only for Core RPA its)
  if (print)
    std::cout << "Filling RPA Diagram matrix ("
              << DiracSpinor::state_config(m_holes) << "/"
              << DiracSpinor::state_config(m_excited) << ") .. " << std::flush;
  m_Wanmb.resize(m_holes.size());
  m_Wabmn.resize(m_holes.size());
  // First set: only use Yhe and Yee
  {
    const auto &Yee = Yhe;
#pragma omp parallel for
    for (std::size_t i = 0; i < m_holes.size(); i++) {
      const auto &Fa = m_holes[i];
      auto &Wa_nmb = m_Wanmb[i];
      auto &Wa_bmn = m_Wabmn[i];
      Wa_nmb.reserve(m_excited.size());
      Wa_bmn.reserve(m_excited.size());
      for (const auto &Fn : m_excited) {
        auto &Wan_mb = Wa_nmb.emplace_back();
        auto &Wab_mn = Wa_bmn.emplace_back();
        Wan_mb.reserve(m_excited.size());
        Wab_mn.reserve(m_excited.size());
        for (const auto &Fm : m_excited) {
          auto &Wanm_b = Wan_mb.emplace_back();
          auto &Wabm_n = Wab_mn.emplace_back();
          Wanm_b.reserve(m_holes.size());
          Wabm_n.reserve(m_holes.size());
          for (const auto &Fb : m_holes) {
            if (h->isZero(Fb.kappa(), Fn.kappa())) {
              Wanm_b.emplace_back(0.0);
              Wabm_n.emplace_back(0.0);
              continue;
            }
            const auto xQ = Yhe.Q(m_rank, Fa, Fn, Fm, Fb);
            const auto xP = Yee.P(m_rank, Fa, Fn, Fm, Fb);
            const auto yQ = Yee.Q(m_rank, Fa, Fb, Fm, Fn);
            const auto yP = Yhe.P(m_rank, Fa, Fb, Fm, Fn);
            // Breit Contribution to core:
            const auto xB = m_Br ? m_Br->BWk_abcd(m_rank, Fa, Fn, Fm, Fb) : 0.0;
            const auto yB = m_Br ? m_Br->BWk_abcd(m_rank, Fa, Fb, Fm, Fn) : 0.0;

            Wanm_b.push_back(xQ + xP + xB);
            Wabm_n.push_back(yQ + yP + yB);
          }
        }
      }
    }
  }

  m_Wmnab.resize(m_excited.size());
  m_Wmban.resize(m_excited.size());
  if (print)
    std::cout << "." << std::flush;
  {
    // Only use Yhe and Yhh here:
    // const Coulomb::YkTable Yhh(holes);
    const auto &Yhh = Yhe;
#pragma omp parallel for
    for (std::size_t i = 0; i < m_excited.size(); i++) {
      const auto &Fm = m_excited[i];
      auto &Wa_nmb = m_Wmnab[i];
      auto &Wa_bmn = m_Wmban[i];
      Wa_nmb.reserve(m_excited.size());
      Wa_bmn.reserve(m_excited.size());
      for (const auto &Fn : m_excited) {
        auto &Wan_mb = Wa_nmb.emplace_back();
        auto &Wab_mn = Wa_bmn.emplace_back();
        Wan_mb.reserve(m_holes.size());
        Wab_mn.reserve(m_holes.size());
        for (const auto &Fa : m_holes) {
          auto &Wanm_b = Wan_mb.emplace_back();
          auto &Wabm_n = Wab_mn.emplace_back();
          Wanm_b.reserve(m_holes.size());
          Wabm_n.reserve(m_holes.size());
          for (const auto &Fb : m_holes) {
            if (h->isZero(Fb.kappa(), Fn.kappa())) {
              // do I need to store the zero's? Or can I skip them below?
              // skipping is dangerous, mis-match of indexes
              // However, we roughly store 2x what we need!
              Wanm_b.emplace_back(0.0);
              Wabm_n.emplace_back(0.0);
              continue;
            }
            const auto xQ = Yhe.Q(m_rank, Fm, Fn, Fa, Fb);
            const auto xP = Yhe.P(m_rank, Fm, Fn, Fa, Fb);
            const auto yQ = Yhe.Q(m_rank, Fm, Fb, Fa, Fn);
            const auto yP = Yhh.P(m_rank, Fm, Fb, Fa, Fn);
            // nb: Breit part not double-checked! XXX
            const auto xB = m_Br ? m_Br->BWk_abcd(m_rank, Fm, Fn, Fa, Fb) : 0.0;
            const auto yB = m_Br ? m_Br->BWk_abcd(m_rank, Fm, Fb, Fa, Fn) : 0.0;
            Wanm_b.push_back(xQ + xP + xB);
            Wabm_n.push_back(yQ + yP + yB);
          }
        }
      }
    }
  }
  if (print)
    std::cout << " done.\n" << std::flush;
}

//==============================================================================
void DiagramRPA::setup_ts(const DiracOperator::TensorOperator *const h) {
  if (m_holes.empty() || m_excited.empty())
    return;

  m_t0am.clear();
  m_t0ma.clear();

  m_t0am.reserve(m_holes.size());
  m_t0ma.reserve(m_excited.size());
  // Calc t0 (and setup t)
  for (const auto &Fa : m_holes) {
    std::vector<double> t0a_m;
    t0a_m.reserve(m_excited.size());
    for (const auto &Fm : m_excited) {
      t0a_m.push_back(h->reducedME(Fa, Fm));
    }
    m_t0am.push_back(t0a_m);
  }
  for (const auto &Fm : m_excited) {
    std::vector<double> t0m_a;
    t0m_a.reserve(m_holes.size());
    for (const auto &Fa : m_holes) {
      t0m_a.push_back(h->reducedME(Fm, Fa));
    }
    m_t0ma.push_back(t0m_a);
  }
  m_tam = m_t0am;
  m_tma = m_t0ma;
}
//==============================================================================
void DiagramRPA::clear() {
  m_tam = m_t0am;
  m_tma = m_t0ma;
}

//==============================================================================
void DiagramRPA::update_t0s(const DiracOperator::TensorOperator *const h) {
  if (h != nullptr) {
    assert(h->rank() == m_rank && "Rank must match in update_t0s");
    assert(h->parity() == m_pi && "Parity must match in update_t0s");
    assert(h->imaginaryQ() == m_imag && "Imaginarity must match in update_t0s");
    m_h = h;
  }

  // on first call, t0 are empty. Just fill them
  if (m_t0am.empty() || m_tam.empty()) {
    setup_ts(h);
    return;
  }

  assert(m_t0am.size() == m_holes.size());
  assert(m_t0ma.size() == m_excited.size());
  if (m_holes.size() > 0) {
    assert(m_t0am.at(0).size() == m_excited.size());
  }
  if (m_excited.size() > 0) {
    assert(m_t0ma.at(0).size() == m_holes.size());
  }
  for (std::size_t ia = 0; ia < m_holes.size(); ++ia) {
    const auto &Fa = m_holes[ia];
    for (std::size_t im = 0; im < m_excited.size(); ++im) {
      const auto &Fm = m_excited[im];
      m_t0am[ia][im] = m_h->reducedME(Fa, Fm);
      m_t0ma[im][ia] = m_h->symm_sign(Fa, Fm) * m_t0am[ia][im];
    }
  }
  // Don't update Ts (faster convergance)
  // m_tam = m_t0am;
  // m_tma = m_t0ma;
}

//==============================================================================
double DiagramRPA::dV(const DiracSpinor &Fw, const DiracSpinor &Fv) const {

  if (m_holes.empty() || m_excited.empty() || m_tam.empty() || m_t0am.empty())
    return 0.0;

  const auto orderOK = true;
  const auto &Fi = orderOK ? Fv : Fw;
  const auto &Ff = orderOK ? Fw : Fv;
  const auto ww = m_core_omega;

  const auto f = (1.0 / (2 * m_rank + 1));

  std::vector<double> sum_a(m_holes.size());
#pragma omp parallel for
  for (std::size_t ia = 0; ia < m_holes.size(); ia++) {
    const auto &Fa = m_holes[ia];
    const auto s1 = ((Fa.twoj() - Ff.twoj() + 2 * m_rank) % 4 == 0) ? 1 : -1;
    for (std::size_t im = 0; im < m_excited.size(); im++) {
      const auto &Fm = m_excited[im];
      if (m_t0am[ia][im] == 0.0)
        continue;
      const auto s2 = ((Fa.twoj() - Fm.twoj()) % 4 == 0) ? 1 : -1;
      // Calculate Wk from scratch here: Fi/Ff may be valence.
      const auto Wwmva = Coulomb::Wk_abcd(m_rank, Ff, Fm, Fi, Fa) +
                         (m_Br ? m_Br->BWk_abcd(m_rank, Ff, Fm, Fi, Fa) : 0.0);
      const auto Wwavm = Coulomb::Wk_abcd(m_rank, Ff, Fa, Fi, Fm) +
                         (m_Br ? m_Br->BWk_abcd(m_rank, Ff, Fa, Fi, Fm) : 0.0);
      const auto ttam = m_tam[ia][im];
      const auto ttma = m_tma[im][ia];
      const auto A = ttam * Wwmva / (Fa.en() - Fm.en() - ww);
      const auto B = Wwavm * ttma / (Fa.en() - Fm.en() + ww);
      sum_a[ia] += s1 * (A + s2 * B);
    }
  }

  return f * std::accumulate(begin(sum_a), end(sum_a), 0.0);
}

//==============================================================================
void DiagramRPA::solve_core(double omega, int max_its, bool print) {

  const auto eps_targ = m_eps;
  const auto a_damp = m_eta; // ? or always 0.5
  const auto b_damp = 1.0 - a_damp;

  m_core_omega = omega;

  if (m_holes.empty() || m_excited.empty())
    return;

  if (print) {
    fmt::print("RPA(D) {:s} (w={:.4f}): ", m_h->name(), m_core_omega);
    std::cout << std::flush;
  }

  // Start at 2:
  int its_performed = 0;
  auto eps = 1.0;
  std::string s_worst;
  const auto f = (1.0 / (2 * m_rank + 1));

  // Assurs that dV=0 if max_its=0
  // 1 iteration is just setup_ts()
  if (max_its != 0) {
    ++its_performed;
    update_t0s(m_h);
  }

  for (int it = 2; it <= max_its; it++) {
    ++its_performed;

    std::vector<std::pair<double, std::string>> eps_m(m_excited.size());
    // Use these internally - should be values from previous iteration!
    const auto Tnb = m_tma;
    const auto Tbn = m_tam;

#pragma omp parallel for
    for (std::size_t im = 0; im < m_excited.size(); im++) {
      const auto &Fm = m_excited[im];

      double eps_worst_a = 0.0;
      std::string t_worst;
      for (std::size_t ia = 0; ia < m_holes.size(); ia++) {
        const auto &Fa = m_holes[ia];

        double sum_am = 0.0;
        double sum_ma = 0.0;

        for (std::size_t ib = 0; ib < m_holes.size(); ib++) {
          const auto &Fb = m_holes[ib];
          const auto s1 =
            ((Fb.twoj() - Fa.twoj() + 2 * m_rank) % 4 == 0) ? 1 : -1;
          const auto s3 =
            ((Fb.twoj() - Fm.twoj() + 2 * m_rank) % 4 == 0) ? 1 : -1;
          for (std::size_t in = 0; in < m_excited.size(); in++) {
            const auto &Fn = m_excited[in];

            const auto tbn = Tbn[ib][in];
            if (tbn == 0.0)
              continue;
            const auto tnb = Tnb[in][ib];
            const auto tdem = tbn / (Fb.en() - Fn.en() - m_core_omega);
            const auto s2 = ((Fb.twoj() - Fn.twoj()) % 4 == 0) ? 1 : -1;
            const auto stdep = s2 * tnb / (Fb.en() - Fn.en() + m_core_omega);
            {
              const auto A = tdem * m_Wanmb[ia][in][im][ib];
              const auto B = stdep * m_Wabmn[ia][in][im][ib];
              const auto C = tdem * m_Wmnab[im][in][ia][ib];
              const auto D = stdep * m_Wmban[im][in][ia][ib];
              sum_am += s1 * (A + B);
              sum_ma += s3 * (C + D);
            }
          }
        }

        // Update core-excited matrix elements, including damping
        const auto prev = m_tam[ia][im];
        m_tam[ia][im] =
          a_damp * m_tam[ia][im] + b_damp * (m_t0am[ia][im] + f * sum_am);
        m_tma[im][ia] =
          a_damp * m_tma[im][ia] + b_damp * (m_t0ma[im][ia] + f * sum_ma);
        const auto delta =
          2.0 * std::abs((m_tam[ia][im] - prev) / (prev + m_tam[ia][im]));

        if (delta > eps_worst_a) {
          eps_worst_a = delta;
          t_worst = Fa.shortSymbol() + "," + Fm.shortSymbol();
        }
      }
      eps_m[im] = {eps_worst_a, t_worst};
    }

    const auto teps =
      *std::max_element(cbegin(eps_m), cend(eps_m),
                        [](auto &a, auto &b) { return a.first < b.first; });
    eps = teps.first;
    s_worst = teps.second;
    if (eps < eps_targ)
      break;
  }

  if (print) {
    printf("%2i %.1e [%s]\n", its_performed, eps, s_worst.c_str());
    std::cout << std::flush;
  }

  m_core_eps = eps;
  m_core_its = its_performed;
}

} // namespace ExternalField
