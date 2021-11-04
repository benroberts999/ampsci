#include "DiagramRPA.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

namespace ExternalField {

template <typename T, typename U> constexpr T MyCast(U x) {
  // Convert between 'Wtype' and double. Wype may be doule or float.
  // Same as static_cast, but don't warn on useless casts
  // (Use when you know it may be a useless cast)
  if constexpr (std::is_same_v<T, U>) {
    return x;
  } else {
    return static_cast<T>(x);
  }
}

//******************************************************************************
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const std::vector<DiracSpinor> &basis,
                       const std::vector<DiracSpinor> &core,
                       const std::string &atom)
    : CorePolarisation(h) {

  // Set up basis:
  for (const auto &Fi : basis) {
    // NB: this makes a huge difference! Try to implement somehow
    // NB: Need to be careful: same basis used when reading W in!
    // if (max_de > 0.0 && std::abs(Fi.en()) > max_de)
    //   continue;
    const bool inCore = std::find(cbegin(core), cend(core), Fi) != cend(core);
    if (inCore) {
      holes.push_back(Fi);
    } else {
      excited.push_back(Fi);
    }
  }

  // Calc t0 (and setup t) [RPA MEs for hole-excited]
  setup_ts(h);

  const auto fname =
      atom + "_" + std::to_string(m_rank) + (m_pi == 1 ? "+" : "-") + ".rpad";

  // Attempt to read W's from a file:
  const auto read_ok = read_write(fname, IO::FRW::read);
  if (!read_ok) {
    // If not, calc W's, and write to file
    fill_W_matrix(h);
    if (!holes.empty() && !excited.empty() && atom != "")
      read_write(fname, IO::FRW::write);
  }
}

//******************************************************************************
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const DiagramRPA *const drpa)
    : CorePolarisation(h) {

  if (m_rank != drpa->m_rank || m_pi != drpa->m_pi) {
    std::cerr << "\nFAIL21 in DiagramRPA: Cannot use 'eat' constructor for "
                 "different rank/parity operators!\n";
    std::abort();
  }

  // Set up basis:
  holes = drpa->holes;
  excited = drpa->excited;

  setup_ts(h);

  // "eat" W matrices from other rpa
  Wanmb = drpa->Wanmb;
  Wabmn = drpa->Wabmn;
  Wmnab = drpa->Wmnab;
  Wmban = drpa->Wmban;
}

//******************************************************************************
bool DiagramRPA::read_write(const std::string &fname, IO::FRW::RoW rw) {
  // Note: only writes W (depends on k/pi, and basis). Do not write t's, since
  // they depend on operator. This makes it very fast when making small changes
  // to operator (don't need to re-calc W)

  const auto readQ = rw == IO::FRW::read;

  if (readQ && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str = !readQ ? "Writing to " : "Reading from ";
  std::cout << rw_str << "RPA(diagram) file: " << fname << " ("
            << DiracSpinor::state_config(holes) << "/"
            << DiracSpinor::state_config(excited) << ") ... " << std::flush;

  if (readQ)
    std::cout
        << "\nNote: still uses Basis for summation (only reads in W matrix)\n";

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  if (holes.empty() || excited.empty()) {
    return false;
  }

  // Note: Basis states must match exactly (since use their index across arrays)
  // Check if same. If not, print status and calc W from scratch

  std::size_t hs = holes.size(), es = excited.size();
  rw_binary(iofs, rw, hs, es);
  if (readQ) {
    if (hs != holes.size() || es != excited.size()) {
      std::cout << "\nCannot read from " << fname << ". Basis mis-match (read "
                << hs << "," << es << "; expected " << holes.size() << ","
                << excited.size() << ").\n"
                << "Will recalculate rpa_Diagram matrix, and overwrite file.\n";
      return false;
    }
  }

  for (const auto porbs : {&holes, &excited}) {
    for (const auto &Fn : *porbs) {
      int n = Fn.n;
      int k = Fn.k;
      rw_binary(iofs, rw, n, k);
      if (readQ) {
        if (Fn.n != n || Fn.k != k) {
          std::cout
              << "\nCannot read from " << fname << ". Basis mis-match (read "
              << n << "," << k << "; expected " << Fn.n << "," << Fn.k << ").\n"
              << "Will recalculate rpa_Diagram matrix, and overwrite file.\n";
          return false;
        }
      }
    }
  }

  // read/write Ws:
  rw_binary(iofs, rw, Wanmb, Wabmn, Wmnab, Wmban);
  std::cout << "done.\n";

  return true;
}

//******************************************************************************
void DiagramRPA::fill_W_matrix(const DiracOperator::TensorOperator *const h) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  if (holes.empty() || excited.empty()) {
    std::cout << "\nWARNING 64 in DiagramRPA: no basis! RPA will be zero\n";
    return;
  }

  Coulomb::YkTable Yhe(holes, excited);
  Yhe.calculate(excited);
  Yhe.calculate(holes);

  // RPA: store W Coulomb integrals (used only for Core RPA its)
  std::cout << "Filling RPA Diagram matrix ("
            << DiracSpinor::state_config(holes) << "/"
            << DiracSpinor::state_config(excited) << ") .. " << std::flush;
  Wanmb.resize(holes.size());
  Wabmn.resize(holes.size());
  // First set: only use Yhe and Yee
  {
    const auto &Yee = Yhe;
#pragma omp parallel for
    for (std::size_t i = 0; i < holes.size(); i++) {
      const auto &Fa = holes[i];
      auto &Wa_nmb = Wanmb[i];
      auto &Wa_bmn = Wabmn[i];
      Wa_nmb.reserve(excited.size());
      Wa_bmn.reserve(excited.size());
      for (const auto &Fn : excited) {
        auto &Wan_mb = Wa_nmb.emplace_back();
        auto &Wab_mn = Wa_bmn.emplace_back();
        Wan_mb.reserve(excited.size());
        Wab_mn.reserve(excited.size());
        for (const auto &Fm : excited) {
          auto &Wanm_b = Wan_mb.emplace_back();
          auto &Wabm_n = Wab_mn.emplace_back();
          Wanm_b.reserve(holes.size());
          Wabm_n.reserve(holes.size());
          for (const auto &Fb : holes) {
            if (h->isZero(Fb.k, Fn.k)) {
              Wanm_b.emplace_back(0.0);
              Wabm_n.emplace_back(0.0);
              continue;
            }
            const auto xQ = Yhe.Q(m_rank, Fa, Fn, Fm, Fb);
            const auto xP = Yee.P(m_rank, Fa, Fn, Fm, Fb);
            const auto yQ = Yee.Q(m_rank, Fa, Fb, Fm, Fn);
            const auto yP = Yhe.P(m_rank, Fa, Fb, Fm, Fn);
            Wanm_b.push_back(MyCast<Wtype>(xQ + xP));
            Wabm_n.push_back(MyCast<Wtype>(yQ + yP));
          }
        }
      }
    }
  }

  Wmnab.resize(excited.size());
  Wmban.resize(excited.size());
  std::cout << "." << std::flush;
  {
    // Only use Yhe and Yhh here:
    // const Coulomb::YkTable Yhh(holes);
    const auto &Yhh = Yhe;
#pragma omp parallel for
    for (std::size_t i = 0; i < excited.size(); i++) {
      const auto &Fm = excited[i];
      auto &Wa_nmb = Wmnab[i];
      auto &Wa_bmn = Wmban[i];
      Wa_nmb.reserve(excited.size());
      Wa_bmn.reserve(excited.size());
      for (const auto &Fn : excited) {
        auto &Wan_mb = Wa_nmb.emplace_back();
        auto &Wab_mn = Wa_bmn.emplace_back();
        Wan_mb.reserve(holes.size());
        Wab_mn.reserve(holes.size());
        for (const auto &Fa : holes) {
          auto &Wanm_b = Wan_mb.emplace_back();
          auto &Wabm_n = Wab_mn.emplace_back();
          Wanm_b.reserve(holes.size());
          Wabm_n.reserve(holes.size());
          for (const auto &Fb : holes) {
            if (h->isZero(Fb.k, Fn.k)) {
              Wanm_b.emplace_back(0.0);
              Wabm_n.emplace_back(0.0);
              continue;
            }
            const auto xQ = Yhe.Q(m_rank, Fm, Fn, Fa, Fb);
            const auto xP = Yhe.P(m_rank, Fm, Fn, Fa, Fb);
            const auto yQ = Yhe.Q(m_rank, Fm, Fb, Fa, Fn);
            const auto yP = Yhh.P(m_rank, Fm, Fb, Fa, Fn);
            Wanm_b.push_back(MyCast<Wtype>(xQ + xP));
            Wabm_n.push_back(MyCast<Wtype>(yQ + yP));
          }
        }
      }
    }
  }
  std::cout << " done.\n" << std::flush;
}

//******************************************************************************
void DiagramRPA::setup_ts(const DiracOperator::TensorOperator *const h) {
  if (holes.empty() || excited.empty())
    return;

  t0am.clear();
  t0ma.clear();

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
  clear();
}
//******************************************************************************
void DiagramRPA::clear() {
  tam = t0am;
  tma = t0ma;
}

//******************************************************************************
double DiagramRPA::dV(const DiracSpinor &Fw, const DiracSpinor &Fv) const {
  return dV_diagram(Fw, Fv, false);
}

//******************************************************************************
double DiagramRPA::dV_diagram(const DiracSpinor &Fw, const DiracSpinor &Fv,
                              const bool first_order) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  if (holes.empty() || excited.empty())
    return 0.0;

  if (Fv.en() > Fw.en()) {
    // ??????
    const auto sj = ((Fv.twoj() - Fw.twoj()) % 4 == 0) ? 1 : -1;
    const auto si = m_imag ? -1 : 1;
    return (sj * si) * dV_diagram(Fv, Fw, first_order);
  }

  const auto orderOK = true; // Fv.en() <= Fw.en();
  // ??????
  const auto &Fi = orderOK ? Fv : Fw;
  const auto &Ff = orderOK ? Fw : Fv;
  const auto ww = m_core_omega;
  // Fv.en() <= Fw.en() ? m_core_omega : -m_core_omega;

  const auto f = (1.0 / (2 * m_rank + 1));

  std::vector<double> sum_a(holes.size());
#pragma omp parallel for
  for (std::size_t ia = 0; ia < holes.size(); ia++) {
    const auto &Fa = holes[ia];
    const auto s1 = ((Fa.twoj() - Ff.twoj() + 2 * m_rank) % 4 == 0) ? 1 : -1;
    for (std::size_t im = 0; im < excited.size(); im++) {
      const auto &Fm = excited[im];
      if (t0am[ia][im] == 0.0)
        continue;
      const auto s2 = ((Fa.twoj() - Fm.twoj()) % 4 == 0) ? 1 : -1;
      const auto Wwmva = Coulomb::Wk_abcd(Ff, Fm, Fi, Fa, m_rank);
      const auto Wwavm = Coulomb::Wk_abcd(Ff, Fa, Fi, Fm, m_rank);
      const auto ttam = first_order ? t0am[ia][im] : tam[ia][im];
      const auto ttma = first_order ? t0ma[im][ia] : tma[im][ia];
      const auto A = ttam * Wwmva / (Fa.en() - Fm.en() - ww);
      const auto B = Wwavm * ttma / (Fa.en() - Fm.en() + ww);
      sum_a[ia] += s1 * (A + s2 * B);
    }
  }
  return f * std::accumulate(begin(sum_a), end(sum_a), 0.0);
}

//******************************************************************************
void DiagramRPA::solve_core(const double omega, int max_its, const bool print) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  m_core_omega = std::abs(omega);

  if (holes.empty() || excited.empty())
    return;

  if (m_h->freqDependantQ) {
    // m_h->updateFrequency(m_core_omega); // Cant, is const. must do outside
    setup_ts(m_h);
  }

  if (print) {
    printf("RPA(D) %s (w=%.3f): ", m_h->name().c_str(), omega);
    std::cout << std::flush;
  }

  int it = 2; // 1 iteration means use t0 inside dV
  auto eps = 0.0;
  const auto f = (1.0 / (2 * m_rank + 1));
  for (; it <= max_its; it++) {
    std::vector<double> eps_m(excited.size()); //"thread-safe" eps..?
// XXX "Small" race condition in here???
#pragma omp parallel for
    for (std::size_t im = 0; im < excited.size(); im++) {
      const auto &Fm = excited[im];
      double eps_worst_a = 0.0;
      for (std::size_t ia = 0; ia < holes.size(); ia++) {
        const auto &Fa = holes[ia];

        double sum_am = 0.0;
        double sum_ma = 0.0;

        // Can replace this with dV?? NO. 1) it calcs W (since valence)
        // 2) not thread safe
        for (std::size_t ib = 0; ib < holes.size(); ib++) {
          const auto &Fb = holes[ib];
          const auto s1 =
              ((Fb.twoj() - Fa.twoj() + 2 * m_rank) % 4 == 0) ? 1 : -1;
          const auto s3 =
              ((Fb.twoj() - Fm.twoj() + 2 * m_rank) % 4 == 0) ? 1 : -1;
          for (std::size_t in = 0; in < excited.size(); in++) {
            const auto &Fn = excited[in];

            const auto tbn = tam[ib][in];
            if (tbn == 0.0)
              continue;
            const auto tnb = tma[in][ib];
            const auto tdem = tbn / (Fb.en() - Fn.en() - m_core_omega);
            const auto s2 = ((Fb.twoj() - Fn.twoj()) % 4 == 0) ? 1 : -1;
            const auto stdep = s2 * tnb / (Fb.en() - Fn.en() + m_core_omega);
            // Cast form 'Wtype' (may be double or float) to double
            const auto A = tdem * MyCast<double>(Wanmb[ia][in][im][ib]);
            const auto B = stdep * MyCast<double>(Wabmn[ia][in][im][ib]);
            const auto C = tdem * MyCast<double>(Wmnab[im][in][ia][ib]);
            const auto D = stdep * MyCast<double>(Wmban[im][in][ia][ib]);
            sum_am += s1 * (A + B);
            sum_ma += s3 * (C + D);
          }
        }

        const auto prev = tam[ia][im];
        // 0.5 factor is for damping. f*sum is dV
        tam[ia][im] = 0.5 * (tam[ia][im] + t0am[ia][im] + f * sum_am);
        tma[im][ia] = 0.5 * (tma[im][ia] + t0ma[im][ia] + f * sum_ma);
        const auto delta = std::abs((tam[ia][im] - prev) / tam[ia][im]);
        if (delta > eps_worst_a)
          eps_worst_a = delta;
      } // a (holes)
      eps_m[im] = eps_worst_a;
    } // m (excited)
    // XXX "small" race condition somewhere regarding eps??
    // The itteraion it converges on always seems to be the same..
    // but the value for eps printed changes slightly each run???
    eps = *std::max_element(cbegin(eps_m), cend(eps_m));
    if (eps < eps_targ)
      break;
    // if (print && it % 25 == 0) {
    //   printf("RPA(D) (w=%.3f): %2i %.1e \r", m_core_omega, it, eps);
    //   std::cout << std::flush;
    // }
  } // its
  if (print) {
    printf("%2i %.1e\n", it, eps);
    std::cout << std::flush;
  }
  m_core_eps = eps;
}

} // namespace ExternalField
