#include "Modules/SRBreakdown.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/GenerateOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace Module {

namespace {
const Register r_SRBreakdown{
  "SRBreakdown",
  "Structure radiation for valence transitions, split into Feynman diagrams",
  &SRBreakdown};
} // namespace

//==============================================================================
void SRBreakdown(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"", "Third-order structure radiation for each valence transition, "
          "split into the six Feynman diagrams: direct (a, b, c) and exchange "
          "(ax, bx, cx). Reduced matrix elements <w||h||v>, with omega = "
          "e_w - e_v. The external states are the basis states matching the "
          "valence states."},
     {"operator", "e.g., E1, hfs [E1]"},
     {"options{}", "options specific to operator; blank by dflt"},
     {"rpa", "RPA for the internal operator line: true(=TDHF), false, TDHF, "
             "basis, diagram. Solved at omega = 0 [false]"},
     {"Qk_file",
      "true/false. If true, Coulomb integrals are read from the QkTable file "
      "<identity>.qk.abf (missing ones are calculated, and the file "
      "updated). If false, they are calculated from the Yk table: slower, "
      "but uses much less memory [true]"},
     {"n_min_core", "Minimum n for core states in the diagrams [1]"},
     {"fk", "list: Coulomb screening factors for each k [none]"},
     {"etak", "list: Hole-particle factors for each k [none]"},
     {"n_digits", "Number of decimal places printed [2]"}});
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timer("SRBreakdown");

  const auto oper = input.get<std::string>("operator", "E1");
  const auto h_options =
    input.getBlock("options").value_or(IO::InputBlock(oper, {}));
  const auto h = DiracOperator::generate(oper, h_options, wf);

  const auto rpa_method = input.get<std::string>("rpa", "false");
  const auto use_Qk = input.get("Qk_file", true);
  const auto n_min_core = input.get("n_min_core", 1);
  const auto fk = input.get("fk", std::vector<double>{});
  const auto etak = input.get("etak", std::vector<double>{});
  const auto n_digits = input.get("n_digits", 2);

  auto [core, excited] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel(), n_min_core);
  auto basis = core;
  basis.insert(basis.end(), excited.begin(), excited.end());

  // External states: the basis states matching the valence states
  std::vector<DiracSpinor> legs;
  for (const auto &v : wf.valence()) {
    const auto v_basis = std::find(excited.cbegin(), excited.cend(), v);
    if (v_basis == excited.cend()) {
      std::cout << "Warning: " << v.shortSymbol()
                << " is not in the basis; skipped\n";
      continue;
    }
    legs.push_back(*v_basis);
  }

  std::cout << "\nStructure radiation, broken up by Feynman diagram: "
            << h->name() << "\n";
  std::cout << "Reduced matrix elements <w||h||v>, omega = e_w - e_v ("
            << h->units() << ")\n";
  std::cout << "Core: " << DiracSpinor::state_config(core)
            << " (n >= " << n_min_core << ")\n";
  std::cout << "Excited: " << DiracSpinor::state_config(excited) << "\n";
  std::cout << "Coulomb integrals: " << (use_Qk ? "QkTable" : "YkTable")
            << "\n";
  if (!fk.empty()) {
    std::cout << "Effective Coulomb screening, fk = ";
    for (const auto f : fk) {
      std::cout << f << ", ";
    }
    std::cout << "\n";
  }
  if (!etak.empty()) {
    std::cout << "Effective hole-particle interaction, etak = ";
    for (const auto eta : etak) {
      std::cout << eta << ", ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  // Coulomb integrals
  Coulomb::YkTable yk;
  yk.extend_angular(DiracSpinor::max_tj(basis));
  yk.calculate(basis);
  Coulomb::QkTable qk;
  if (use_Qk) {
    const auto Qk_fname = wf.identity() + ".qk.abf";
    qk.read(Qk_fname);
    const auto existing = qk.count();
    qk.fill(basis, yk);
    if (qk.count() != existing) {
      qk.write(Qk_fname);
    }
  }

  // Operator matrix elements (with RPA, if requested)
  const auto rpa = ExternalField::make_rpa(rpa_method, h.get(), wf.vHF(), true,
                                           wf.basis(), wf.identity());
  if (rpa) {
    rpa->solve_core(0.0);
  }
  const auto h_table =
    Amplitudes::me_table(basis, h.get(), rpa.get(), nullptr, 0.0);
  std::cout << "\n";

  // Values are printed as +d.[n_digits]e-xx; columns at least as wide as the
  // "exchange" header
  const auto width = std::max(n_digits + 7, 8);
  fmt::print("{:5s} {:5s}", "w", "v");
  for (const auto label :
       {"a", "b", "c", "ax", "bx", "cx", "direct", "exchange", "total"}) {
    fmt::print(" {:>{}s}", label, width);
  }
  fmt::print("\n");
  for (std::size_t iv = 0; iv < legs.size(); ++iv) {
    for (std::size_t iw = iv; iw < legs.size(); ++iw) {
      const auto &v = legs[iv];
      const auto &w = legs[iw];
      if (h->isZero(w, v))
        continue;
      const auto omega = w.en() - v.en();
      const auto [a, c, ax, cx] = use_Qk ?
                                    SR_a_c_ax_cx(w, v, omega, core, excited, qk,
                                                 h_table, h->rank(), fk, etak) :
                                    SR_a_c_ax_cx(w, v, omega, core, excited, yk,
                                                 h_table, h->rank(), fk, etak);
      const auto [b, bx] = use_Qk ? SR_b_bx(w, v, omega, core, excited, qk,
                                            h_table, h->rank(), fk, etak) :
                                    SR_b_bx(w, v, omega, core, excited, yk,
                                            h_table, h->rank(), fk, etak);
      const auto direct = a + b + c;
      const auto exchange = ax + bx + cx;
      fmt::print("{:5s} {:5s}", w.shortSymbol(), v.shortSymbol());
      for (const auto value :
           {a, b, c, ax, bx, cx, direct, exchange, direct + exchange}) {
        fmt::print(" {:+{}.{}e}", value, width, n_digits);
      }
      fmt::print("\n");
      std::cout << std::flush;
    }
  }
  std::cout << "\n";
}

//==============================================================================
double factor_k(const std::vector<double> &f, int k) {
  const auto sk = std::size_t(k);
  return sk < f.size() ? f[sk] : 1.0;
}

//==============================================================================
template <typename CoulombTable>
double Q_screened(const CoulombTable &qk, const std::vector<double> &fk, int k,
                  const DiracSpinor &a, const DiracSpinor &b,
                  const DiracSpinor &c, const DiracSpinor &d) {
  return std::sqrt(factor_k(fk, k)) * qk.Q(k, a, b, c, d);
}

template <typename CoulombTable>
double P_screened(const CoulombTable &qk, const std::vector<double> &fk, int k,
                  const DiracSpinor &a, const DiracSpinor &b,
                  const DiracSpinor &c, const DiracSpinor &d) {
  // P^k_abcd = [k] sum_l {a c k; b d l} Q^l_abdc
  const auto [l_min, l_max] = Coulomb::k_minmax_Q(a, b, d, c);
  double pk = 0.0;
  for (int l = l_min; l <= l_max; l += 2) {
    const auto sj =
      Angular::sixj_2(a.twoj(), c.twoj(), 2 * k, b.twoj(), d.twoj(), 2 * l);
    if (sj == 0.0)
      continue;
    pk += sj * Q_screened(qk, fk, l, a, b, d, c);
  }
  return (2 * k + 1) * pk;
}

//==============================================================================
// Each term: sum over internal states (core a,b,c; excited m,n,r) and
// multipolarity mu of
//   (-1)^(mu + k + j_w - j_x) / [mu] {v w k; x y mu} W^mu_(v) W^mu_(w)
//   <x||h||y> / (energy denominators),
// with <x||h||y> the operator line, W_(v) and W_(w) the Coulomb lines
// attached to v and w, and k the operator rank. The final-state energy is
// e_v + omega.
template <typename CoulombTable>
std::array<double, 4>
SR_a_c_ax_cx(const DiracSpinor &w, const DiracSpinor &v, double omega,
             const std::vector<DiracSpinor> &core,
             const std::vector<DiracSpinor> &excited, const CoulombTable &qk,
             const Coulomb::meTable<double> &h, int K,
             const std::vector<double> &fk, const std::vector<double> &etak) {

  double sr_a{0.0}, sr_c{0.0}, sr_ax{0.0}, sr_cx{0.0};

  // + <a||h||m> W_mnvb W_bwna / [(e_vb - e_mn)(e_a - e_m - omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_a, sr_c, sr_ax, sr_cx)
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    for (std::size_t im = 0; im < excited.size(); ++im) {
      const auto &a = core[ia];
      const auto &m = excited[im];
      const auto h_am = h.getv(a, m);
      if (h_am == 0.0)
        continue;
      const auto f_am = Angular::neg1pow_2(2 * K + w.twoj() - a.twoj()) * h_am /
                        (a.en() - m.en() - omega);
      for (const auto &n : excited) {
        for (const auto &b : core) {
          const auto f = f_am / (v.en() + b.en() - m.en() - n.en());
          const auto [min_v, max_v] = Coulomb::k_minmax_W(m, n, v, b);
          const auto [min_w, max_w] = Coulomb::k_minmax_W(b, w, n, a);
          for (int mu = std::max(min_v, min_w); mu <= std::min(max_v, max_w);
               ++mu) {
            const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, a.twoj(),
                                            m.twoj(), 2 * mu);
            if (sj == 0.0)
              continue;
            const auto f_mu =
              factor_k(etak, mu) * f * Angular::neg1pow(mu) * sj / (2 * mu + 1);
            const auto q_v = Q_screened(qk, fk, mu, m, n, v, b);
            const auto p_v = P_screened(qk, fk, mu, m, n, v, b);
            const auto q_w = Q_screened(qk, fk, mu, b, w, n, a);
            const auto p_w = P_screened(qk, fk, mu, b, w, n, a);
            sr_a += f_mu * q_v * q_w;
            sr_ax += f_mu * p_v * q_w;
            sr_cx += f_mu * q_v * p_w;
            sr_c += f_mu * p_v * p_w;
          }
        }
      }
    }
  }

  // - <n||h||m> W_mrva W_awrn / [(e_va - e_mr)(e_va - e_nr + omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_a, sr_c, sr_ax, sr_cx)
  for (std::size_t in = 0; in < excited.size(); ++in) {
    for (std::size_t im = 0; im < excited.size(); ++im) {
      const auto &n = excited[in];
      const auto &m = excited[im];
      const auto h_nm = h.getv(n, m);
      if (h_nm == 0.0)
        continue;
      const auto f_nm = -Angular::neg1pow_2(2 * K + w.twoj() - n.twoj()) * h_nm;
      for (const auto &r : excited) {
        for (const auto &a : core) {
          const auto f = f_nm / ((v.en() + a.en() - m.en() - r.en()) *
                                 (v.en() + a.en() - n.en() - r.en() + omega));
          const auto [min_v, max_v] = Coulomb::k_minmax_W(m, r, v, a);
          const auto [min_w, max_w] = Coulomb::k_minmax_W(a, w, r, n);
          for (int mu = std::max(min_v, min_w); mu <= std::min(max_v, max_w);
               ++mu) {
            const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, n.twoj(),
                                            m.twoj(), 2 * mu);
            if (sj == 0.0)
              continue;
            const auto f_mu =
              factor_k(etak, mu) * f * Angular::neg1pow(mu) * sj / (2 * mu + 1);
            const auto q_v = Q_screened(qk, fk, mu, m, r, v, a);
            const auto p_v = P_screened(qk, fk, mu, m, r, v, a);
            const auto q_w = Q_screened(qk, fk, mu, a, w, r, n);
            const auto p_w = P_screened(qk, fk, mu, a, w, r, n);
            sr_a += f_mu * q_v * q_w;
            sr_ax += f_mu * p_v * q_w;
            sr_cx += f_mu * q_v * p_w;
            sr_c += f_mu * p_v * p_w;
          }
        }
      }
    }
  }

  // - <a||h||m> W_mbvn W_nwba / [(e_vn - e_ab + omega)(e_a - e_m - omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_a, sr_c, sr_ax, sr_cx)
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    for (std::size_t im = 0; im < excited.size(); ++im) {
      const auto &a = core[ia];
      const auto &m = excited[im];
      const auto h_am = h.getv(a, m);
      if (h_am == 0.0)
        continue;
      const auto f_am = -Angular::neg1pow_2(2 * K + w.twoj() - a.twoj()) *
                        h_am / (a.en() - m.en() - omega);
      for (const auto &n : excited) {
        for (const auto &b : core) {
          const auto f = f_am / (v.en() + n.en() - a.en() - b.en() + omega);
          const auto [min_v, max_v] = Coulomb::k_minmax_W(m, b, v, n);
          const auto [min_w, max_w] = Coulomb::k_minmax_W(n, w, b, a);
          for (int mu = std::max(min_v, min_w); mu <= std::min(max_v, max_w);
               ++mu) {
            const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, a.twoj(),
                                            m.twoj(), 2 * mu);
            if (sj == 0.0)
              continue;
            const auto f_mu =
              factor_k(etak, mu) * f * Angular::neg1pow(mu) * sj / (2 * mu + 1);
            const auto q_v = Q_screened(qk, fk, mu, m, b, v, n);
            const auto p_v = P_screened(qk, fk, mu, m, b, v, n);
            const auto q_w = Q_screened(qk, fk, mu, n, w, b, a);
            const auto p_w = P_screened(qk, fk, mu, n, w, b, a);
            sr_a += f_mu * q_v * q_w;
            sr_ax += f_mu * p_v * q_w;
            sr_cx += f_mu * q_v * p_w;
            sr_c += f_mu * p_v * p_w;
          }
        }
      }
    }
  }

  // - <b||h||a> W_acvm W_mwcb / [(e_vm - e_bc + omega)(e_vm - e_ac)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_a, sr_c, sr_ax, sr_cx)
  for (std::size_t ib = 0; ib < core.size(); ++ib) {
    for (std::size_t ia = 0; ia < core.size(); ++ia) {
      const auto &b = core[ib];
      const auto &a = core[ia];
      const auto h_ba = h.getv(b, a);
      if (h_ba == 0.0)
        continue;
      const auto f_ba = -Angular::neg1pow_2(2 * K + w.twoj() - b.twoj()) * h_ba;
      for (const auto &m : excited) {
        for (const auto &c : core) {
          const auto f = f_ba / ((v.en() + m.en() - b.en() - c.en() + omega) *
                                 (v.en() + m.en() - a.en() - c.en()));
          const auto [min_v, max_v] = Coulomb::k_minmax_W(a, c, v, m);
          const auto [min_w, max_w] = Coulomb::k_minmax_W(m, w, c, b);
          for (int mu = std::max(min_v, min_w); mu <= std::min(max_v, max_w);
               ++mu) {
            const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, b.twoj(),
                                            a.twoj(), 2 * mu);
            if (sj == 0.0)
              continue;
            const auto f_mu =
              factor_k(etak, mu) * f * Angular::neg1pow(mu) * sj / (2 * mu + 1);
            const auto q_v = Q_screened(qk, fk, mu, a, c, v, m);
            const auto p_v = P_screened(qk, fk, mu, a, c, v, m);
            const auto q_w = Q_screened(qk, fk, mu, m, w, c, b);
            const auto p_w = P_screened(qk, fk, mu, m, w, c, b);
            sr_a += f_mu * q_v * q_w;
            sr_ax += f_mu * p_v * q_w;
            sr_cx += f_mu * q_v * p_w;
            sr_c += f_mu * p_v * p_w;
          }
        }
      }
    }
  }

  // + <m||h||a> W_anvb W_bwnm / [(e_a - e_m + omega)(e_vb - e_mn + omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_a, sr_c, sr_ax, sr_cx)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    for (std::size_t ia = 0; ia < core.size(); ++ia) {
      const auto &m = excited[im];
      const auto &a = core[ia];
      const auto h_ma = h.getv(m, a);
      if (h_ma == 0.0)
        continue;
      const auto f_ma = Angular::neg1pow_2(2 * K + w.twoj() - m.twoj()) * h_ma /
                        (a.en() - m.en() + omega);
      for (const auto &n : excited) {
        for (const auto &b : core) {
          const auto f = f_ma / (v.en() + b.en() - m.en() - n.en() + omega);
          const auto [min_v, max_v] = Coulomb::k_minmax_W(a, n, v, b);
          const auto [min_w, max_w] = Coulomb::k_minmax_W(b, w, n, m);
          for (int mu = std::max(min_v, min_w); mu <= std::min(max_v, max_w);
               ++mu) {
            const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, m.twoj(),
                                            a.twoj(), 2 * mu);
            if (sj == 0.0)
              continue;
            const auto f_mu =
              factor_k(etak, mu) * f * Angular::neg1pow(mu) * sj / (2 * mu + 1);
            const auto q_v = Q_screened(qk, fk, mu, a, n, v, b);
            const auto p_v = P_screened(qk, fk, mu, a, n, v, b);
            const auto q_w = Q_screened(qk, fk, mu, b, w, n, m);
            const auto p_w = P_screened(qk, fk, mu, b, w, n, m);
            sr_a += f_mu * q_v * q_w;
            sr_ax += f_mu * p_v * q_w;
            sr_cx += f_mu * q_v * p_w;
            sr_c += f_mu * p_v * p_w;
          }
        }
      }
    }
  }

  // - <m||h||a> W_abvn W_nwbm / [(e_a - e_m + omega)(e_vn - e_ab)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_a, sr_c, sr_ax, sr_cx)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    for (std::size_t ia = 0; ia < core.size(); ++ia) {
      const auto &m = excited[im];
      const auto &a = core[ia];
      const auto h_ma = h.getv(m, a);
      if (h_ma == 0.0)
        continue;
      const auto f_ma = -Angular::neg1pow_2(2 * K + w.twoj() - m.twoj()) *
                        h_ma / (a.en() - m.en() + omega);
      for (const auto &n : excited) {
        for (const auto &b : core) {
          const auto f = f_ma / (v.en() + n.en() - a.en() - b.en());
          const auto [min_v, max_v] = Coulomb::k_minmax_W(a, b, v, n);
          const auto [min_w, max_w] = Coulomb::k_minmax_W(n, w, b, m);
          for (int mu = std::max(min_v, min_w); mu <= std::min(max_v, max_w);
               ++mu) {
            const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, m.twoj(),
                                            a.twoj(), 2 * mu);
            if (sj == 0.0)
              continue;
            const auto f_mu =
              factor_k(etak, mu) * f * Angular::neg1pow(mu) * sj / (2 * mu + 1);
            const auto q_v = Q_screened(qk, fk, mu, a, b, v, n);
            const auto p_v = P_screened(qk, fk, mu, a, b, v, n);
            const auto q_w = Q_screened(qk, fk, mu, n, w, b, m);
            const auto p_w = P_screened(qk, fk, mu, n, w, b, m);
            sr_a += f_mu * q_v * q_w;
            sr_ax += f_mu * p_v * q_w;
            sr_cx += f_mu * q_v * p_w;
            sr_c += f_mu * p_v * p_w;
          }
        }
      }
    }
  }

  return {sr_a, sr_c, sr_ax, sr_cx};
}

//==============================================================================
// Each term: sum over internal states (core a,b,c; excited m,n,r) and
// multipolarities mu, lambda of
//   (phase) {6j} {6j} Q^mu W^lambda <x||h||y> / (energy denominators),
// with <x||h||y> the operator line, Q the direct Coulomb line, W the
// antisymmetrised one, and k the operator rank. The final-state energy is
// e_v + omega.
template <typename CoulombTable>
std::array<double, 2>
SR_b_bx(const DiracSpinor &w, const DiracSpinor &v, double omega,
        const std::vector<DiracSpinor> &core,
        const std::vector<DiracSpinor> &excited, const CoulombTable &qk,
        const Coulomb::meTable<double> &h, int K, const std::vector<double> &fk,
        const std::vector<double> &etak) {

  double sr_b{0.0}, sr_bx{0.0};

  // + (-1)^(k+v-m+r-n) {v w k; lambda mu m} {k a n; r lambda mu}
  //   <a||h||n> Q_mrva W_nwrm / [(e_va - e_mr)(e_a - e_n - omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_b, sr_bx)
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    for (std::size_t in = 0; in < excited.size(); ++in) {
      const auto &a = core[ia];
      const auto &n = excited[in];
      const auto h_an = h.getv(a, n);
      if (h_an == 0.0)
        continue;
      const auto f_an = h_an / (a.en() - n.en() - omega);
      for (const auto &m : excited) {
        for (const auto &r : excited) {
          const auto f = Angular::neg1pow_2(2 * K + v.twoj() - m.twoj() +
                                            r.twoj() - n.twoj()) *
                         f_an / (v.en() + a.en() - m.en() - r.en());
          const auto [min_mu, max_mu] = Coulomb::k_minmax_Q(m, r, v, a);
          const auto [min_l, max_l] = Coulomb::k_minmax_W(n, w, r, m);
          for (int mu = min_mu; mu <= max_mu; mu += 2) {
            const auto q_mu = Q_screened(qk, fk, mu, m, r, v, a);
            if (q_mu == 0.0)
              continue;
            const auto f_mu = factor_k(etak, mu) * f * q_mu;
            for (int lambda = min_l; lambda <= max_l; ++lambda) {
              const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K,
                                              2 * lambda, 2 * mu, m.twoj()) *
                              Angular::sixj_2(2 * K, a.twoj(), n.twoj(),
                                              r.twoj(), 2 * lambda, 2 * mu);
              if (sj == 0.0)
                continue;
              sr_b += f_mu * sj * Q_screened(qk, fk, lambda, n, w, r, m);
              sr_bx += f_mu * sj * P_screened(qk, fk, lambda, n, w, r, m);
            }
          }
        }
      }
    }
  }

  // - (-1)^(k+v-m+n-b) {v w k; lambda mu m} {k a b; n lambda mu}
  //   <a||h||b> Q_mnva W_bwnm / [(e_va - e_mn)(e_vb - e_mn + omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_b, sr_bx)
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    for (std::size_t ib = 0; ib < core.size(); ++ib) {
      const auto &a = core[ia];
      const auto &b = core[ib];
      const auto h_ab = h.getv(a, b);
      if (h_ab == 0.0)
        continue;
      for (const auto &m : excited) {
        for (const auto &n : excited) {
          const auto f = -Angular::neg1pow_2(2 * K + v.twoj() - m.twoj() +
                                             n.twoj() - b.twoj()) *
                         h_ab /
                         ((v.en() + a.en() - m.en() - n.en()) *
                          (v.en() + b.en() - m.en() - n.en() + omega));
          const auto [min_mu, max_mu] = Coulomb::k_minmax_Q(m, n, v, a);
          const auto [min_l, max_l] = Coulomb::k_minmax_W(b, w, n, m);
          for (int mu = min_mu; mu <= max_mu; mu += 2) {
            const auto q_mu = Q_screened(qk, fk, mu, m, n, v, a);
            if (q_mu == 0.0)
              continue;
            const auto f_mu = factor_k(etak, mu) * f * q_mu;
            for (int lambda = min_l; lambda <= max_l; ++lambda) {
              const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K,
                                              2 * lambda, 2 * mu, m.twoj()) *
                              Angular::sixj_2(2 * K, a.twoj(), b.twoj(),
                                              n.twoj(), 2 * lambda, 2 * mu);
              if (sj == 0.0)
                continue;
              sr_b += f_mu * sj * Q_screened(qk, fk, lambda, b, w, n, m);
              sr_bx += f_mu * sj * P_screened(qk, fk, lambda, b, w, n, m);
            }
          }
        }
      }
    }
  }

  // - (-1)^(k+v-m+c-a) {v w k; mu lambda a} {k m b; c lambda mu}
  //   <b||h||m> Q_mwca W_acvb / [(e_vm - e_ac + omega)(e_b - e_m - omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_b, sr_bx)
  for (std::size_t ib = 0; ib < core.size(); ++ib) {
    for (std::size_t im = 0; im < excited.size(); ++im) {
      const auto &b = core[ib];
      const auto &m = excited[im];
      const auto h_bm = h.getv(b, m);
      if (h_bm == 0.0)
        continue;
      const auto f_bm = -h_bm / (b.en() - m.en() - omega);
      for (const auto &c : core) {
        for (const auto &a : core) {
          const auto f = Angular::neg1pow_2(2 * K + v.twoj() - m.twoj() +
                                            c.twoj() - a.twoj()) *
                         f_bm / (v.en() + m.en() - a.en() - c.en() + omega);
          const auto [min_mu, max_mu] = Coulomb::k_minmax_Q(m, w, c, a);
          const auto [min_l, max_l] = Coulomb::k_minmax_W(a, c, v, b);
          for (int mu = min_mu; mu <= max_mu; mu += 2) {
            const auto q_mu = Q_screened(qk, fk, mu, m, w, c, a);
            if (q_mu == 0.0)
              continue;
            const auto f_mu = factor_k(etak, mu) * f * q_mu;
            for (int lambda = min_l; lambda <= max_l; ++lambda) {
              const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, 2 * mu,
                                              2 * lambda, a.twoj()) *
                              Angular::sixj_2(2 * K, m.twoj(), b.twoj(),
                                              c.twoj(), 2 * lambda, 2 * mu);
              if (sj == 0.0)
                continue;
              sr_b += f_mu * sj * Q_screened(qk, fk, lambda, a, c, v, b);
              sr_bx += f_mu * sj * P_screened(qk, fk, lambda, a, c, v, b);
            }
          }
        }
      }
    }
  }

  // - (-1)^(k+v-n+b-a) {v w k; mu lambda a} {k n m; b lambda mu}
  //   <m||h||n> Q_nwba W_abvm / [(e_vn - e_ab + omega)(e_vm - e_ab)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_b, sr_bx)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    for (std::size_t in = 0; in < excited.size(); ++in) {
      const auto &m = excited[im];
      const auto &n = excited[in];
      const auto h_mn = h.getv(m, n);
      if (h_mn == 0.0)
        continue;
      for (const auto &b : core) {
        for (const auto &a : core) {
          const auto f = -Angular::neg1pow_2(2 * K + v.twoj() - n.twoj() +
                                             b.twoj() - a.twoj()) *
                         h_mn /
                         ((v.en() + n.en() - a.en() - b.en() + omega) *
                          (v.en() + m.en() - a.en() - b.en()));
          const auto [min_mu, max_mu] = Coulomb::k_minmax_Q(n, w, b, a);
          const auto [min_l, max_l] = Coulomb::k_minmax_W(a, b, v, m);
          for (int mu = min_mu; mu <= max_mu; mu += 2) {
            const auto q_mu = Q_screened(qk, fk, mu, n, w, b, a);
            if (q_mu == 0.0)
              continue;
            const auto f_mu = factor_k(etak, mu) * f * q_mu;
            for (int lambda = min_l; lambda <= max_l; ++lambda) {
              const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, 2 * mu,
                                              2 * lambda, a.twoj()) *
                              Angular::sixj_2(2 * K, n.twoj(), m.twoj(),
                                              b.twoj(), 2 * lambda, 2 * mu);
              if (sj == 0.0)
                continue;
              sr_b += f_mu * sj * Q_screened(qk, fk, lambda, a, b, v, m);
              sr_bx += f_mu * sj * P_screened(qk, fk, lambda, a, b, v, m);
            }
          }
        }
      }
    }
  }

  // + (-1)^(k+v-a+m-r) {v w k; lambda mu m} {k n a; r lambda mu}
  //   <n||h||a> Q_mrvn W_awrm / [(e_a - e_n + omega)(e_va - e_mr + omega)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_b, sr_bx)
  for (std::size_t in = 0; in < excited.size(); ++in) {
    for (std::size_t ia = 0; ia < core.size(); ++ia) {
      const auto &n = excited[in];
      const auto &a = core[ia];
      const auto h_na = h.getv(n, a);
      if (h_na == 0.0)
        continue;
      const auto f_na = h_na / (a.en() - n.en() + omega);
      for (const auto &m : excited) {
        for (const auto &r : excited) {
          const auto f = Angular::neg1pow_2(2 * K + v.twoj() - a.twoj() +
                                            m.twoj() - r.twoj()) *
                         f_na / (v.en() + a.en() - m.en() - r.en() + omega);
          const auto [min_mu, max_mu] = Coulomb::k_minmax_Q(m, r, v, n);
          const auto [min_l, max_l] = Coulomb::k_minmax_W(a, w, r, m);
          for (int mu = min_mu; mu <= max_mu; mu += 2) {
            const auto q_mu = Q_screened(qk, fk, mu, m, r, v, n);
            if (q_mu == 0.0)
              continue;
            const auto f_mu = factor_k(etak, mu) * f * q_mu;
            for (int lambda = min_l; lambda <= max_l; ++lambda) {
              const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K,
                                              2 * lambda, 2 * mu, m.twoj()) *
                              Angular::sixj_2(2 * K, n.twoj(), a.twoj(),
                                              r.twoj(), 2 * lambda, 2 * mu);
              if (sj == 0.0)
                continue;
              sr_b += f_mu * sj * Q_screened(qk, fk, lambda, a, w, r, m);
              sr_bx += f_mu * sj * P_screened(qk, fk, lambda, a, w, r, m);
            }
          }
        }
      }
    }
  }

  // - (-1)^(k+v-b+c-a) {v w k; mu lambda a} {k m b; c mu lambda}
  //   <m||h||b> Q_bwca W_acvm / [(e_b - e_m + omega)(e_vm - e_ac)]
#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(+ : sr_b, sr_bx)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    for (std::size_t ib = 0; ib < core.size(); ++ib) {
      const auto &m = excited[im];
      const auto &b = core[ib];
      const auto h_mb = h.getv(m, b);
      if (h_mb == 0.0)
        continue;
      const auto f_mb = -h_mb / (b.en() - m.en() + omega);
      for (const auto &c : core) {
        for (const auto &a : core) {
          const auto f = Angular::neg1pow_2(2 * K + v.twoj() - b.twoj() +
                                            c.twoj() - a.twoj()) *
                         f_mb / (v.en() + m.en() - a.en() - c.en());
          const auto [min_mu, max_mu] = Coulomb::k_minmax_Q(b, w, c, a);
          const auto [min_l, max_l] = Coulomb::k_minmax_W(a, c, v, m);
          for (int mu = min_mu; mu <= max_mu; mu += 2) {
            const auto q_mu = Q_screened(qk, fk, mu, b, w, c, a);
            if (q_mu == 0.0)
              continue;
            const auto f_mu = factor_k(etak, mu) * f * q_mu;
            for (int lambda = min_l; lambda <= max_l; ++lambda) {
              const auto sj = Angular::sixj_2(v.twoj(), w.twoj(), 2 * K, 2 * mu,
                                              2 * lambda, a.twoj()) *
                              Angular::sixj_2(2 * K, m.twoj(), b.twoj(),
                                              c.twoj(), 2 * mu, 2 * lambda);
              if (sj == 0.0)
                continue;
              sr_b += f_mu * sj * Q_screened(qk, fk, lambda, a, c, v, m);
              sr_bx += f_mu * sj * P_screened(qk, fk, lambda, a, c, v, m);
            }
          }
        }
      }
    }
  }

  return {sr_b, sr_bx};
}

//==============================================================================
template std::array<double, 4> SR_a_c_ax_cx(
  const DiracSpinor &w, const DiracSpinor &v, double omega,
  const std::vector<DiracSpinor> &core, const std::vector<DiracSpinor> &excited,
  const Coulomb::YkTable &qk, const Coulomb::meTable<double> &h, int K,
  const std::vector<double> &fk, const std::vector<double> &etak);
template std::array<double, 4> SR_a_c_ax_cx(
  const DiracSpinor &w, const DiracSpinor &v, double omega,
  const std::vector<DiracSpinor> &core, const std::vector<DiracSpinor> &excited,
  const Coulomb::QkTable &qk, const Coulomb::meTable<double> &h, int K,
  const std::vector<double> &fk, const std::vector<double> &etak);
template std::array<double, 2>
SR_b_bx(const DiracSpinor &w, const DiracSpinor &v, double omega,
        const std::vector<DiracSpinor> &core,
        const std::vector<DiracSpinor> &excited, const Coulomb::YkTable &qk,
        const Coulomb::meTable<double> &h, int K, const std::vector<double> &fk,
        const std::vector<double> &etak);
template std::array<double, 2>
SR_b_bx(const DiracSpinor &w, const DiracSpinor &v, double omega,
        const std::vector<DiracSpinor> &core,
        const std::vector<DiracSpinor> &excited, const Coulomb::QkTable &qk,
        const Coulomb::meTable<double> &h, int K, const std::vector<double> &fk,
        const std::vector<double> &etak);
template double Q_screened(const Coulomb::YkTable &qk,
                           const std::vector<double> &fk, int k,
                           const DiracSpinor &a, const DiracSpinor &b,
                           const DiracSpinor &c, const DiracSpinor &d);
template double Q_screened(const Coulomb::QkTable &qk,
                           const std::vector<double> &fk, int k,
                           const DiracSpinor &a, const DiracSpinor &b,
                           const DiracSpinor &c, const DiracSpinor &d);
template double P_screened(const Coulomb::YkTable &qk,
                           const std::vector<double> &fk, int k,
                           const DiracSpinor &a, const DiracSpinor &b,
                           const DiracSpinor &c, const DiracSpinor &d);
template double P_screened(const Coulomb::QkTable &qk,
                           const std::vector<double> &fk, int k,
                           const DiracSpinor &a, const DiracSpinor &b,
                           const DiracSpinor &c, const DiracSpinor &d);

} // namespace Module
