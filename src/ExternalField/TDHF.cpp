#include "TDHF.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <fstream>
#include <memory>
#include <vector>

namespace ExternalField {

//==============================================================================
TDHF::TDHF(const DiracOperator::TensorOperator *const h,
           const HF::HartreeFock *const hf)
    : CorePolarisation((assert(h != nullptr), h)),
      p_hf((assert(hf != nullptr), hf)),
      m_core(hf->core()),
      m_alpha(hf->alpha()),
      p_VBr(hf->vBreit()) {
  initialise_dPsi();
}

//==============================================================================
void TDHF::initialise_dPsi() {
  // Initialise dPsi vectors, accounting for selection rules
  constexpr bool print = false;
  m_X.resize(m_core.size());
  for (auto ic = 0u; ic < m_core.size(); ic++) {
    const auto &Fc = m_core[ic];
    const auto pi_ch = Fc.parity() * m_pi;
    const auto tj_c = Fc.twoj();
    const auto tjmin_tmp = tj_c - 2 * m_rank;
    const auto tjmin = tjmin_tmp < 1 ? 1 : tjmin_tmp;
    const auto tjmax = tj_c + 2 * m_rank;
    if constexpr (print) {
      std::cout << "|" << Fc.symbol() << ">  -->  ";
    }
    for (int tj = tjmin; tj <= tjmax; tj += 2) {
      const auto l_minus = (tj - 1) / 2;
      const auto pi_chla = Angular::parity_l(l_minus) * pi_ch;

      if (Angular::triangle(tj_c, tj, 2 * m_rank) == 0)
        continue;

      const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
      const auto kappa = Angular::kappa_twojl(tj, l);
      m_X[ic].emplace_back(0, kappa, Fc.grid_sptr());
      m_X[ic].back().max_pt() = Fc.max_pt();
      if constexpr (print) {
        std::cout << "|" << m_X[ic].back().symbol() << "> + ";
      }
    }
    if constexpr (print) {
      std::cout << "\n";
    }
  }
  m_Y = m_X;
}

//==============================================================================
void TDHF::clear() {
  using namespace qip::overloads;
  m_X *= 0.0;
  m_Y *= 0.0;
}

//==============================================================================
const std::vector<DiracSpinor> &TDHF::get_dPsis(const DiracSpinor &Fc,
                                                dPsiType XorY) const {

  const auto index = static_cast<std::size_t>(
      std::find(m_core.cbegin(), m_core.cend(), Fc) - m_core.cbegin());
  assert(index < m_X.size());
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

//==============================================================================
const DiracSpinor &TDHF::get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                    const int kappa_x) const {
  const auto &dPsis = get_dPsis(Fc, XorY);
  auto match_kappa_x = [=](const auto &Fa) { return Fa.kappa() == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//==============================================================================
std::vector<DiracSpinor>
TDHF::solve_dPsis(const DiracSpinor &Fv, const double omega, dPsiType XorY,
                  const MBPT::CorrelationPotential *const Sigma, StateType st,
                  bool incl_dV) const {
  std::vector<DiracSpinor> dFvs;
  const auto tjmin = std::max(1, Fv.twoj() - 2 * m_rank);
  const auto tjmax = Fv.twoj() + 2 * m_rank;
  for (int tjbeta = tjmin; tjbeta <= tjmax; tjbeta += 2) {
    const auto kappa = Angular::kappa_twojpi(tjbeta, Fv.parity() * m_pi);
    dFvs.push_back(solve_dPsi(Fv, omega, XorY, kappa, Sigma, st, incl_dV));
  }
  return dFvs;
}
//==============================================================================
DiracSpinor TDHF::solve_dPsi(const DiracSpinor &Fv, const double omega,
                             dPsiType XorY, const int kappa_x,
                             const MBPT::CorrelationPotential *const Sigma,
                             StateType st, bool incl_dV) const {
  // Solves (H + Sigma - e - w)X = -(h + dV - de)Psi
  // or     (H + Sigma - e + w)Y = -(h^dag + dV^dag - de)Psi

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;

  const auto imag = m_h->imaginaryQ();

  auto rhs = m_h->reduced_rhs(kappa_x, Fv);
  // if (imag && conj)
  //   rhs *= -1;
  /// XXX Think this ^ shouldn't be here. Haven't checked though!

  if (incl_dV)
    rhs += dV_rhs(kappa_x, Fv, conj);
  if (kappa_x == Fv.kappa() && !imag) {
    auto de = m_h->reducedME(Fv, Fv);
    if (incl_dV)
      de += dV(Fv, Fv, conj);
    rhs -= de * Fv;
  }

  // Do we want |Y> or <Y| ?
  auto s2 = 1;
  if (st == StateType::bra) {
    // "left-hand-side" : "reduced" ket, so has factor (+ confugate)
    const auto sj =
        Angular::evenQ_2(Fv.twoj() - Angular::twoj_k(kappa_x)) ? 1 : -1;
    // if conj, extra * => +1
    const auto si = imag && !conj ? -1 : 1;
    s2 = sj * si;
  }

  const auto vl = p_hf->vlocal(Angular::l_k(Fv.kappa()));
  const auto &Hmag = p_hf->Hmag(Angular::l_k(Fv.kappa()));
  // The l from X ? or from Fv ?
  return s2 * ExternalField::solveMixedState(Fv, ww, vl, m_alpha, m_core, rhs,
                                             1.0e-9, Sigma, p_VBr, Hmag);
}

//==============================================================================
void TDHF::solve_ms_core(std::vector<DiracSpinor> &dFb, const DiracSpinor &Fb,
                         const std::vector<DiracSpinor> &hFbs,
                         const double omega, dPsiType XorY,
                         double eps_ms) const {
  // Solves (H - e - w)Xb = -(h + dV - de)Psi
  // or     (H - e + w)Y = -(h^dag + dV^dag - de)Psi

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;

  const auto imag = m_h->imaginaryQ();
  for (auto ibeta = 0ul; ibeta < dFb.size(); ibeta++) {

    auto &dF_beta = dFb[ibeta];
    const int kappa_beta = dF_beta.kappa();

    const auto &hFb = hFbs[ibeta];
    const auto s = (imag && conj) ? -1.0 : 1.0;
    auto rhs = s * hFb + dV_rhs(kappa_beta, Fb, conj);
    if (kappa_beta == Fb.kappa() && !imag) {
      const auto de = Fb * rhs;
      rhs -= de * Fb;
    }

    const auto vl = p_hf->vlocal(Angular::l_k(Fb.kappa()));
    const auto &Hmag = p_hf->Hmag(Angular::l_k(Fb.kappa()));
    // The l from X ? or from Fv ?
    ExternalField::solveMixedState(dF_beta, Fb, ww, vl, m_alpha, m_core, rhs,
                                   eps_ms, nullptr, p_VBr, Hmag);
  }
}

//==============================================================================
std::vector<std::vector<DiracSpinor>> TDHF::form_hFcore() const {
  std::vector<std::vector<DiracSpinor>> hFcore(m_core.size());
  for (auto ic = 0ul; ic < m_X.size(); ic++) {
    const auto &Fc = m_core[ic];
    hFcore.reserve(m_X[ic].size()); // each h projection
    for (auto beta = 0ul; beta < m_X[ic].size(); beta++) {
      const auto &Xx = m_X[ic][beta];
      hFcore[ic].push_back(m_h->reduced_rhs(Xx.kappa(), Fc));
    }
  }
  return hFcore;
}

//==============================================================================
std::pair<double, std::string> TDHF::tdhf_core_it(double omega,
                                                  double eta_damp) {
  // solve TDHF equation for core - single iteration
  using namespace qip::overloads;
  const bool staticQ = std::abs(omega) < 1.0e-10;
  const auto s = m_imag ? -1 : 1;

  // make copies (instead of blank spinors) since the solve_mixed_states
  // is faster if it starts from existing solution
  auto Xs = m_X;
  auto Ys = m_Y;
  std::vector<std::pair<double, std::string>> epss(m_core.size());

  const auto eps_ms = 1.0e-10;

  double DdF2 = 0.0;
  double dF2 = 0.0;
#pragma omp parallel for reduction(+ : DdF2) reduction(+ : dF2)
  for (auto ib = 0ul; ib < m_core.size(); ib++) {
    const auto &Fb = m_core.at(ib);

    const auto &hFbs = m_hFcore[ib];

    // solve for dF, and damp
    solve_ms_core(Xs[ib], Fb, hFbs, omega, dPsiType::X, eps_ms);
    Xs[ib] = eta_damp * m_X[ib] + (1.0 - eta_damp) * Xs[ib];
    if (staticQ) {
      Ys[ib] = s * Xs[ib];
    } else {
      solve_ms_core(Ys[ib], Fb, hFbs, omega, dPsiType::Y, eps_ms);
      Ys[ib] = eta_damp * m_Y[ib] + (1.0 - eta_damp) * Ys[ib];
    }

    // find eps: of each orbital (just so I know which was the 'worst')
    // and for total core (for overal convergance)
    epss[ib] = {0.0, ""};
    for (std::size_t ibeta = 0; ibeta < Xs[ib].size(); ++ibeta) {
      const auto &Xx = Xs[ib][ibeta];
      const auto &oldX = m_X[ib][ibeta];

      const auto t_DdF2 = (Xx - oldX).norm2();
      const auto t_dF2 = 0.5 * (Xx + oldX).norm2();
      DdF2 += t_DdF2;
      dF2 += t_dF2;

      // find worst orbital
      const auto t_eps = t_dF2 == 0.0 ? 1.0 : t_DdF2 / t_dF2;
      if (t_eps > epss[ib].first) {
        epss[ib].first = t_eps;
        epss[ib].second = Fb.shortSymbol() + "," + Xx.shortSymbol();
      }
    }
  }
  const auto total_eps = DdF2 / dF2;

  // Find the worst orbital (just for reporting)
  const auto comp_first = [](const auto &a, const auto &b) {
    return a.first < b.first;
  };
  const auto worst_orbital =
      *std::max_element(epss.begin(), epss.end(), comp_first);

  m_X = std::move(Xs);
  m_Y = std::move(Ys);

  return {total_eps, worst_orbital.second};
}

//==============================================================================
void TDHF::solve_core(double omega, int max_its, const bool print) {
  const double converge_targ = m_eps;
  const auto eta_damp = m_eta;
  omega = std::abs(omega);

  if (print) {
    printf("TDHF %s (w=%.4f): ", m_h->name().c_str(), omega);
    std::cout << std::flush;
  }

  m_hFcore = form_hFcore();

  std::pair<double, std::string> eps{};
  double best_eps{1.0};
  int count_worse = 0;
  int it{0};
  for (; it < max_its; it++) {
    const auto eta = it == 0 ? 0.0 : eta_damp;
    eps = tdhf_core_it(omega, eta);

    // Check for a "platau" in convergance (count # of 'worse' iterations)
    if (it > 15) {
      if (eps.first > best_eps) {
        ++count_worse;
      } else {
        best_eps = eps.first;
        count_worse = 0;
      }
    }

    if (it > 1 && eps.first < converge_targ)
      break;
    if (count_worse > 5 && eps.first < 1.0e4 * converge_targ)
      break;
    if (std::isnan(eps.first)) {
      break;
    }
  }

  if (print) {
    printf("%2i %.1e [%s]", it, eps.first, eps.second.c_str());
    if (eps.first > 1.0e-8 && max_its > 1)
      std::cout << "  *";
    if (eps.first > 1.0e-6 && max_its > 1)
      std::cout << "*";
    if (eps.first > 1.0e-4 && max_its > 1)
      std::cout << "*";
    std::cout << "\n" << std::flush;
  }

  // set last eps (convergance) and frequency (omega)
  m_core_eps = eps.first;
  m_core_its = it;
  m_core_omega = omega;
}

//==============================================================================
// does it matter if a or b is in the core?
double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm, bool conj) const {
  const auto s = conj && m_h->imaginaryQ() ? -1 : 1; // careful. OK?
  return s * Fn * dV_rhs(Fn.kappa(), Fm, conj);
}

//==============================================================================
double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm) const {
  const auto conj = Fm.en() > Fn.en();
  return dV(Fn, Fm, conj);
}

//==============================================================================
DiracSpinor TDHF::dV_rhs(const int kappa_n, const DiracSpinor &Fa,
                         bool conj) const {

  auto dVFa = DiracSpinor(0, kappa_n, Fa.grid_sptr());
  dVFa.max_pt() = Fa.max_pt();

  const auto ChiType = !conj ? dPsiType::X : dPsiType::Y;
  const auto EtaType = !conj ? dPsiType::Y : dPsiType::X;

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);
  const auto tjn = Angular::twoj_k(kappa_n);

  // nb: faster to not //ize this one
  for (const auto &Fb : m_core) {

    const auto &X_b = get_dPsis(Fb, ChiType);
    const auto &Y_b = get_dPsis(Fb, EtaType);

    for (auto ibeta = 0ul; ibeta < X_b.size(); ++ibeta) {
      const auto &X_beta = X_b[ibeta];
      const auto &Y_beta = Y_b[ibeta];

      const auto sQ = Angular::neg1pow_2(tjn - X_beta.twoj() + 2 * k);
      if (sQ == 1) { // faster than multiplying!
        dVFa += Coulomb::Wkv_bcd(k, kappa_n, Fb, Fa, X_beta);
        dVFa += Coulomb::Wkv_bcd(k, kappa_n, Y_beta, Fa, Fb);
      } else {
        dVFa -= Coulomb::Wkv_bcd(k, kappa_n, Fb, Fa, X_beta);
        dVFa -= Coulomb::Wkv_bcd(k, kappa_n, Y_beta, Fa, Fb);
      }

      // Breit part:
      if (p_VBr) {
        // Makes small difference for E1, huge difference for HFS
        dVFa += tkp1 * p_VBr->dV_Br(kappa_n, k, Fa, Fb, X_beta, Y_beta);
      }
    }
  }

  dVFa *= (1.0 / tkp1);

  return dVFa;
}

} // namespace ExternalField
