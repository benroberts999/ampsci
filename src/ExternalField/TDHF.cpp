#include "TDHF.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <fstream>
#include <memory>
#include <vector>

namespace ExternalField {

//******************************************************************************
TDHF::TDHF(const DiracOperator::TensorOperator *const h,
           const HF::HartreeFock *const hf)
    : CorePolarisation(h),
      p_hf(hf),
      m_core(hf->get_core()),
      // m_vl(hf->get_vlocal(0)), //udpated, to use l-dependent QED [no change]
      m_Hmag(hf->get_Hrad_mag(0)), //(same each l)
      m_alpha(hf->get_alpha()),
      p_VBr(hf->get_Breit())
// Add check for null hf?
{
  initialise_dPsi();
}

//******************************************************************************
void TDHF::initialise_dPsi() {
  // Initialise dPsi vectors, accounting for selection rules
  const bool print = false;
  m_X.resize(m_core.size());
  for (auto ic = 0u; ic < m_core.size(); ic++) {
    const auto &Fc = m_core[ic];
    const auto pi_ch = Fc.parity() * m_pi;
    const auto tj_c = Fc.twoj();
    const auto tjmin_tmp = tj_c - 2 * m_rank;
    const auto tjmin = tjmin_tmp < 1 ? 1 : tjmin_tmp;
    const auto tjmax = tj_c + 2 * m_rank;
    if (print)
      std::cout << "|" << Fc.symbol() << ">  -->  ";
    for (int tj = tjmin; tj <= tjmax; tj += 2) {
      const auto l_minus = (tj - 1) / 2;
      const auto pi_chla = Angular::parity_l(l_minus) * pi_ch;
      const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
      const auto kappa = Angular::kappa_twojl(tj, l);
      m_X[ic].emplace_back(0, kappa, Fc.rgrid);
      m_X[ic].back().set_max_pt() = Fc.max_pt();
      if (print)
        std::cout << "|" << m_X[ic].back().symbol() << "> + ";
    }
    if (print)
      std::cout << "\n";
  }
  m_Y = m_X;
  m_X0 = m_X;
  m_Y0 = m_X;
}

//******************************************************************************
void TDHF::clear() {
  // re-set p0/pinf? no need.
  for (auto &mx : m_X) {
    for (auto &m : mx) {
      m *= 0.0;
    }
  }
  m_Y = m_X;
  m_X0 = m_X;
  m_Y0 = m_X;
}

//******************************************************************************
const std::vector<DiracSpinor> &TDHF::get_dPsis(const DiracSpinor &Fc,
                                                dPsiType XorY) const {

  const auto index = static_cast<std::size_t>(
      std::find(m_core.cbegin(), m_core.cend(), Fc) - m_core.cbegin());
  // Note: no bounds checking here! Used in critical loop. Better way?
  assert(index < m_X.size());
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

const std::vector<DiracSpinor> &TDHF::get_dPsis_0(const DiracSpinor &Fc,
                                                  dPsiType XorY) const {

  // return solve_dPsis(Fc, m_core_omega, XorY, nullptr, StateType::ket, false);

  const auto index = static_cast<std::size_t>(
      std::find(m_core.cbegin(), m_core.cend(), Fc) - m_core.cbegin());
  // Note: no bounds checking here! Used in critical loop. Better way?
  assert(index < m_X.size());
  return XorY == dPsiType::X ? m_X0[index] : m_Y0[index];
}

//******************************************************************************
const DiracSpinor &TDHF::get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                    const int kappa_x) const {
  const auto &dPsis = get_dPsis(Fc, XorY);
  auto match_kappa_x = [=](const auto &Fa) { return Fa.k == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//******************************************************************************
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
//******************************************************************************
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
  // const auto s = (imag && conj) ? -1 : 1;
  if (imag && conj)
    rhs *= -1;

  // auto rhs = s * hFv + dV_rhs(kappa_x, Fv, conj);
  if (incl_dV)
    rhs += dV_rhs(kappa_x, Fv, conj);
  if (kappa_x == Fv.k && !imag) {
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

  const auto vl = p_hf->get_vlocal(Angular::l_k(kappa_x));
  // The l from X ? or from Fv ?
  return s2 * ExternalField::solveMixedState(kappa_x, Fv, ww, vl, m_alpha,
                                             m_core, rhs, 1.0e-9, Sigma, p_VBr,
                                             m_Hmag);
}

//******************************************************************************
void TDHF::solve_core(const double omega, const int max_its, const bool print) {

  const double converge_targ = 1.0e-8;
  const auto damper = HF::rampedDamp(0.75, 0.25, 1, 20);

  const bool staticQ = std::abs(omega) < 1.0e-10;

  const auto imag = m_h->imaginaryQ();
  const auto has_de = !imag && m_h->parity() == 1;

  // The h*Fc terms don't change, so calculate them just once
  // auto hFcore = m_X;
  std::vector<std::vector<DiracSpinor>> hFcore(m_core.size());
  for (auto ic = 0ul; ic < m_X.size(); ic++) {
    const auto &Fc = m_core[ic];
    hFcore.reserve(m_X[ic].size()); // each h projection
    for (auto beta = 0ul; beta < m_X[ic].size(); beta++) {
      const auto &Xx = m_X[ic][beta];
      hFcore[ic].push_back(m_h->reduced_rhs(Xx.k, Fc));
    }
  }

  if (print) {
    printf("TDHF %s (w=%.3f): ", m_h->name().c_str(), omega);
    std::cout << std::flush;
  }
  auto eps = 0.0;
  double ceiling_eps = 1.0;
  int worse_count = 0;
  double extra_damp = 0.0;
  int it = 0;
  for (; it < max_its; it++) {
    eps = 0.0;
    const auto a_damp = (it == 0) ? 0.0 : damper(it) + extra_damp;

    // eps for solveMixedState - doesn't need to be small!
    // When have de though, equations unstable, so start from scratch
    const auto eps_ms = (it == 0) ? 1.0e-9 : has_de ? 1.0e-9 : 1.0e-3;

    std::vector<double> eps_vec(m_core.size(), 0.0);

    auto tmp_X = m_X; // "temp" allows parallelisation
    auto tmp_Y = m_Y;
#pragma omp parallel for
    for (auto ic = 0ul; ic < m_core.size(); ic++) {
      const auto &Fc = m_core[ic];
      auto eps_c = 0.0;

      // Note: we could, but do not, use solve_dPsi() here.
      // Though it would be cleaner in the code, it is much more efficient
      // To solve manually here, since we don't need to start from scratch

      // delta_en: always same, usually zero; move above!
      const auto de0 = m_h->reducedME(Fc, Fc);
      const auto de1 = dV(Fc, Fc, false);
      const auto de1_dag = dV(Fc, Fc, true);

      for (auto j = 0ul; j < tmp_X[ic].size(); j++) {
        auto &Xx = tmp_X[ic][j];
        const auto &oldX = m_X[ic][j];
        const auto &hFc = hFcore[ic][j];
        auto rhs = hFc + dV_rhs(Xx.k, Fc, false);
        if (Xx.k == Fc.k && !imag)
          rhs -= (de0 + de1) * Fc;
        if (has_de) {
          // Force solveMixedState to start from scratch
          Xx *= 0.0;
        }
        const auto vl = p_hf->get_vlocal(Xx.l()); // to include l-dep QED
        ExternalField::solveMixedState(Xx, Fc, omega, vl, m_alpha, m_core, rhs,
                                       eps_ms, nullptr, p_VBr, m_Hmag);
        Xx = a_damp * oldX + (1.0 - a_damp) * Xx;
        const auto delta = (Xx - oldX) * (Xx - oldX) / (Xx * Xx);
        if (delta > eps_c)
          eps_c = delta;
      }
      if (!staticQ) {
        for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
          auto &Yx = tmp_Y[ic][j];
          const auto &oldY = m_Y[ic][j];
          const auto &hFc = hFcore[ic][j];
          const auto s = imag ? -1 : 1;
          auto rhs = s * hFc + dV_rhs(Yx.k, Fc, true);
          if (Yx.k == Fc.k && !imag)
            rhs -= (de0 + de1_dag) * Fc;
          if (has_de) {
            Yx *= 0.0;
          }
          const auto vl = p_hf->get_vlocal(Yx.l());
          ExternalField::solveMixedState(Yx, Fc, -omega, vl, m_alpha, m_core,
                                         rhs, eps_ms, nullptr, p_VBr, m_Hmag);
          Yx = a_damp * oldY + (1.0 - a_damp) * Yx;
        }
      } else {
        const auto s = imag ? -1 : 1;
        for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
          tmp_Y[ic][j] = s * tmp_X[ic][j];
        }
      }
      eps_vec[ic] = eps_c;
    }
    m_X = tmp_X;
    m_Y = tmp_Y;

    if (it == 0) {
      // On first iteration, store dPsi (for first-order dV)
      m_X0 = m_X;
      m_Y0 = m_Y;
    }

    eps = *std::max_element(cbegin(eps_vec), cend(eps_vec));

    // Work out if converging, or getting worse (early quit)
    if (it > 30 && eps >= ceiling_eps) {
      ++worse_count;
      extra_damp = (it % 2) ? 0.3 : 0.2;
    } else if (eps < ceiling_eps) {
      worse_count = 0;
    }
    ceiling_eps = std::min(eps, ceiling_eps);

    if ((it > 1 && eps < converge_targ) || worse_count > 3)
      break;
  }
  if (print) {
    printf("%2i %.1e\n", it, eps);
    std::cout << std::flush;
  }
  // set last eps (convergance) and frequency (omega)
  m_core_eps = eps;
  m_core_omega = omega;
}

//******************************************************************************
// does it matter if a or b is in the core?
double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm, bool conj,
                const DiracSpinor *const Fexcl, bool incl_dV) const {
  const auto s = conj && m_h->imaginaryQ() ? -1 : 1; // careful. OK?
  return s * Fn * dV_rhs(Fn.k, Fm, conj, Fexcl, incl_dV);
}

double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm) const {
  const auto conj = Fm.en() > Fn.en();
  return dV(Fn, Fm, conj);
}

double TDHF::dV1(const DiracSpinor &Fn, const DiracSpinor &Fm) const {
  const auto conj = Fm.en() > Fn.en();
  return dV(Fn, Fm, conj, nullptr, false);
}

//******************************************************************************
DiracSpinor TDHF::dV_rhs(const int kappa_n, const DiracSpinor &Fa, bool conj,
                         const DiracSpinor *const Fexcl, bool incl_dV) const {

  auto dVFa = DiracSpinor(0, kappa_n, Fa.rgrid);
  dVFa.set_max_pt() = Fa.max_pt();

  const auto ChiType = !conj ? dPsiType::X : dPsiType::Y;
  const auto EtaType = !conj ? dPsiType::Y : dPsiType::X;

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);
  const auto tjn = Angular::twoj_k(kappa_n);

  // nb: faster to not //ize this one
  for (const auto &Fb : m_core) {

    const auto &X_b =
        incl_dV ? get_dPsis(Fb, ChiType) : get_dPsis_0(Fb, ChiType);
    const auto &Y_b =
        incl_dV ? get_dPsis(Fb, EtaType) : get_dPsis_0(Fb, EtaType);

    // only for testing: exclude certain (core) states from dV sum
    if (Fexcl && (*Fexcl) == Fb)
      continue;

    for (auto ibeta = 0ul; ibeta < X_b.size(); ++ibeta) {
      const auto &X_beta = X_b[ibeta];
      const auto &Y_beta = Y_b[ibeta];

      const auto sQ = Angular::neg1pow_2(tjn + X_beta.twoj() + 2 * k + 2);
      if (sQ == 1) { // faster than multiplying!
        dVFa += Coulomb::Wkv_bcd(kappa_n, Fb, Fa, X_beta, k);
        dVFa += Coulomb::Wkv_bcd(kappa_n, Y_beta, Fa, Fb, k);
      } else {
        dVFa -= Coulomb::Wkv_bcd(kappa_n, Fb, Fa, X_beta, k);
        dVFa -= Coulomb::Wkv_bcd(kappa_n, Y_beta, Fa, Fb, k);
      }

      // Breit part:
      if (p_VBr) {
        // Note: Not perfectly symmetric for E1 - some issue??
        dVFa += tkp1 * p_VBr->dVbrD_Fa(kappa_n, k, Fa, Fb, X_beta, Y_beta);
        dVFa += tkp1 * p_VBr->dVbrX_Fa(kappa_n, k, Fa, Fb, X_beta, Y_beta);
      }
    }
  }

  dVFa *= (1.0 / tkp1);

  return dVFa;
}

//******************************************************************************
void TDHF::print(const std::string &ofname) const {
  std::ofstream of(ofname);
  const auto &gr = *((m_core.front()).rgrid);
  for (auto i = 0ul; i < gr.num_points(); ++i) {
    of << gr.r(i) << " ";
    for (auto ic = 0ul; ic < m_core.size(); ic++) {
      const auto &Fc = m_core[ic];
      of << Fc.f(i) << " ";
      for (auto j = 0ul; j < m_X[ic].size(); j++) {
        const auto &Xx = m_X[ic][j];
        const auto &Yx = m_Y[ic][j];
        of << Xx.f(i) << " " << Yx.f(i) << " ";
      }
    }
    of << "\n";
  }
}

} // namespace ExternalField
