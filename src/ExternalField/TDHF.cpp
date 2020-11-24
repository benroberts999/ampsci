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
      m_core(hf->get_core()),
      m_vl(hf->get_vlocal(0)),     // for now, const Vl (same each l)
      m_Hmag(hf->get_Hrad_mag(0)), //(same each l)
      m_alpha(hf->get_alpha()),
      p_VBr(hf->get_Breit())
// XXX Add check for null hf?
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
      m_X[ic].back().pinf = Fc.pinf;
      if (print)
        std::cout << "|" << m_X[ic].back().symbol() << "> + ";
    }
    if (print)
      std::cout << "\n";
  }
  m_Y = m_X;
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
                  const MBPT::CorrelationPotential *const Sigma,
                  StateType st) const {
  std::vector<DiracSpinor> dFvs;
  const auto tjmin = std::max(1, Fv.twoj() - 2 * m_rank);
  const auto tjmax = Fv.twoj() + 2 * m_rank;
  for (int tjbeta = tjmin; tjbeta <= tjmax; tjbeta += 2) {
    const auto kappa = Angular::kappa_twojpi(tjbeta, Fv.parity() * m_pi);
    dFvs.push_back(solve_dPsi(Fv, omega, XorY, kappa, Sigma, st));
  }
  return dFvs;
}
//******************************************************************************
DiracSpinor TDHF::solve_dPsi(const DiracSpinor &Fv, const double omega,
                             dPsiType XorY, const int kappa_x,
                             const MBPT::CorrelationPotential *const Sigma,
                             StateType st) const {
  // Solves (H + Sigma - e - w)X = -(h + dV - de)Psi
  // or     (H + Sigma - e + w)Y = -(h^dag + dV^dag - de)Psi

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;

  const auto imag = m_h->imaginaryQ();

  const auto hPsic = m_h->reduced_rhs(kappa_x, Fv);
  const auto s = (imag && conj) ? -1 : 1;

  auto rhs = s * hPsic + dV_rhs(kappa_x, Fv, conj);
  if (kappa_x == Fv.k && !imag) {
    const auto de = m_h->reducedME(Fv, Fv) + dV(Fv, Fv, conj);
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

  return s2 * ExternalField::solveMixedState(kappa_x, Fv, ww, m_vl, m_alpha,
                                             m_core, rhs, 1.0e-9, Sigma, p_VBr,
                                             m_Hmag);
}

// //******************************************************************************
// void TDHF::solve_core_basis(const std::vector<DiracSpinor> &basis,
//                             const double omega, const int max_its,
//                             const bool print) {
//   //
//   const double converge_targ = 1.0e-8;
//   const auto damper = HF::rampedDamp(0.75, 0.25, 1, 20);
//
//   const bool staticQ = std::abs(omega) < 1.0e-10;
//
//   auto tmp_X = m_X;
//   auto tmp_Y = m_Y;
//
//   auto eps = 0.0;
//   int it = 0;
//   if (print) {
//     printf("TDHFb (w=%.3f): ", omega);
//     std::cout << std::flush;
//   }
//   for (; it < max_its; it++) {
//     eps = 0.0;
//     std::vector<double> eps_vec(m_core.size(), 0.0);
//     const auto a_damp = (it == 0) ? 0.0 : damper(it);
// #pragma omp parallel for
//     for (auto ic = 0ul; ic < m_core.size(); ic++) {
//       double eps_c = 0.0;
//       const auto &Fc = m_core[ic];
//
//       for (auto j = 0ul; j < m_X[ic].size(); j++) {
//         auto &Xx = tmp_X[ic][j];
//         const auto oldX = Xx;
//         Xx *= 0.0;
//
//         for (const auto &Fn : basis) {
//           if (Fn.k != Xx.k || m_h->isZero(Fn.k, Fc.k) || Fc == Fn)
//             continue;
//           const auto hnc = m_h->reducedME(Fn, Fc) + dV(Fn, Fc, false);
//           Xx += (hnc / (Fc.en - Fn.en + omega)) * Fn;
//         }
//
//         const auto delta = (Xx - oldX) * (Xx - oldX) / (Xx * Xx);
//         if (delta > eps_c)
//           eps_c = delta;
//
//         Xx = a_damp * oldX + (1.0 - a_damp) * Xx;
//         if (!staticQ) {
//           auto &Yx = tmp_Y[ic][j];
//           const auto oldY = Yx;
//           Yx *= 0.0;
//           for (const auto &Fn : basis) {
//             if (Fn.k != Yx.k || m_h->isZero(Fn.k, Fc.k) || Fc == Fn)
//               continue;
//             const auto s = m_h->imaginaryQ() ? -1 : 1;
//             const auto hnc =
//                 // s * m_h->reducedME(Fc, Fn) + dV(Fn, Fc, true); //???
//                 s * m_h->reducedME(Fn, Fc) + dV(Fn, Fc, true); //???
//             Yx += (hnc / (Fc.en - Fn.en - omega)) * Fn;
//           }
//           Yx = a_damp * oldY + (1.0 - a_damp) * Yx;
//         } else {
//           tmp_Y[ic][j] = Xx;
//         }
//       }
//       eps_vec[ic] = eps_c;
//     }
//
//     m_X = tmp_X;
//     m_Y = tmp_Y;
//     eps = *std::max_element(cbegin(eps_vec), cend(eps_vec));
//
//     if ((it > 1 && eps < converge_targ))
//       break;
//   }
//   if (print) {
//     printf("%2i %.1e\n", it, eps);
//   }
//   std::cout << std::flush;
//   m_core_eps = eps;
//   m_core_omega = omega;
// }

//******************************************************************************
void TDHF::solve_core(const double omega, const int max_its, const bool print) {

  const double converge_targ = 1.0e-8;
  const auto damper = HF::rampedDamp(0.75, 0.25, 1, 20);

  const bool staticQ = std::abs(omega) < 1.0e-10;

  const auto imag = m_h->imaginaryQ();
  const auto has_de = !imag && m_h->parity() == 1;

  // The h Psi_c terms don't change, so calculate them just once
  auto hPsi = m_X;
  for (auto ic = 0ul; ic < m_core.size(); ic++) {
    const auto &Fc = m_core[ic];
    for (auto j = 0ul; j < hPsi[ic].size(); j++) {
      const auto &Xx = m_X[ic][j];
      hPsi[ic][j] = m_h->reduced_rhs(Xx.k, Fc);
    }
  }

  auto eps = 0.0;
  double ceiling_eps = 1.0;
  int worse_count = 0;
  double extra_damp = 0.0;
  int it = 0;
  if (print) {
    printf("TDHF (w=%.3f): ", omega);
    std::cout << std::flush;
  }
  for (; it < max_its; it++) {
    eps = 0.0;
    const auto a_damp = (it == 0) ? 0.0 : damper(it) + extra_damp;

    // eps for solveMixedState - doesn't need to be small!
    // When have de though, equations unstable, so start from scratch
    const auto eps_ms = (it == 0) ? 1.0e-9 : has_de ? 1.0e-9 : 1.0e-3;

    std::vector<double> eps_vec(m_core.size(), 0.0);

    auto tmp_X = m_X;
    auto tmp_Y = m_Y;
#pragma omp parallel for
    for (auto ic = 0ul; ic < m_core.size(); ic++) {
      const auto &Fc = m_core[ic];
      auto eps_c = 0.0;

      // delta_en: always same, usually zero; move above!
      const auto de0 = m_h->reducedME(Fc, Fc);
      const auto de1 = dV(Fc, Fc, false);
      const auto de1_dag = dV(Fc, Fc, true);

      for (auto j = 0ul; j < tmp_X[ic].size(); j++) {
        auto &Xx = tmp_X[ic][j];
        const auto &oldX = m_X[ic][j];
        const auto &hPsic = hPsi[ic][j];
        auto rhs = hPsic + dV_rhs(Xx.k, Fc, false);
        if (Xx.k == Fc.k && !imag)
          rhs -= (de0 + de1) * Fc;
        if (has_de) {
          // Force solveMixedState to start from scratch
          Xx *= 0.0;
        }
        ExternalField::solveMixedState(Xx, Fc, omega, m_vl, m_alpha, m_core,
                                       rhs, eps_ms, nullptr, p_VBr, m_Hmag);
        Xx = a_damp * oldX + (1.0 - a_damp) * Xx;
        const auto delta = (Xx - oldX) * (Xx - oldX) / (Xx * Xx);
        // use <a|dV|b> instead? But, for which <a|?
        // const auto delta = std::abs(Fc * (Xx - oldX) / (Fc * Xx));
        if (delta > eps_c)
          eps_c = delta;
      }
      if (!staticQ) {
        for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
          auto &Yx = tmp_Y[ic][j];
          const auto &oldY = m_Y[ic][j];
          const auto &hPsic = hPsi[ic][j];
          const auto s = imag ? -1 : 1;
          auto rhs = s * hPsic + dV_rhs(Yx.k, Fc, true);
          if (Yx.k == Fc.k && !imag)
            rhs -= (de0 + de1_dag) * Fc;
          if (has_de) {
            Yx *= 0.0;
          }
          ExternalField::solveMixedState(Yx, Fc, -omega, m_vl, m_alpha, m_core,
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

    eps = *std::max_element(cbegin(eps_vec), cend(eps_vec));

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
  }
  std::cout << std::flush;
  m_core_eps = eps;
  m_core_omega = omega;
}

//******************************************************************************
// does it matter if a or b is in the core?
double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm, bool conj,
                const DiracSpinor *const Fexcl) const {
  auto s = conj && m_h->imaginaryQ() ? -1 : 1; // careful, not always needed
  return s * Fn * dV_rhs(Fn.k, Fm, conj, Fexcl);
}

double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm) const {
  auto conj = Fm.en > Fn.en;
  return dV(Fn, Fm, conj);
}

//******************************************************************************
DiracSpinor TDHF::dV_rhs(const int kappa_n, const DiracSpinor &Fa, bool conj,
                         const DiracSpinor *const Fexcl) const {

  auto dVFa = DiracSpinor(0, kappa_n, Fa.rgrid);
  dVFa.pinf = Fa.pinf;

  const auto ChiType = !conj ? dPsiType::X : dPsiType::Y;
  const auto EtaType = !conj ? dPsiType::Y : dPsiType::X;

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);
  const auto tjn = Angular::twoj_k(kappa_n);

  // nb: faster to not //ize this one
  for (const auto &Fb : m_core) {
    const auto &X_b = get_dPsis(Fb, ChiType);
    const auto &Y_b = get_dPsis(Fb, EtaType);

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

//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************
void TDHF::solve_TDHFcore_matrix(const Wavefunction &wf, const double omega,
                                 const int max_its) {
  // This is just for testing?? Very slow. Should give same as reg method!

  IO::ChronoTimer timer("solve_TDHFcore_matrix");
  const bool staticQ = std::abs(omega) < 1.0e-10;

  const std::size_t nspl = 50;
  const std::size_t kspl = 4;
  const double rmin = 1.0e-4; //?
  const double rmax = 40.0;   //?

  // A := H - (e-w)S
  // b := -h -dV +deS_c

  auto max_ki = 0;
  for (const auto &Px : m_X) {
    for (const auto &x : Px) {
      auto ki = Angular::indexFromKappa(x.k);
      if (ki > max_ki)
        max_ki = ki;
    }
  }
  std::vector<std::vector<DiracSpinor>> basis_kappa;
  basis_kappa.reserve(std::size_t(max_ki + 1));
  for (int ki = 0; ki <= max_ki; ki++) {
    auto k = Angular::kappaFromIndex(ki);
    basis_kappa.push_back(SplineBasis::form_spline_basis(
        k, nspl, kspl, rmin, rmax, m_core.front().rgrid, m_alpha));
  }

  const auto imag = m_h->imaginaryQ();

  const double converge_targ = 1.0e-4;
  const auto damper = HF::rampedDamp(0.5, 0.25, 1, 10);

  auto eps = 0.0;
  for (int it = 0; it < max_its; it++) {
    IO::ChronoTimer timer2("solve_core: iterations");
    eps = 0.0;
    const auto a_damp = (it == 0) ? 0.0 : damper(it);

    auto tmp_X = m_X;
    auto tmp_Y = m_Y;
#pragma omp parallel for
    for (std::size_t ic = 0; ic < m_core.size(); ic++) {
      const auto &Fc = m_core[ic];
      auto &dPsiX = tmp_X[ic];
      auto &dPsiY = tmp_Y[ic];
      const auto de0 = m_h->reducedME(Fc, Fc);
      const auto de1 = dV(Fc, Fc, false);
      const auto de1_dag = dV(Fc, Fc, true);
      for (auto ibeta = 0ul; ibeta < dPsiX.size(); ++ibeta) {
        auto &Xx = dPsiX[ibeta];
        auto &Yx = dPsiY[ibeta];
        const auto ki = std::size_t(Angular::indexFromKappa(Xx.k));
        const auto &basis = basis_kappa[ki];

        LinAlg::Vector bi_X(basis.size());
        LinAlg::Vector bi_Y(basis.size());
        for (auto i = 0ul; i < basis.size(); ++i) {
          const auto &xi = basis[i];
          // fill LHS vector, b
          const auto hi = m_h->reducedME(xi, Fc);
          const auto hidag = m_h->reducedME(Fc, xi); //??
          const auto dVic = dV(xi, Fc, false);
          const auto dV_dag = dV(xi, Fc, true);
          const auto s = imag ? -1 : 1;
          const auto Sic = (xi.k == Fc.k && !imag) ? (xi * Fc) : 0.0;
          const auto deS = (de0 + de1) * Sic;
          const auto deS_dag = (de0 + de1_dag) * Sic;
          bi_X[int(i)] = -s * hi - dVic + deS; // why s here? check above??
          // bi_Y[i] = -s * (s * hi + dV_dag) + deS_dag;
          bi_Y[int(i)] = -s * hidag - s * dV_dag + deS_dag; //???
        }
        const auto [Hij, Sij] = SplineBasis::fill_Hamiltonian_matrix(basis, wf);

        auto Aij_X = Hij - (Fc.en + omega) * Sij;
        auto Aij_Y = Hij - (Fc.en - omega) * Sij;

        auto s = imag ? -1 : 1;
        const auto c_X = LinAlg::solve_Axeqb(Aij_X, bi_X);
        const auto c_Y = staticQ ? s * c_X : LinAlg::solve_Axeqb(Aij_Y, bi_Y);

        Xx.scale(a_damp);
        Yx.scale(a_damp);
        for (auto i = 0ul; i < basis.size(); ++i) {
          Xx += (1.0 - a_damp) * c_X[int(i)] * basis[i];
          Yx += (1.0 - a_damp) * c_Y[int(i)] * basis[i];
        }
      }
    }

    for (std::size_t ic = 0; ic < m_core.size(); ic++) {
      for (auto ibeta = 0ul; ibeta < tmp_X[ic].size(); ++ibeta) {
        const auto &dF = tmp_X[ic][ibeta];
        const auto &dF0 = m_X[ic][ibeta];
        const auto eps_c = (dF - dF0) * (dF - dF0) / (dF * dF);
        if (eps_c > eps)
          eps = eps_c;
      }
    }

    m_Y = tmp_Y;
    m_X = tmp_X;
    printf("TDHF [matrix] (w=%.3f): %2i  %.1e\n", omega, it + 1, eps);
    if (it > 0 && eps < converge_targ)
      break;
  }
}

//******************************************************************************
void TDHF::print(const std::string &ofname) const {
  std::ofstream of(ofname);
  const auto &gr = *((m_core.front()).rgrid);
  for (auto i = 0ul; i < gr.num_points; ++i) {
    of << gr.r[i] << " ";
    for (auto ic = 0ul; ic < m_core.size(); ic++) {
      const auto &Fc = m_core[ic];
      of << Fc.f[i] << " ";
      for (auto j = 0ul; j < m_X[ic].size(); j++) {
        const auto &Xx = m_X[ic][j];
        const auto &Yx = m_Y[ic][j];
        of << Xx.f[i] << " " << Yx.f[i] << " ";
      }
    }
    of << "\n";
  }
}

//******************************************************************************

//******************************************************************************

//******************************************************************************
//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************
//******************************************************************************

//******************************************************************************
double TDHF::dX_nm_bbe_rhs(const DiracSpinor &Fn, const DiracSpinor &Fm,
                           const DiracSpinor &Fb,
                           const DiracSpinor &X_beta) const {

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  double dX_nm_bbe = 0.0;

  const auto tjn = Fn.twoj();
  const auto tjm = Fm.twoj();
  const auto Ckala = Angular::Ck_kk(k, Fn.k, Fm.k);

  const auto tjb = Fb.twoj();
  const auto tjbeta = X_beta.twoj();
  const auto Ckbeb = Angular::Ck_kk(k, X_beta.k, Fb.k);

  if (Ckala != 0 && Ckbeb != 0) {
    const auto Rkabcd = Coulomb::Rk_abcd(Fn, Fb, Fm, X_beta, k);
    dX_nm_bbe += (Ckala * Ckbeb / tkp1) * Rkabcd;
  }

  auto s = Angular::evenQ_2(tjn + tjbeta + 2) ? 1 : -1;

  // exchange part (X):
  const auto l_min_X =
      std::max(std::abs(tjn - tjbeta), std::abs(tjm - tjb)) / 2;
  const auto l_max_X = std::min((tjn + tjbeta), (tjm + tjb)) / 2;

  for (int l = l_min_X; l <= l_max_X; ++l) {
    const auto sixj = Angular::sixj_2(tjm, tjn, 2 * k, tjbeta, tjb, 2 * l);
    if (sixj == 0)
      continue;
    const auto m1kpl = Angular::evenQ(k + l) ? 1 : -1;
    const auto Ckba = Angular::Ck_kk(l, Fm.k, Fb.k);
    const auto Ckalbe = Angular::Ck_kk(l, Fn.k, X_beta.k);
    if (Ckba == 0 || Ckalbe == 0)
      continue;
    const auto Rk = Coulomb::Rk_abcd(Fn, Fm, X_beta, Fb, l);
    dX_nm_bbe += (s * m1kpl * Ckba * Ckalbe * sixj) * Rk;
  }

  return dX_nm_bbe;
}

//******************************************************************************
double TDHF::dY_nm_bbe_rhs(const DiracSpinor &Fn, const DiracSpinor &Fm,
                           const DiracSpinor &Fb,
                           const DiracSpinor &Y_beta) const {

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  double dY_nm_bbe = 0.0;

  const auto tjn = Fn.twoj();
  const auto tjm = Fm.twoj();
  const auto Ckala = Angular::Ck_kk(k, Fn.k, Fm.k);

  const auto tjb = Fb.twoj();
  const auto tjbeta = Y_beta.twoj();
  const auto Ckbeb = Angular::Ck_kk(k, Y_beta.k, Fb.k);

  if (Ckala != 0 && Ckbeb != 0) {
    const auto Rkabcd = Coulomb::Rk_abcd(Fn, Fb, Fm, Y_beta, k);
    dY_nm_bbe += (Ckala * Ckbeb / tkp1) * Rkabcd;
  }

  auto s = Angular::evenQ_2(tjn + tjbeta + 2) ? 1 : -1;

  // exchange part (Y):
  auto l_min_Y = std::max(std::abs(tjn - tjb), std::abs(tjm - tjbeta)) / 2;
  auto l_max_Y = std::min((tjn + tjb), (tjm + tjbeta)) / 2;

  // l_min_Y = 0;
  // l_max_Y = 20;

  for (int l = l_min_Y; l <= l_max_Y; ++l) {
    const auto sixj = Angular::sixj_2(tjm, tjn, 2 * k, tjb, tjbeta, 2 * l);
    if (sixj == 0)
      continue;
    const auto m1kpl = Angular::evenQ(k + l) ? 1 : -1;
    const auto Ckbea = Angular::Ck_kk(l, Fm.k, Y_beta.k);
    const auto Ckbal = Angular::Ck_kk(l, Fn.k, Fb.k);
    if (Ckbea == 0 || Ckbal == 0)
      continue;
    const auto Rk = Coulomb::Rk_abcd(Fn, Fm, Fb, Y_beta, l);
    dY_nm_bbe += (s * m1kpl * Ckbea * Ckbal * sixj) * Rk;
  }

  return dY_nm_bbe;
}

} // namespace ExternalField
